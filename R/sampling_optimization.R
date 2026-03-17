#' Sampling optimization for a target species
#'
#' Runs bootstrap-based sampling optimization for a single target species using
#' COST catch sampling (`CS`) and landing (`CL`) objects. The function evaluates
#' the coefficient of variation (CV) of raised length distributions across
#' strata and derives an "optimal" sample size range.
#'
#' @param CS A COST catch sampling object.
#' @param CL A COST landing object.
#' @param metier_shortcode A data.frame with metier recoding information.
#' @param past_situations A data.frame produced by [data_investigation()].
#' @param space Space stratification. One of `"GSA"`, `"Country"`, `"ALL"`.
#' @param time Time stratification. One of `"Y"`, `"Q"`, `"S"`.
#' @param technical Technical stratification. One of `"lev4"`, `"lev6"`, `"NONE"`.
#' @param target_SPECIES Character string with the target species scientific name.
#' @param species_code Optional short code used in output file names.
#' @param nIter Number of bootstrap iterations.
#' @param threshold_for_recycling_rate Threshold used in recycling-rate plots.
#' @param output_dir Directory where outputs will be saved.
#' @param table_min_accepted_sample_size table containing for each stratum the minimum number of trips to be considered
#' @return A list containing bootstrap results, selected sample size ranges,
#'   available strata, and output file paths.
#' @export
#' @import COSTcore
#' @import COSTdbe
#' @import COSTeda
#' @import ggplot2
sampling_optimization <- function(CS,
                                  CL,
                                  metier_shortcode,
                                  past_situations,
                                  space = c("GSA", "Country", "ALL"),
                                  time = c("Y", "Q", "S"),
                                  technical = c("lev4", "lev6", "NONE"),
                                  target_SPECIES,
                                  species_code = NULL,
                                  nIter = 50,
                                  threshold_for_recycling_rate = 0.25,
                                  output_dir = tempdir(),
                                  table_min_accepted_sample_size) {

    if (FALSE){
        CS = costCS_from_CS
        CL = costCL_from_CL
        metier_shortcode = metier_shortcode
        past_situations = past_situations
        space = "GSA"
        time = "Q"
        technical = "lev4"
        target_SPECIES = target_species
        nIter = 20      # basso per test veloce
        output_dir = out_dir
    }

    ## -------------------------------------------------------------------------
    ## Checks
    ## -------------------------------------------------------------------------
    if (!inherits(CS, c("csData", "csDataVal"))) {
        stop("`CS` must be a csData or csDataVal object.")
    }

    if (!inherits(CL, c("clData", "clDataVal", "csData", "csDataVal"))) {
        stop("`CL` must be a valid COST landing object.")
    }

    if (!is.character(target_SPECIES) || length(target_SPECIES) != 1L || !nzchar(target_SPECIES)) {
        stop("`target_SPECIES` must be a single non-empty character string.")
    }


    needed_metier_cols <- c("Fishing_activity_category_level6", "lev4", "short_code")
    if (!all(needed_metier_cols %in% names(metier_shortcode))) {
        stop("`metier_shortcode` must contain: ", paste(needed_metier_cols, collapse = ", "))
    }

    if (!is.data.frame(past_situations)) {
        stop("`past_situations` must be a data.frame.")
    }

    if (is.null(species_code) || !nzchar(species_code)) {
        species_code <- gsub(" ", "_", target_SPECIES)
    }

    res_dir <- file.path(output_dir, paste0("Optimization_", species_code))
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

    ## -------------------------------------------------------------------------
    ## Harmonise COST objects
    ## -------------------------------------------------------------------------
    CS@hh$date <- to_date(CS@hh$date)

    CS@hh$foVal <- "V"
    CS@sl$landCat <- "HUC"
    CS@hl$landCat <- "HUC"
    CS@sl$subSampCat <- ""
    CS@hl$subSampCat <- ""

    CS@sl$catchCat <- toupper(CS@sl$catchCat)
    CS@hl$catchCat <- toupper(CS@hl$catchCat)

    if ("catchCat" %in% names(CS@ca)) {
        CS@ca$catchCat <- toupper(CS@ca$catchCat)
    }

    ## In the original Rmd this was used to avoid issues in bootstrap steps
    CS@ca <- COSTcore::csData()@ca

    recode_metier <- function(df, col6 = "foCatEu6", col5 = "foCatEu5") {
        idx <- match(as.character(df[[col6]]),
                     as.character(metier_shortcode$Fishing_activity_category_level6))

        if (!col5 %in% names(df)) {
            df[[col5]] <- NA_character_
        }

        found <- !is.na(idx)
        df[[col5]][found] <- as.character(metier_shortcode$lev4[idx[found]])
        df[[col6]][found] <- as.character(metier_shortcode$short_code[idx[found]])

        df
    }

    CL@cl <- recode_metier(CL@cl, col6 = "foCatEu6", col5 = "foCatEu5")
    CS@hh <- recode_metier(CS@hh, col6 = "foCatEu6", col5 = "foCatEu5")

    ## -------------------------------------------------------------------------
    ## Species subset
    ## -------------------------------------------------------------------------
    ## -------------------------------------------------------------------------
    ## Species subsets
    ## -------------------------------------------------------------------------
    CL_boot <- CL
    if ("taxon" %in% names(CL_boot@cl)) {
        CL_boot@cl <- CL_boot@cl[CL_boot@cl$taxon %in% target_SPECIES, , drop = FALSE]
    }

    CS_boot <- CS
    if ("spp" %in% names(CS_boot@hl)) {
        CS_boot@hl <- CS_boot@hl[CS_boot@hl$spp %in% target_SPECIES, , drop = FALSE]
    }
    if ("spp" %in% names(CS_boot@sl)) {
        CS_boot@sl <- CS_boot@sl[CS_boot@sl$spp %in% target_SPECIES, , drop = FALSE]
    }

    if (nrow(CS_boot@hl) == 0L || nrow(CL_boot@cl) == 0L) {
        stop("No data available for target species after subsetting.")
    }

    ## -------------------------------------------------------------------------
    ## Force bootstrap year as in original script
    ## -------------------------------------------------------------------------
    CS_boot@tr$year <- 2999
    CS_boot@hh$year <- 2999
    CS_boot@sl$year <- 2999
    CS_boot@hl$year <- 2999

    availSamp <- merge(CS_boot@hh, CS_boot@hl)
    availLand <- CL_boot@cl

    availSamp$SampQuarter <- get_quarter(availSamp$date)
    availSamp$SampSemester <- get_semester(availSamp$date)

    availLand$semester <- ifelse(availLand$quarter %in% c(3, 4), 2L, 1L)

    ## -------------------------------------------------------------------------
    ## Stratification settings
    ## -------------------------------------------------------------------------
    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    stratification_table <- data.frame(
        space = space,
        time = time,
        technical = technical,
        stringsAsFactors = FALSE
    )

    col_ele_var1 <- as.character(stratification_table[1, "space"])
    col_ele_var2 <- as.character(stratification_table[1, "time"])
    col_ele_var3 <- as.character(stratification_table[1, "technical"])

    Var1 <- switch(
        col_ele_var1,
        GSA = availSamp$area,
        Country = availSamp$vslFlgCtry,
        ALL = rep("ALL", nrow(availSamp)),
        stop("Unsupported value in `space`: ", col_ele_var1)
    )

    Var1_lan <- switch(
        col_ele_var1,
        GSA = availLand$area,
        Country = availLand$vslFlgCtry,
        ALL = rep("ALL", nrow(availLand)),
        stop("Unsupported value in `space`: ", col_ele_var1)
    )

    Var2 <- switch(
        col_ele_var2,
        Q = availSamp$SampQuarter,
        S = availSamp$SampSemester,
        Y = rep(1, nrow(availSamp)),
        stop("Unsupported value in `time`: ", col_ele_var2)
    )

    Var2_lan <- switch(
        col_ele_var2,
        Q = availLand$quarter,
        S = availLand$semester,
        Y = rep(1, nrow(availLand)),
        stop("Unsupported value in `time`: ", col_ele_var2)
    )

    Var3 <- switch(
        col_ele_var3,
        lev6 = availSamp$foCatEu6,
        lev4 = availSamp$foCatEu5,
        NONE = rep("1", nrow(availSamp)),
        stop("Unsupported value in `technical`: ", col_ele_var3)
    )

    Var3_lan <- switch(
        col_ele_var3,
        lev6 = availLand$foCatEu6,
        lev4 = availLand$foCatEu5,
        NONE = rep("1", nrow(availLand)),
        stop("Unsupported value in `technical`: ", col_ele_var3)
    )

    availStratSamp <- unique(data.frame(Var1 = Var1, Var2 = Var2, Var3 = Var3))
    availStratLand <- unique(data.frame(Var1 = Var1_lan, Var2 = Var2_lan, Var3 = Var3_lan))
    tabStrat <- merge(availStratSamp, availStratLand)

    if (nrow(tabStrat) == 0L) {
        stop("No common strata found between sampling and landing data.")
    }

    ## -------------------------------------------------------------------------
    ## Mean number of measured individuals by trip within stratum
    ## NOTE: assumes meanNmeasByTripByStr() exists in the package
    ## -------------------------------------------------------------------------
    availStrat <- Reduce(
        "rbind",
        lapply(seq_len(nrow(tabStrat)), function(i) {
            meanNmeasByTripByStr(
                i,
                tabStrat,
                CS_boot,
                space = space,
                time = time,
                technical = technical
            )
        })
    )

    if (nrow(availStrat) == 0L) {
        stop("No available strata after computing mean numbers by trip.")
    }

    availStrat$maxSampSize <- round(5000 / availStrat$meanNmeasByTrip)
    availStrat$maxSampSize[availStrat$maxSampSize > 1000] <- 1000
    availStrat$maxSampSize[availStrat$maxSampSize < 100] <- 100

    ## -------------------------------------------------------------------------
    ## Bootstrap CV estimation
    ## NOTE: assumes cvFunStr() fills/returns the needed output
    ## -------------------------------------------------------------------------
    LengthStr_dataframe <- data.frame(
        iteration = integer(0),
        Var1 = character(0),
        Var2 = character(0),
        Var3 = character(0),
        length = numeric(0),
        value = numeric(0)
    )

    resOptCV_list <- vector("list", nrow(availStrat))

    LengthStr_dataframe_all <- data.frame(
        iteration = integer(0),
        Var1 = character(0),
        Var2 = character(0),
        Var3 = character(0),
        length = numeric(0),
        value = numeric(0)
    )

    for (nr in seq_len(nrow(availStrat))) {

        resOptCV <- try(
            cvFunStr(
                caseInd = nr,
                stratification = COSTcore::strIni(),
                availStrat = availStrat,
                availSamp = availSamp,
                CSobject = CS_boot,
                CLobjectList = CL_boot,
                sppName = target_SPECIES,
                nIter = nIter,
                space = space,
                time = time,
                technical = technical,
                lengthStr_dataframe = LengthStr_dataframe
            ),
            silent = TRUE
        )

        if (inherits(resOptCV, "try-error")) {
            next
        }

        ## prendi solo il data frame delle iterazioni
        this_iter <- resOptCV$res_all_iter

        if (!is.null(this_iter) && is.data.frame(this_iter) && nrow(this_iter) > 0L) {
            this_iter$Var1 <- availStrat$Var1[nr]
            this_iter$Var2 <- availStrat$Var2[nr]
            this_iter$Var3 <- availStrat$Var3[nr]
            resOptCV_list[[nr]] <- this_iter
        }

        ## accumula anche le length distributions, se presenti
        this_len <- resOptCV$LengthStr_dataframe
        if (!is.null(this_len) && is.data.frame(this_len) && nrow(this_len) > 0L) {
            LengthStr_dataframe_all <- rbind(LengthStr_dataframe_all, this_len)
        }
    }

    ## tieni solo gli elementi validi
    resOptCV_list <- Filter(Negate(is.null), resOptCV_list)

    resOptCV_final <- data.table::rbindlist(resOptCV_list, fill = TRUE)
    if (nrow(resOptCV_final) == 0L) {
        stop("Bootstrap CV estimation failed for all strata.")
    }

    cv_boot_file <- file.path(res_dir, paste0(species_code, "_CV_from_bootstrap.csv"))
    write.table(resOptCV_final, cv_boot_file, sep = ";", row.names = FALSE, na = "")

    raised_length_file <- file.path(res_dir, paste0(species_code, "_Raised_length_distribution_from_bootstrap.csv"))
    write.table(LengthStr_dataframe, raised_length_file, sep = ";", row.names = FALSE, na = "")

    ## -------------------------------------------------------------------------
    ## Bootstrap CV plots
    ## -------------------------------------------------------------------------
    bootstrap_cv_plot_files <- character(0)

    for (i in seq_len(nrow(availStrat))) {
        this_str_CV <- resOptCV_final[
            as.character(resOptCV_final$Var1) == as.character(availStrat$Var1[i]) &
                as.character(resOptCV_final$Var2) == as.character(availStrat$Var2[i]) &
                as.character(resOptCV_final$Var3) == as.character(availStrat$Var3[i]),
            ,
            drop = FALSE
        ]

        if (nrow(this_str_CV) == 0L) next

        plot_title <- paste(
            target_SPECIES,
            availStrat$Var1[i], availStrat$Var2[i], availStrat$Var3[i],
            "- Bootstrap CVs"
        )

        gg <- ggplot2::ggplot(this_str_CV, ggplot2::aes(x = sampSize, y = cv)) +
            ggplot2::geom_point() +
            ggplot2::geom_smooth() +
            ggplot2::labs(
                title = plot_title,
                x = "# trips",
                y = "coefficient of variation"
            )

        file_name <- file.path(
            res_dir,
            paste0(
                species_code, "_Bootstrap_CVs_",
                availStrat$Var1[i], "_", availStrat$Var2[i], "_", availStrat$Var3[i], ".png"
            )
        )

        ggplot2::ggsave(file_name, gg, width = 8, height = 6)
        bootstrap_cv_plot_files <- c(bootstrap_cv_plot_files, file_name)
    }

    ## -------------------------------------------------------------------------
    ## Optimal sample size range
    ## NOTE: assumes selectSampRange2() exists in the package
    ## -------------------------------------------------------------------------
    sampSizeRange_list <- vector("list", nrow(availStrat))

    for (nr_str in seq_len(nrow(availStrat))) {
        sampSizeRange_list[[nr_str]] <- selectSampRange2(nr_str, availStrat, resOptCV_final,table_min_accepted_sample_size,output_dir,species_code)
    }

    sampSizeRange <- data.table::rbindlist(sampSizeRange_list, fill = TRUE)

    optimized_sample_size_file <- file.path(res_dir, paste0(species_code, "_Optimized_sample_size.csv"))
    write.table(sampSizeRange, optimized_sample_size_file, sep = ";", row.names = FALSE, na = "")

    ## -------------------------------------------------------------------------
    ## Recycling rate plots
    ## -------------------------------------------------------------------------
    recycling_rate_plot_files <- character(0)

    add_solution_lines <- function(p, blue_lines) {
        if (length(blue_lines) > 0L) {
            for (bl in blue_lines) {
                p <- p + ggplot2::geom_vline(xintercept = bl, colour = "blue", linewidth = 0.2)
            }
        }
        p
    }

    for (i in seq_len(nrow(availStrat))) {
        this_str_CV <- resOptCV_final[
            as.character(resOptCV_final$Var1) == as.character(availStrat$Var1[i]) &
                as.character(resOptCV_final$Var2) == as.character(availStrat$Var2[i]) &
                as.character(resOptCV_final$Var3) == as.character(availStrat$Var3[i]),
            ,
            drop = FALSE
        ]

        if (nrow(this_str_CV) == 0L || !"recyclingRate" %in% names(this_str_CV)) next

        sampSizeRange_this_str <- sampSizeRange[
            as.character(sampSizeRange$Var1) == as.character(availStrat$Var1[i]) &
                as.character(sampSizeRange$Var2) == as.character(availStrat$Var2[i]) &
                as.character(sampSizeRange$Var3) == as.character(availStrat$Var3[i]),
            ,
            drop = FALSE
        ]

        blue_lines <- suppressWarnings(
            as.numeric(sampSizeRange_this_str$solutions[!is.na(as.numeric(sampSizeRange_this_str$solutions))])
        )

        plot_title <- paste(
            target_SPECIES,
            availStrat$Var1[i], availStrat$Var2[i], availStrat$Var3[i],
            "- Recycling Rate bootstrap"
        )

        gg <- ggplot2::ggplot(this_str_CV) +
            ggplot2::geom_point(ggplot2::aes(x = sampSize, y = recyclingRate), size = 0.4) +
            ggplot2::geom_smooth(
                ggplot2::aes(x = sampSize, y = recyclingRate),
                fill = NA,
                linewidth = 0.2,
                colour = "black"
            ) +
            ggplot2::geom_hline(
                yintercept = threshold_for_recycling_rate,
                colour = "red",
                linewidth = 0.2
            ) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 7),
                axis.text.y = ggplot2::element_text(size = 7),
                legend.position = "top"
            ) +
            ggplot2::labs(
                title = plot_title,
                x = "Number of trips",
                y = "Data Recycling Rate"
            )

        gg <- add_solution_lines(gg, blue_lines)

        file_name <- file.path(
            res_dir,
            paste0(
                species_code, "_RecyclingRate_vs_trips_",
                availStrat$Var1[i], "_", availStrat$Var2[i], "_", availStrat$Var3[i], ".png"
            )
        )

        ggplot2::ggsave(file_name, gg, width = 10, height = 8)
        recycling_rate_plot_files <- c(recycling_rate_plot_files, file_name)
    }

    ## -------------------------------------------------------------------------
    ## Comparison with past situation
    ## -------------------------------------------------------------------------
    comparison_plot_files <- character(0)

    if (nrow(past_situations) > 0L) {
        pastData <- merge(past_situations, availStrat, by = c("Var1", "Var2", "Var3"))
        resOptCV_final$cv[is.na(resOptCV_final$cv)] <- 0

        for (i in seq_len(nrow(availStrat))) {
            resOptCV_final_str <- resOptCV_final[
                as.character(resOptCV_final$Var1) == as.character(availStrat$Var1[i]) &
                    as.character(resOptCV_final$Var2) == as.character(availStrat$Var2[i]) &
                    as.character(resOptCV_final$Var3) == as.character(availStrat$Var3[i]),
                ,
                drop = FALSE
            ]

            pastData_str <- pastData[
                as.character(pastData$Var1) == as.character(availStrat$Var1[i]) &
                    as.character(pastData$Var2) == as.character(availStrat$Var2[i]) &
                    as.character(pastData$Var3) == as.character(availStrat$Var3[i]),
                ,
                drop = FALSE
            ]

            if (nrow(resOptCV_final_str) == 0L) next

            sampSizeRange_str <- sampSizeRange[
                as.character(sampSizeRange$Var1) == as.character(availStrat$Var1[i]) &
                    as.character(sampSizeRange$Var2) == as.character(availStrat$Var2[i]) &
                    as.character(sampSizeRange$Var3) == as.character(availStrat$Var3[i]),
                ,
                drop = FALSE
            ]

            blue_lines <- suppressWarnings(
                as.numeric(sampSizeRange_str$solutions[!is.na(as.numeric(sampSizeRange_str$solutions))])
            )

            plot_title <- paste(
                target_SPECIES,
                availStrat$Var1[i], availStrat$Var2[i], availStrat$Var3[i],
                "- CV vs no. of trips"
            )

            gg <- ggplot2::ggplot(resOptCV_final_str) +
                ggplot2::geom_point(
                    data = pastData_str,
                    mapping = ggplot2::aes(x = sampSize, y = cv, colour = as.factor(year)),
                    size = 1.8
                ) +
                ggplot2::geom_point(
                    mapping = ggplot2::aes(x = sampSize, y = cv),
                    size = 0.6
                ) +
                ggplot2::geom_smooth(
                    ggplot2::aes(x = sampSize, y = cv),
                    fill = NA,
                    linewidth = 0.2,
                    colour = "black"
                ) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 7),
                    axis.text.y = ggplot2::element_text(size = 7),
                    legend.position = "top"
                ) +
                ggplot2::labs(
                    title = plot_title,
                    x = "number of trips",
                    y = "Coefficient of Variation"
                )

            gg <- add_solution_lines(gg, blue_lines)

            file_name <- file.path(
                res_dir,
                paste0(
                    species_code, "_CV_vs_trips_",
                    availStrat$Var1[i], "_", availStrat$Var2[i], "_", availStrat$Var3[i], ".png"
                )
            )

            ggplot2::ggsave(file_name, gg, width = 10, height = 8)
            comparison_plot_files <- c(comparison_plot_files, file_name)
        }
    }

    ## -------------------------------------------------------------------------
    ## Return
    ## -------------------------------------------------------------------------
    invisible(list(
        species = target_SPECIES,
        species_code = species_code,
        output_dir = res_dir,
        stratification = list(
            space = col_ele_var1,
            time = col_ele_var2,
            technical = col_ele_var3
        ),
        available_strata = availStrat,
        bootstrap_cv = resOptCV_final,
        optimized_sample_size = sampSizeRange,
        files = list(
            bootstrap_cv_csv = cv_boot_file,
            raised_length_csv = raised_length_file,
            optimized_sample_size_csv = optimized_sample_size_file,
            bootstrap_cv_plots = bootstrap_cv_plot_files,
            recycling_rate_plots = recycling_rate_plot_files,
            comparison_plots = comparison_plot_files
        )
    ))
}
