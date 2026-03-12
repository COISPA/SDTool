#' Data investigation on COST tables (CS and CL)
#'
#' Performs exploratory analyses on COST catch sampling (`CS`) and landing (`CL`)
#' objects for a single target species.
#'
#' @param CS A COST catch sampling object.
#' @param CL A COST landing object.
#' @param metier_shortcode A data.frame with metier recoding information.
#' @param space Space stratification. One of `"GSA"`, `"Country"`, `"ALL"`.
#' @param time Time stratification. One of `"Y"`, `"Q"`, `"S"`.
#' @param technical Technical stratification. One of `"lev4"`, `"lev6"`, `"NONE"`.
#' @param target_SPECIES Character string with the target species scientific name.
#' @param output_dir Directory where outputs will be saved.
#'
#' @return A list containing summary tables and output file paths.
#' @export
#' @import COSTcore
#' @import COSTdbe
#' @import COSTeda
#' @import ggplot2
data_investigation <- function(CS,
                               CL,
                               metier_shortcode,
                               space = c("GSA", "Country", "ALL"),
                               time = c("Y", "Q", "S"),
                               technical = c("lev4", "lev6", "NONE"),
                               target_SPECIES,
                               output_dir = tempdir()) {

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    if (!inherits(CS, c("csData", "csDataVal"))) {
        stop("`CS` must be a csData or csDataVal object.")
    }
    if (!inherits(CL, c("clData", "clDataVal", "csData", "csDataVal"))) {
        stop("`CL` must be a valid COST landing object.")
    }
    if (!is.character(target_SPECIES) || length(target_SPECIES) != 1L || !nzchar(target_SPECIES)) {
        stop("`target_SPECIES` must be a single non-empty character string.")
    }

    needed_cols <- c("Fishing_activity_category_level6", "lev4", "short_code")
    if (!all(needed_cols %in% names(metier_shortcode))) {
        stop("`metier_shortcode` must contain: ", paste(needed_cols, collapse = ", "))
    }

    dir.create("DataInvestigation", recursive = TRUE, showWarnings = FALSE)

    res_dir <- paste(output_dir,"\\DataInvestigation",sep="")

    ## ---------------------------------------------------------------------------
    ## Harmonise COST objects
    ## ---------------------------------------------------------------------------
    CS@hh$date <- to_date(CS@hh$date)

    CS@hh$foVal <- "V"
    CS@sl$landCat <- "HUC"
    CS@hl$landCat <- "HUC"
    CS@sl$subSampCat <- ""
    CS@hl$subSampCat <- ""

    CS@sl$catchCat <- toupper(CS@sl$catchCat)
    CS@hl$catchCat <- toupper(CS@hl$catchCat)
    CS@ca$catchCat <- toupper(CS@ca$catchCat)

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

    years <- sort(unique(CS@tr$year))
    gsas <- unique(CS@hh$area)

    sppName <- target_SPECIES
    sppShortCode <- gsub(" ", "_", target_SPECIES)

    stratification <- COSTcore::strIni()

    ## ---------------------------------------------------------------------------
    ## Species subsets
    ## ---------------------------------------------------------------------------
    cl_spe <- CL
    if ("taxon" %in% names(cl_spe@cl)) {
        cl_spe@cl <- cl_spe@cl[cl_spe@cl$taxon %in% sppName, , drop = FALSE]
    }

    cs_spe <- CS
    if ("spp" %in% names(cs_spe@hl)) {
        cs_spe@hl <- cs_spe@hl[cs_spe@hl$spp %in% sppName, , drop = FALSE]
    }
    if ("spp" %in% names(cs_spe@sl)) {
        cs_spe@sl <- cs_spe@sl[cs_spe@sl$spp %in% sppName, , drop = FALSE]
    }

    ## ---------------------------------------------------------------------------
    ## 1. Past situation
    ## ---------------------------------------------------------------------------
    past_list <- vector("list", length(years))
    names(past_list) <- as.character(years)

    for (ye in years) {
        print(paste("Investigate data for", ye, sppName, "in progress..."), quote = FALSE)

        cs_ye <- cs_spe
        cs_ye@hh <- cs_ye@hh[cs_ye@hh$year %in% ye, , drop = FALSE]
        cs_ye@hl <- cs_ye@hl[cs_ye@hl$year %in% ye, , drop = FALSE]
        cs_ye@sl <- cs_ye@sl[cs_ye@sl$year %in% ye, , drop = FALSE]

        cl_ye <- cl_spe
        cl_ye@cl <- cl_ye@cl[cl_ye@cl$year %in% ye, , drop = FALSE]

        if (nrow(cs_ye@hh) == 0L || nrow(cs_ye@hl) == 0L || nrow(cl_ye@cl) == 0L) {
            next
        }

        hsl <- merge(cs_ye@hh, cs_ye@hl)
        hsl$SampQuarter <- get_quarter(hsl$date)
        hsl$SampSemester <- get_semester(hsl$date)

        ## restore foCatEu5 in merged table if needed
        if (!"foCatEu5" %in% names(hsl) && "foCatEu6" %in% names(hsl)) {
            idx_lev4 <- match(as.character(hsl$foCatEu6), as.character(metier_shortcode$short_code))
            hsl$foCatEu5 <- as.character(metier_shortcode$lev4[idx_lev4])
        }

        cl_temp <- cl_ye@cl
        cl_temp$semester <- ifelse(cl_temp$quarter %in% c(3, 4), 2L, 1L)

        hsl <- setVar123(
            x = hsl,
            space = space,
            time = time,
            technical = technical,
            is_landings = FALSE
        )

        cl_temp <- setVar123(
            x = cl_temp,
            space = space,
            time = time,
            technical = technical,
            is_landings = TRUE
        )

        hsl <- hsl[!is.na(hsl$Var1) & !is.na(hsl$Var2) & !is.na(hsl$Var3), , drop = FALSE]
        cl_temp <- cl_temp[!is.na(cl_temp$Var1) & !is.na(cl_temp$Var2) & !is.na(cl_temp$Var3), , drop = FALSE]

        if (nrow(hsl) == 0L || nrow(cl_temp) == 0L) {
            next
        }

        tabsamp <- aggregate(lenNum ~ Var1 + Var2 + Var3, data = hsl, sum, na.rm = TRUE)
        tabland <- aggregate(landWt ~ Var1 + Var2 + Var3, data = cl_temp, sum, na.rm = TRUE)

        tabsamp$sampling <- tabsamp$lenNum / sum(tabsamp$lenNum) * 100
        tabland$landings <- tabland$landWt / sum(tabland$landWt) * 100

        tab <- merge(tabsamp, tabland, by = c("Var1", "Var2", "Var3"), all = TRUE)
        tab$year <- ye
        tab$sampling[is.na(tab$sampling)] <- 0
        tab$landings[is.na(tab$landings)] <- 0

        poinByStr <- tab[, c("Var1", "Var2", "Var3"), drop = FALSE]
        poinByStr$sampSize <- 0
        poinByStr$nMeasTotal <- 0
        poinByStr$MeanNmeasByStr <- 0
        poinByStr$cv <- NA_real_

        for (nr in seq_len(nrow(tab))) {
            loc_case <- tab[nr, c("Var1", "Var2", "Var3"), drop = FALSE]

            cs_loc <- cs_ye
            if (space == "GSA") {
                cs_loc <- subsetCOST(cs_loc, area %in% loc_case$Var1, table = "hh")
            } else if (space == "Country") {
                cs_loc <- subsetCOST(cs_loc, vslFlgCtry %in% loc_case$Var1, table = "hh")
            }

            if (time == "Q") {
                cs_loc <- subsetCOST(cs_loc, get_quarter(cs_loc@hh$date) %in% loc_case$Var2, table = "hh")
            } else if (time == "S") {
                cs_loc <- subsetCOST(cs_loc, get_semester(cs_loc@hh$date) %in% loc_case$Var2, table = "hh")
            }

            if (technical == "lev6") {
                cs_loc <- subsetCOST(cs_loc, foCatEu6 %in% loc_case$Var3, table = "hh")
            } else if (technical == "lev4") {
                cs_loc <- subsetCOST(cs_loc, foCatEu5 %in% loc_case$Var3, table = "hh")
            }

            if (nrow(cs_loc@hl) == 0L) {
                next
            }

            poinByStr$sampSize[nr] <- length(unique(cs_loc@hl$trpCode))
            poinByStr$nMeasTotal[nr] <- sum(cs_loc@hl$lenNum, na.rm = TRUE)

            if (poinByStr$sampSize[nr] > 0) {
                poinByStr$MeanNmeasByStr[nr] <- round(
                    poinByStr$nMeasTotal[nr] / poinByStr$sampSize[nr], 0
                )
            }

            CSvalstr <- csDataVal(cs_loc)
            CLvaly <- clDataVal(cl_ye)
            CSconsstr <- csDataCons(CSvalstr, stratification)
            CLconsy <- clDataCons(CLvaly, stratification)

            dbeOutputy <- dbeObject(
                desc = sppName,
                species = sppName,
                catchCat = "Lan",
                strataDesc = stratification,
                methodDesc = "analytical"
            )

            raisingy <- try(
                RaiseLgth(
                    dbeOutputy,
                    CSconsstr,
                    CLconsy,
                    spp = sppName,
                    taxon = sppName,
                    strDesc = stratification
                ),
                silent = TRUE
            )

            if (!inherits(raisingy, "try-error")) {
                poinByStr$cv[nr] <- raisingy@lenNum$DCRcvIndicator
            }
        }

        tab_year <- merge(tab, poinByStr, by = c("Var1", "Var2", "Var3"))
        tab_year$Species <- sppName
        past_list[[as.character(ye)]] <- tab_year
    }

    pastSituations <- data.table::rbindlist(past_list, fill = TRUE)
    past_file <- file.path(res_dir, "Past_situation.csv")
    write.table(pastSituations, file = past_file, sep = ";", row.names = FALSE, na = "")

    ## ---------------------------------------------------------------------------
    ## 2. Sampling vs landings plots
    ## ---------------------------------------------------------------------------
    sampling_plot_files <- character(0)

    if (nrow(pastSituations) > 0L) {
        vect_met <- unique(pastSituations$Var3)

        for (met in vect_met) {
            dat_met <- pastSituations[pastSituations$Var3 == met, , drop = FALSE]
            lim_max <- ceiling(max(c(dat_met$sampling, dat_met$landings), na.rm = TRUE))

            p <- ggplot2::ggplot(dat_met) +
                ggplot2::geom_point(
                    ggplot2::aes(x = landings, y = sampling, colour = as.factor(year)),
                    size = 2
                ) +
                ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.2) +
                ggplot2::facet_grid(Var1 ~ Var2, scales = "fixed") +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(size = 7),
                    axis.text.y = ggplot2::element_text(size = 7),
                    legend.title = ggplot2::element_blank(),
                    legend.position = "top"
                ) +
                ggplot2::labs(
                    title = paste(sppName, "-", met, "- Sampling vs Landings"),
                    x = "% Landings",
                    y = "% Sampling"
                ) +
                ggplot2::xlim(0, lim_max) +
                ggplot2::ylim(0, lim_max)

            file_name <- file.path(res_dir, paste0(sppShortCode, "_", met, "_Sampling_vs_Landing.jpg"))
            ggplot2::ggsave(file_name, p, width = 8, height = 10)
            sampling_plot_files <- c(sampling_plot_files, file_name)
        }
    }

    ## ---------------------------------------------------------------------------
    ## 3. Landings summary
    ## ---------------------------------------------------------------------------
    tab_land <- aggregate(
        CL@cl$landWt,
        by = list(CL@cl$area, CL@cl$taxon, CL@cl$quarter, CL@cl$year, CL@cl$foCatEu6),
        FUN = sum,
        na.rm = TRUE
    )
    colnames(tab_land) <- c("Area", "Species", "Quarter", "Year", "foCatEu6", "Value")

    landings_file <- file.path(res_dir, "Landings_by_year_quarter_foCatEu6.csv")
    write.table(tab_land, landings_file, sep = ";", row.names = FALSE, na = "")

    landings_plot_files <- character(0)

    for (gs in gsas) {
        tab_this_gsa <- tab_land[tab_land$Area == gs, , drop = FALSE]
        if (nrow(tab_this_gsa) == 0L) next

        tab_this_gsa$Value <- tab_this_gsa$Value / 1000

        p <- ggplot2::ggplot(tab_this_gsa, ggplot2::aes(x = Year, y = Value, fill = as.factor(Quarter))) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::ggtitle(paste("Landings (tons) by foCatEu6", gs)) +
            ggplot2::facet_grid(Species ~ foCatEu6, scales = "fixed")

        file_name <- file.path(res_dir, paste0("Landings_by_year_quarter_foCatEu6_", gs, ".jpg"))
        ggplot2::ggsave(file_name, p, width = 10, height = 6)
        landings_plot_files <- c(landings_plot_files, file_name)
    }

    ## ---------------------------------------------------------------------------
    ## 4. Trips / samples / measures summary
    ## ---------------------------------------------------------------------------
    hsl_sl <- merge(cs_spe@hh, cs_spe@sl)
    hsl_hl <- merge(cs_spe@hh, cs_spe@hl)

    T_final <- data.frame()
    for (gs in unique(hsl_sl$area)) {
        hslgs <- hsl_sl[hsl_sl$area == gs, , drop = FALSE]
        for (yy in sort(unique(hslgs$year))) {
            hsly <- hslgs[hslgs$year == yy, , drop = FALSE]
            t_1 <- data.frame(tapply(hsly$commCat, list(hsly$foCatEu6, get_quarter(hsly$date)), length))
            if (nrow(t_1) == 0L) next

            colnames(t_1) <- sort(unique(get_quarter(hsly$date)))
            rownames(t_1) <- paste(as.character(yy), rownames(t_1))

            to_add <- setdiff(1:4, as.integer(colnames(t_1)))
            if (length(to_add) > 0L) {
                ad <- data.frame(matrix(0, nrow = nrow(t_1), ncol = length(to_add)))
                colnames(ad) <- to_add
                t_1 <- cbind(t_1, ad)
                t_1 <- t_1[, order(as.integer(colnames(t_1))), drop = FALSE]
            }

            t_1 <- data.frame(
                Area = as.character(gs),
                Year = as.character(yy),
                Metier = substring(rownames(t_1), 6, nchar(rownames(t_1))),
                Species = sppName,
                t_1,
                row.names = NULL
            )
            T_final <- rbind(T_final, t_1)
        }
    }

    tab1 <- data.frame(
        Area = T_final$Area,
        Year = T_final$Year,
        Metier = T_final$Metier,
        Species = T_final$Species,
        SQ1 = T_final[, 5],
        SQ2 = T_final[, 6],
        SQ3 = T_final[, 7],
        SQ4 = T_final[, 8]
    )

    T_final2 <- data.frame()
    for (gs in unique(hsl_hl$area)) {
        hslgs <- hsl_hl[hsl_hl$area == gs, , drop = FALSE]
        for (yy in sort(unique(hslgs$year))) {
            hsl2y <- hslgs[hslgs$year == yy, , drop = FALSE]
            t_2 <- data.frame(tapply(hsl2y$lenNum, list(hsl2y$foCatEu6, get_quarter(hsl2y$date)), sum))
            if (nrow(t_2) == 0L) next

            colnames(t_2) <- sort(unique(get_quarter(hsl2y$date)))
            rownames(t_2) <- paste(as.character(yy), rownames(t_2))

            to_add <- setdiff(1:4, as.integer(colnames(t_2)))
            if (length(to_add) > 0L) {
                ad <- data.frame(matrix(0, nrow = nrow(t_2), ncol = length(to_add)))
                colnames(ad) <- to_add
                t_2 <- cbind(t_2, ad)
                t_2 <- t_2[, order(as.integer(colnames(t_2))), drop = FALSE]
            }

            t_2 <- data.frame(
                Area = as.character(gs),
                Year = as.character(yy),
                Metier = substring(rownames(t_2), 6, nchar(rownames(t_2))),
                Species = sppName,
                t_2,
                row.names = NULL
            )
            T_final2 <- rbind(T_final2, t_2)
        }
    }

    tab2 <- data.frame(
        Area = T_final2$Area,
        Year = T_final2$Year,
        Metier = T_final2$Metier,
        Species = T_final2$Species,
        NQ1 = T_final2[, 5],
        NQ2 = T_final2[, 6],
        NQ3 = T_final2[, 7],
        NQ4 = T_final2[, 8]
    )

    tab_TRIP <- tab2
    colnames(tab_TRIP)[5:8] <- c("TQ1", "TQ2", "TQ3", "TQ4")
    tab_TRIP[, 5:8] <- 0

    no_trips_mat <- CS@hh
    no_trips_mat$quarter <- get_quarter(no_trips_mat$date)

    no_trips_mat <- with(
        no_trips_mat,
        aggregate(year, by = list(area, year, foCatEu6, quarter, trpCode), FUN = "length")
    )
    colnames(no_trips_mat) <- c("area", "year", "foCatEu6", "quarter", "trpCode", "x")

    no_trips_mat_2 <- with(
        no_trips_mat,
        aggregate(year, by = list(area, year, foCatEu6, quarter), FUN = "length")
    )
    colnames(no_trips_mat_2) <- c("area", "year", "foCatEu6", "quarter", "no_trips")

    for (nr in seq_len(nrow(no_trips_mat_2))) {
        qq <- no_trips_mat_2$quarter[nr]
        idx <- tab_TRIP$Area == no_trips_mat_2$area[nr] &
            tab_TRIP$Year == no_trips_mat_2$year[nr] &
            tab_TRIP$Metier == no_trips_mat_2$foCatEu6[nr]

        if (qq == 1) tab_TRIP[idx, "TQ1"] <- no_trips_mat_2$no_trips[nr]
        if (qq == 2) tab_TRIP[idx, "TQ2"] <- no_trips_mat_2$no_trips[nr]
        if (qq == 3) tab_TRIP[idx, "TQ3"] <- no_trips_mat_2$no_trips[nr]
        if (qq == 4) tab_TRIP[idx, "TQ4"] <- no_trips_mat_2$no_trips[nr]
    }

    tab0 <- merge(tab_TRIP, merge(tab1, tab2), all = TRUE)
    tab0[is.na(tab0)] <- 0
    tab0[, 5:16] <- round(tab0[, 5:16], 0)

    tab0$annualT <- rowSums(tab0[, 5:8], na.rm = TRUE)
    tab0$annualS <- rowSums(tab0[, 9:12], na.rm = TRUE)
    tab0$annualN <- rowSums(tab0[, 13:16], na.rm = TRUE)

    tsm_file <- file.path(res_dir, "No_of_trips_samples_measures.csv")
    write.table(tab0, file = tsm_file, row.names = FALSE, sep = ";", na = "")

    ## ---------------------------------------------------------------------------
    ## Return
    ## ---------------------------------------------------------------------------
    invisible(list(
        species = target_SPECIES,
        output_dir = res_dir,
        past_situation = pastSituations,
        landings_summary = tab_land,
        trips_samples_measures = tab0,
        files = list(
            past_situation_csv = past_file,
            landings_csv = landings_file,
            trips_samples_measures_csv = tsm_file,
            sampling_vs_landings_plots = sampling_plot_files,
            landings_plots = landings_plot_files
        )
    ))
}
