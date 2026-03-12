#' Convert an RCG Catch Sampling table to COST csData (TR/HH/SL/HL/CA)
#'
#' Converts a catch sampling table (RCG data call format) into SDEF COST CS tables
#' (TR, HH, SL, HL, CA) and returns a COST `csData` object.
#'
#' The input CS can be provided in two modalities:
#' 1) aggregated by length class (fishId and individual weight missing)
#' 2) individual data (fishId provided)
#'
#' @param CS A data.frame containing the RCG Catch Sampling table.
#' @param dataset_proj Character. Project name stored in the output (default: "STREAM project").
#' @param log_var Logical. If TRUE, writes a log CSV marking duplicated keys.
#' @param log_varFilePath Optional path for the log file (CSV, `;` separated). If NULL and `log_var=TRUE`,
#'   a temporary file is used.
#' @param bad_rm Logical. If TRUE, remove duplicated records (where implemented); if FALSE, duplicates trigger an error.
#' @param write_csv Logical. If TRUE, writes SDEF CS-TR/HH/SL/HL/CA CSV files in `out_dir`.
#' @param out_dir Output directory for CSV/RDS exports (created if missing).
#' @param write_rds Logical. If TRUE, saves the resulting `csData` object as RDS to `rds_path`.
#' @param rds_path Path for the RDS output (used if `write_rds=TRUE`).
#' @param fill_na Value used to replace NA in numeric/character fields before processing (default: -1).
#'
#' @return A COST `csData` object.
#' @export
#'
#' @examples
#' data(cs_gsa99_dummy)
#' costCS <- convert_cs_to_costCS(cs_gsa99_dummy, bad_rm = TRUE, log_var = FALSE)
#' costCS

convert_cs_to_costCS <- function(
        CS,
        dataset_proj = "STREAM project",
        log_var = TRUE,
        log_varFilePath = NULL,
        bad_rm = FALSE,
        write_csv = FALSE,
        out_dir = ".",
        write_rds = TRUE,
        rds_path = file.path(out_dir, "costCS.rds"),
        fill_na = -1
) {
    stopifnot(is.data.frame(CS))

    dir.create("COSTobjects", recursive = FALSE, showWarnings = FALSE)

    out_dir=file.path(out_dir, "COSTobjects")
    # ---- helper: primary key fingerprint
    fpKey <- function(tab, colIndex, sep = ":-:") {
        key <- tab[, colIndex, drop = FALSE]
        key <- apply(key, 1, paste, collapse = sep)
        key <- gsub("[[:space:]]", "", key)
        key
    }
    # DEBUG
    type_mismatch <- function(df, template_df, slot_name = "slot") {
        common <- intersect(names(df), names(template_df))
        yours <- vapply(df[common], function(x) class(x)[1], character(1))
        templ <- vapply(template_df[common], function(x) class(x)[1], character(1))
        bad <- common[yours != templ]

        list(
            slot = slot_name,
            mismatches = data.frame(
                column = bad,
                yours = yours[bad],
                template = templ[bad],
                stringsAsFactors = FALSE
            )
        )
    }


    align_to_template <- function(df, template_df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)

        # 1) add missing columns with correct NA type
        for (nm in names(template_df)) {
            if (!nm %in% names(df)) {
                df[[nm]] <- template_df[[nm]][NA_integer_]
            }
        }

        # 2) coerce existing columns to template types
        for (nm in names(template_df)) {
            tgt <- class(template_df[[nm]])[1]
            df[[nm]] <- switch(
                tgt,
                integer = as.integer(df[[nm]]),
                numeric = as.numeric(df[[nm]]),
                character = as.character(df[[nm]]),
                logical = as.logical(df[[nm]]),
                factor = as.factor(df[[nm]]),
                Date = as.Date(df[[nm]]),
                POSIXct = as.POSIXct(df[[nm]]),
                df[[nm]]
            )
        }

        # 3) order columns as template
        df <- df[, names(template_df), drop = FALSE]
        df
    }

    # FINE DEBUG
    # ---- init log
    if (isTRUE(log_var)) {
        log_varCs <- CS
        log_varCs$duplicated_TR <- FALSE
        log_varCs$duplicated_HH <- FALSE
        log_varCs$duplicated_SL <- FALSE
        log_varCs$duplicated_HL <- FALSE
        log_varCs$duplicated_CA <- FALSE
        log_varMsg <- ""
    } else {
        log_varMsg <- " Consider using log_var=TRUE."
    }

    error <- FALSE

    # ---- standardize column names (case-insensitive mapping)
    nms <- names(CS)
    low <- tolower(nms)

    rename_if_present <- function(from_low, to) {
        idx <- which(low == from_low)
        if (length(idx) > 0) names(CS)[idx[1]] <<- to
    }

    # core
    rename_if_present("sampling.type", "sampType")
    rename_if_present("sampling_type", "sampType")
    rename_if_present("samplingtype", "sampType")

    rename_if_present("flag.country", "vslFlgCtry")
    rename_if_present("flag_country", "vslFlgCtry")
    rename_if_present("flagcountry", "vslFlgCtry")

    rename_if_present("year", "year")

    rename_if_present("trip.code", "trpCode")
    rename_if_present("trip_code", "trpCode")
    rename_if_present("tripcode", "trpCode")

    rename_if_present("number.of.sets...hauls.on.trip", "foNum")
    rename_if_present("number_of_setshauls", "foNum")
    rename_if_present("nsets", "foNum")

    rename_if_present("days.at.sea", "daysAtSea")
    rename_if_present("days_at_sea", "daysAtSea")
    rename_if_present("daysatsea", "daysAtSea")

    rename_if_present("sampling.method", "sampMeth")
    rename_if_present("sampling_method", "sampMeth")

    rename_if_present("aggregation.level", "aggLev")
    rename_if_present("aggregation_level", "aggLev")
    rename_if_present("aggregationlevel", "aggLev")

    rename_if_present("station.number", "staNum")
    rename_if_present("station_number", "staNum")
    rename_if_present("stationnumber", "staNum")

    rename_if_present("catch.registration", "catReg")
    rename_if_present("catch_registration", "catReg")
    rename_if_present("catchregistration", "catReg")

    rename_if_present("species.registration", "sppReg")
    rename_if_present("species_registration", "sppReg")
    rename_if_present("speciesregistration", "sppReg")

    rename_if_present("date", "date")
    rename_if_present("area", "area")

    rename_if_present("fishing.activity.category.national", "foCatNat")
    rename_if_present("fishingactivitycategorynational", "foCatNat")
    rename_if_present("fac_national", "foCatNat")

    rename_if_present("fishing.activity.category.european.lvl.6", "foCatEu6")
    rename_if_present("metier", "foCatEu6")
    rename_if_present("fac_ec_lvl6", "foCatEu6")

    rename_if_present("species", "spp")

    rename_if_present("catch.category", "catchCat")
    rename_if_present("catch_category", "catchCat")
    rename_if_present("catchcategory", "catchCat")

    rename_if_present("commercial.size.category", "commCat")
    rename_if_present("commercial_size_category", "commCat")

    rename_if_present("commercial.size.category.scale", "commCatScl")
    rename_if_present("commercial_size_category_scale", "commCatScl")
    rename_if_present("commercialsizecategoryscale", "commCatScl")

    rename_if_present("weight", "wt")

    rename_if_present("subsample.weight", "subSampWt")
    rename_if_present("subsample_weight", "subSampWt")
    rename_if_present("subsampleweight", "subSampWt")

    rename_if_present("length.code", "lenCode")
    rename_if_present("length_code", "lenCode")
    rename_if_present("lengthcode", "lenCode")

    rename_if_present("length.class", "lenCls")
    rename_if_present("length_class", "lenCls")
    rename_if_present("lengthclass", "lenCls")

    rename_if_present("number.at.length", "lenNum")
    rename_if_present("number_at_length", "lenNum")
    rename_if_present("numberatlength", "lenNum")

    # 2018+ extra fields
    rename_if_present("fish.id", "fishId")
    rename_if_present("fish_id", "fishId")
    rename_if_present("fishid", "fishId")

    rename_if_present("individual.weight", "indWt")
    rename_if_present("individual_weight", "indWt")
    rename_if_present("individualweight", "indWt")

    rename_if_present("harbour", "harbour")

    rename_if_present("duration.of.fishing.operation", "foDur")
    rename_if_present("duration_of_fishing_operation", "foDur")
    rename_if_present("durationoffishing_operation", "foDur")

    rename_if_present("initial.latitude", "latIni")
    rename_if_present("initial_latitude", "latIni")
    rename_if_present("initiallatitude", "latIni")

    rename_if_present("initial.longitude", "lonIni")
    rename_if_present("initial_longitude", "lonIni")

    rename_if_present("final.latitude", "latFin")
    rename_if_present("final_latitude", "latFin")
    rename_if_present("finallatitude", "latFin")

    rename_if_present("final.longitude", "lonFin")
    rename_if_present("final_longitude", "lonFin")
    rename_if_present("finallongitude", "lonFin")

    rename_if_present("depth.of.fishing.operation", "foDep")
    rename_if_present("depth_of_fishing_operation", "foDep")
    rename_if_present("depthoffishingoperation", "foDep")

    rename_if_present("water.depth", "waterDep")
    rename_if_present("water_depth", "waterDep")
    rename_if_present("waterdepth", "waterDep")

    rename_if_present("sex", "sex")

    rename_if_present("maturity.method", "matMeth")
    rename_if_present("maturity_method", "matMeth")
    rename_if_present("maturitymethod", "matMeth")

    rename_if_present("maturity.scale", "matScale")
    rename_if_present("maturity_scale", "matScale")
    rename_if_present("maturityscale", "matScale")

    rename_if_present("maturity.stage", "matStage")
    rename_if_present("maturity_stage", "matStage")
    rename_if_present("maturitystage", "matStage")

    rename_if_present("ageing.method", "ageMeth")
    rename_if_present("ageing_method", "ageMeth")
    rename_if_present("ageingmethod", "ageMeth")

    rename_if_present("age", "age")

    # ---- fill NA (matches your original behavior)
    CS[is.na(CS)] <- fill_na

    # ---- normalize aggLev and species capitalization
    CS$aggLev <- as.character(CS$aggLev)
    CS$aggLev[CS$aggLev %in% c("t", "TRUE", TRUE)] <- "T"

    CS$spp <- vapply(CS$spp, function(x) {
        x <- as.character(x)
        paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    }, character(1))

    # ---- add project & landing country
    CS$proj <- dataset_proj
    CS$landCtry <- CS$vslFlgCtry

    # ---- define PKs
    trPk <- c("sampType", "vslFlgCtry", "year", "trpCode", "proj", "landCtry")
    trOther <- c("foNum", "daysAtSea", "sampMeth", "harbour")

    hhPk <- c(trPk, "staNum")
    hhOther <- c(
        "aggLev", "catReg", "sppReg", "date", "area",
        "foCatNat", "foCatEu6", "foDur", "latIni", "lonIni", "latFin", "lonFin", "foDep", "waterDep"
    )

    slPk <- c(hhPk, "spp", "catchCat", "commCat", "commCatScl", "sex")
    slOther <- c("wt", "subSampWt", "lenCode")

    hlPk <- c(slPk, "lenCls")
    hlOther <- c("lenNum")

    caPk <- c(hlPk, "age", "area", "fishId", "matStage")
    caOther <- c("matMeth", "matScale", "ageMeth", "indWt", "lenCode")

    allFields <- unique(c(caPk, trOther, hhOther, slOther, hlOther, caOther))
    missingFields <- allFields[!allFields %in% names(CS)]
    if (length(missingFields) > 0) {
        stop("Missing fields: ", paste(missingFields, collapse = ", "))
    }

    # ---- TR integrity
    csTr <- unique(CS[, c(trPk, trOther), drop = FALSE])
    trPkV <- fpKey(csTr, trPk)
    trPkVDup <- trPkV %in% trPkV[duplicated(trPkV)]
    if (any(trPkVDup)) {
        if (isTRUE(log_var)) {
            log_varCs$duplicated_TR <- fpKey(CS, trPk) %in% trPkV[trPkVDup]
        }
        if (isTRUE(bad_rm)) {
            message("Integrity problem for CS/TR, ", sum(trPkVDup), " row(s) removed.", log_varMsg)
            csTr <- csTr[!trPkVDup, , drop = FALSE]
        } else {
            error <- TRUE
            message("Integrity problem for CS/TR, ", sum(trPkVDup), " row(s) concerned.", log_varMsg)
        }
    }

    # ---- HH integrity
    csHh <- unique(CS[, c(hhPk, hhOther), drop = FALSE])
    hhPkV <- fpKey(csHh, hhPk)
    hhPkVDup <- hhPkV %in% hhPkV[duplicated(hhPkV)]
    if (any(hhPkVDup)) {
        if (isTRUE(log_var)) {
            log_varCs$duplicated_HH <- fpKey(CS, hhPk) %in% hhPkV[hhPkVDup]
        }
        if (isTRUE(bad_rm)) {
            message("Integrity problem for CS/HH, ", sum(hhPkVDup), " row(s) removed.", log_varMsg)
            csHh <- csHh[!hhPkVDup, , drop = FALSE]
        } else {
            error <- TRUE
            message("Integrity problem for CS/HH, ", sum(hhPkVDup), " row(s) concerned.", log_varMsg)
        }
    }

    # ---- SL integrity
    csSl <- unique(CS[, c(slPk, slOther), drop = FALSE])
    slPkV <- fpKey(csSl, slPk)
    slPkVDup <- slPkV %in% slPkV[duplicated(slPkV)]
    if (any(slPkVDup)) {
        if (isTRUE(log_var)) {
            log_varCs$duplicated_SL <- fpKey(CS, slPk) %in% slPkV[slPkVDup]
        }
        if (isTRUE(bad_rm)) {
            message("Integrity problem for CS/SL, ", sum(slPkVDup), " row(s) removed.", log_varMsg)
            csSl <- csSl[!slPkVDup, , drop = FALSE]
        } else {
            error <- TRUE
            message("Integrity problem for CS/SL, ", sum(slPkVDup), " row(s) concerned.", log_varMsg)
        }
    }

    # ---- Aggregate by length (HL-level) like your script
    CS_aggregated_by_length <- stats::aggregate(
        CS$lenNum,
        by = list(
            CS$sampType, CS$vslFlgCtry, CS$year, CS$trpCode, CS$proj, CS$landCtry, CS$harbour,
            CS$foNum, CS$daysAtSea, CS$sampMeth, CS$aggLev, CS$staNum,
            CS$foDur, CS$foDep, CS$catReg, CS$sppReg, CS$date, CS$area,
            CS$foCatNat, CS$foCatEu6, CS$spp, CS$catchCat, CS$wt, CS$subSampWt,
            CS$sex, CS$lenCode, CS$lenCls, CS$commCatScl, CS$commCat
        ),
        FUN = "sum"
    )

    colnames(CS_aggregated_by_length) <- c(
        "sampType", "vslFlgCtry", "year", "trpCode", "proj", "landCtry", "harbour",
        "foNum", "daysAtSea", "sampMeth", "aggLev", "staNum", "foDur", "foDep",
        "catReg", "sppReg", "date", "area", "foCatNat", "foCatEu6",
        "spp", "catchCat", "wt", "subSampWt", "sex", "lenCode", "lenCls",
        "commCatScl", "commCat", "lenNum"
    )

    csHl <- unique(CS_aggregated_by_length[, c(hlPk, hlOther), drop = FALSE])
    hlPkV <- fpKey(csHl, hlPk)
    hlPkVDup <- hlPkV %in% hlPkV[duplicated(hlPkV)]
    if (any(hlPkVDup)) {
        if (isTRUE(log_var)) {
            log_varCs$duplicated_HL <- fpKey(CS, hlPk) %in% hlPkV[hlPkVDup]
        }
        # HL duplicates are treated as error in your original script
        error <- TRUE
        message("Integrity problem for CS/HL, ", sum(hlPkVDup), " row(s) concerned.", log_varMsg)
    }

    # ---- CA integrity (from raw CS)
    csCa <- unique(CS[, c(caPk, caOther), drop = FALSE])
    caPkV <- fpKey(csCa, caPk)
    caPkVDup <- caPkV %in% caPkV[duplicated(caPkV)]
    if (any(caPkVDup)) {
        if (isTRUE(log_var)) {
            log_varCs$duplicated_CA <- fpKey(CS, caPk) %in% caPkV[caPkVDup]
        }
        if (isTRUE(bad_rm)) {
            message("Integrity problem for CS/CA, ", sum(caPkVDup), " row(s) removed.", log_varMsg)
            csCa <- csCa[!caPkVDup, , drop = FALSE]
        } else {
            error <- TRUE
            message("Integrity problem for CS/CA, ", sum(caPkVDup), " row(s) concerned.", log_varMsg)
        }
    }

    # ---- write log
    if (isTRUE(log_var)) {
        if (is.null(log_varFilePath)) log_varFilePath <- tempfile(fileext = ".csv")
        utils::write.table(log_varCs, file = log_varFilePath, row.names = FALSE, sep = ";")
        message("Log file: ", log_varFilePath)
    }

    if (isTRUE(error)) {
        stop("An error occurred in the transformation. See log file if enabled.")
    }

    # ---- Get COST CS table column names from csData class
    # ---- Get COST CS table column names from csData class (robust across COST versions)
    cs_class <- NULL

    for (pkg in c("COSTdbe", "COSTcore")) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ns <- asNamespace(pkg)
            if (methods::isClass("csData", where = ns)) {
                cs_class <- methods::getClass("csData", where = ns)
                break
            }
        }
    }

    if (is.null(cs_class)) {
        stop("S4 class 'csData' not found. Check COST packages installed/loaded (COSTdbe/COSTcore).")
    }

    obj <- methods::new(cs_class)

    TR.col <- names(obj@tr)
    HH.col <- names(obj@hh)
    SL.col <- names(obj@sl)
    HL.col <- names(obj@hl)
    CA.col <- names(obj@ca)
# DEBUG
    tr_template <- obj@tr
    hh_template <- obj@hh
    sl_template <- obj@sl
    hl_template <- obj@hl
    ca_template <- obj@ca
# FINE DEBUG
# ---- Align tables to COST schemas (typed + coerced)
    costCS.TR <- align_to_template(csTr, tr_template)
    costCS.HH <- align_to_template(csHh, hh_template)
    costCS.SL <- align_to_template(csSl, sl_template)
    costCS.HL <- align_to_template(csHl, hl_template)

    # ---- Build CA
    if (all(CS$fishId == fill_na)) {
        # aggregated modality: expand individuals from HL/SL
        ca_cols_extra <- c(caPk, caOther, "subSampWt", "lenNum")
        ca_sl_hl <- unique(CS[, ca_cols_extra, drop = FALSE])

        # indWt = subSampWt / lenNum
        ca_sl_hl$indWt <- ifelse(ca_sl_hl$lenNum == 0, fill_na, ca_sl_hl$subSampWt / ca_sl_hl$lenNum)

        # replicate rows lenNum times
        reps <- as.integer(ca_sl_hl$lenNum)
        reps[is.na(reps) | reps < 0] <- 0
        ca_expanded <- ca_sl_hl[rep(seq_len(nrow(ca_sl_hl)), reps), , drop = FALSE]

        # assign fishId sequence
        if (nrow(ca_expanded) > 0) {
            ca_expanded$fishId <- seq_len(nrow(ca_expanded))
        }

        # keep only CA columns
        costCS.CA <- ca_expanded[, intersect(names(ca_expanded), CA.col), drop = FALSE]
        missCA <- CA.col[!CA.col %in% names(costCS.CA)]
        costCS.CA[, missCA] <- NA
        costCS.CA <- costCS.CA[, CA.col, drop = FALSE]

        # match original script: indWt set to -1 in aggregated case
        costCS.CA$indWt <- fill_na
    } else {
        missCA <- CA.col[!CA.col %in% names(csCa)]
        costCS.CA <- csCa
        costCS.CA[, missCA] <- NA

        costCS.CA$proj <- dataset_proj
        costCS.CA$landCtry <- costCS.CA$vslFlgCtry
        costCS.CA <- costCS.CA[, CA.col, drop = FALSE]

        # ensure fishId is sequential
        costCS.CA$fishId <- seq_len(nrow(costCS.CA))
    }
    costCS.CA <- align_to_template(costCS.CA, ca_template)


    # ---- write outputs
    if (isTRUE(write_csv)) {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        utils::write.table(costCS.TR, file.path(out_dir, "SDEF CS-TR data.csv"), sep = ";", row.names = FALSE)
        utils::write.table(costCS.HH, file.path(out_dir, "SDEF CS-HH data.csv"), sep = ";", row.names = FALSE)
        utils::write.table(costCS.SL, file.path(out_dir, "SDEF CS-SL data.csv"), sep = ";", row.names = FALSE)
        utils::write.table(costCS.HL, file.path(out_dir, "SDEF CS-HL data.csv"), sep = ";", row.names = FALSE)
        utils::write.table(costCS.CA, file.path(out_dir, "SDEF CS-CA data.csv"), sep = ";", row.names = FALSE)
    }

    # DEBUG

    mm_tr <- type_mismatch(costCS.TR, tr_template, "tr")
    if (nrow(mm_tr$mismatches) > 0) {
        print(mm_tr$mismatches)
        stop("Type mismatch in TR slot (see printed table above).")
    }
    mm_hh <- type_mismatch(costCS.HH, hh_template, "hh")
    if (nrow(mm_hh$mismatches) > 0) {
        print(mm_hh$mismatches)
        stop("Type mismatch in HH slot (see printed table above).")
    }

        # FINE DEBUG

    # ---- build csData object without calling unexported constructor
    costCS <- methods::new(
        "csData",
        tr = costCS.TR,
        hh = costCS.HH,
        sl = costCS.SL,
        hl = costCS.HL,
        ca = costCS.CA
    )

    methods::validObject(costCS)

    if (isTRUE(write_rds)) {
        dir.create(dirname(rds_path), showWarnings = FALSE, recursive = TRUE)
        saveRDS(costCS, rds_path)
    }

    return(costCS)


}
