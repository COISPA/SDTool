#' Convert a CL landing table to a COST clData object
#'
#' This function standardizes column names, checks key integrity (duplicates),
#' formats taxon capitalization, builds the SDEF CL data.frame and returns a
#' COST `clData` object.
#'
#' @param CL A data.frame containing the CL landing table.
#' @param log Logical. If TRUE, a log table is written (see `logFilePath`).
#' @param logFilePath Optional path for the log CSV file. If NULL and `log=TRUE`,
#'   a temporary file is used.
#' @param bad_rm Logical. If TRUE, duplicated primary keys are removed.
#'   If FALSE, duplicates trigger an error.
#' @param write_csv Logical. If TRUE, writes the SDEF CL CSV to `csv_path`.
#' @param csv_path Path where the SDEF CL CSV is written (used if `write_csv=TRUE`).
#' @param write_rds Logical. If TRUE, saves the resulting `costCL` to `rds_path`.
#' @param rds_path Path for the RDS output (used if `write_rds=TRUE`).
#'
#' @return A `clData` object.
#' @export
#'
#' @examples
#' # costCL <- convert_cl_to_costCL(cl_gsa99_dummy, log=TRUE, bad_rm=TRUE)
convert_cl_to_costCL <- function(
        CL,
        log = TRUE,
        logFilePath = NULL,
        bad_rm = FALSE,
        out_dir = ".",
        write_csv = FALSE,
        csv_path = file.path(out_dir,"SDEF CL data.csv"),
        write_rds = FALSE,
        rds_path = file.path(out_dir,"costCL.rds")

) {
    stopifnot(is.data.frame(CL))

    # ---- helper: primary key fingerprint
    fpKey <- function(tab, colIndex, sep = ":-:") {
        key <- tab[, colIndex, drop = FALSE]
        key <- apply(key, 1, paste, collapse = sep)
        key <- gsub("[[:space:]]", "", key)
        key
    }

    align_to_template <- function(df, template_df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)

        # add missing columns with correct NA type
        for (nm in names(template_df)) {
            if (!nm %in% names(df)) {
                df[[nm]] <- template_df[[nm]][NA_integer_]
            }
        }

        # coerce existing columns to template types
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

        # order columns
        df <- df[, names(template_df), drop = FALSE]
        df
    }


    # ---- log init
    if (isTRUE(log)) {
        logCl <- CL
        logCl$duplicated <- FALSE
        logMsg <- ""
    } else {
        logMsg <- " Consider using the log=TRUE parameter."
    }
    error <- FALSE

    # ---- standardize column names (case-insensitive mapping)
    nms <- names(CL)
    low <- tolower(nms)

    rename_if_present <- function(from_low, to) {
        idx <- which(low == from_low)
        if (length(idx) > 0) names(CL)[idx[1]] <<- to
    }

    rename_if_present("flag.country", "vslFlgCtry")
    rename_if_present("flag_country", "vslFlgCtry")

    rename_if_present("year", "year")
    rename_if_present("quarter", "quarter")
    rename_if_present("month", "month")
    rename_if_present("area", "area")
    rename_if_present("species", "taxon")

    rename_if_present("fishing.activity.category.national", "foCatNat")
    rename_if_present("fac_national", "foCatNat")

    rename_if_present("fishing.activity.category.european.lvl.6", "foCatEu6")
    rename_if_present("fac_ec_lvl6", "foCatEu6")

    rename_if_present("harbour", "harbour")

    rename_if_present("official.landings.weight", "landWt")
    rename_if_present("official_landings_weight", "landWt")

    rename_if_present("official.landings.value", "landValue")
    rename_if_present("official_landings_value", "landValue")

    # ---- primary keys & required fields
    clPk <- c("vslFlgCtry", "year", "quarter", "month", "area", "taxon", "foCatNat", "foCatEu6")
    clOther <- c("landWt", "landValue")
    allFields <- c(clPk, clOther)

    missingFields <- allFields[!allFields %in% names(CL)]
    if (length(missingFields) > 0) {
        stop("Missing fields: ", paste(missingFields, collapse = ", "))
    }

    clPkV <- fpKey(CL, clPk)
    clPkVDup <- clPkV %in% clPkV[duplicated(clPkV)]

    # ---- integrity check
    if (any(clPkVDup)) {
        if (isTRUE(log)) {
            # segnala righe duplicate nel log
            logCl$duplicated <- clPkVDup
        }

        if (isTRUE(bad_rm)) {
            message("Integrity problem for CL, ", sum(clPkVDup),
                    " row(s) removed.", logMsg)
            CL <- CL[!clPkVDup, , drop = FALSE]
        } else {
            error <- TRUE
            message("Integrity problem for CL, ", sum(clPkVDup),
                    " row(s) concerned.", logMsg)
        }
    }

    # ---- write log file
    if (isTRUE(log)) {
        if (is.null(logFilePath)) logFilePath <- tempfile(fileext = ".csv")
        utils::write.table(logCl, file = logFilePath, row.names = FALSE, sep = ";")
        message("Log file: ", logFilePath)
    }

    if (isTRUE(error)) stop("Stop on reported errors.")

    # ---- formatting: taxon capitalization
    CL$taxon <- vapply(CL$taxon, function(x) {
        x <- as.character(x)
        paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    }, character(1))

    # ---- build SDEF CL data.frame
    clDf <- data.frame(
        landCtry = NA,
        vslFlgCtry = CL$vslFlgCtry,
        year = CL$year,
        quarter = CL$quarter,
        month = CL$month,
        area = CL$area,
        rect = NA,
        subRect = NA,
        taxon = CL$taxon,
        landCat = NA,
        commCatScl = NA,
        commCat = NA,
        foCatNat = CL$foCatNat,
        foCatEu5 = NA,
        foCatEu6 = CL$foCatEu6,
        harbour = CL$harbour,
        vslLenCat = NA,
        unallocCatchWt = NA,
        misRepCatchWt = NA,
        landWt = CL$landWt,
        landMult = NA,
        landValue = CL$landValue,
        stringsAsFactors = FALSE
    )

    if (isTRUE(write_csv)) {
        utils::write.table(clDf, file = csv_path, sep = ";", row.names = FALSE)
    }

    # ---- create COST object

    # build clData object without calling unexported constructor
    cl_class <- NULL
    for (pkg in c("COSTdbe", "COSTcore")) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ns <- asNamespace(pkg)
            if (methods::isClass("clData", where = ns)) {
                cl_class <- methods::getClass("clData", where = ns)
                break
            }
        }
    }
    if (is.null(cl_class)) {
        stop("S4 class 'clData' not found. Check COST packages installed/loaded (COSTdbe/COSTcore).")
    }

    # template of the cl slot (correct names + types)
    cl_template <- methods::new(cl_class)@cl

    # align + coerce
    clDf <- align_to_template(clDf, cl_template)

    costCL <- methods::new(cl_class)
    costCL@cl <- clDf
    methods::validObject(costCL)



    if (isTRUE(write_rds)) {
        saveRDS(costCL, rds_path)
    }

    costCL
}
