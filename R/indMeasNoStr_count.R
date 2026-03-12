#' Count measured individuals by stratum in a bootstrap iteration
#'
#' For each stratum in \code{availStrat}, the function selects a pool of trip
#' codes according to \code{selectSampStr()}, rebuilds a reduced COST catch
#' sampling object, and computes the total number of measured individuals
#' (\code{lenNum}) in the selected subset.
#'
#' @param iter Integer identifying the bootstrap iteration. Currently not used
#'   internally, but retained for compatibility with iterative workflows.
#' @param availStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, and \code{Var3}.
#' @param availSamp A data.frame containing the available samples, used by
#'   \code{selectSampStr()}.
#' @param sampSizeRange A data.frame defining the admissible sample-size ranges
#'   by stratum, used by \code{selectSampStr()}.
#' @param CSobject A \code{csData} or \code{csDataVal} object containing the
#'   sampling tables.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#'
#' @return A data.frame with one row per stratum and the columns:
#' \describe{
#'   \item{Var1}{space stratum}
#'   \item{Var2}{time stratum}
#'   \item{Var3}{technical stratum}
#'   \item{nInd}{total number of measured individuals in the selected samples}
#' }
#'
#' If no valid trip-code selection is available for a stratum, \code{nInd} is
#' set to \code{NA}.
#'
#' @examples
#' \dontrun{
#' out <- indMeasNoStr_count(
#'   iter = 1,
#'   availStrat = availStrat,
#'   availSamp = availSamp,
#'   sampSizeRange = sampSizeRange,
#'   CSobject = CS_boot,
#'   space = "GSA",
#'   time = "Q",
#'   technical = "lev6"
#' )
#' }
#'
#' @keywords internal
indMeasNoStr_count <- function(iter,
                               availStrat,
                               availSamp,
                               sampSizeRange,
                               CSobject,
                               space = c("GSA", "Country", "ALL"),
                               time = c("Q", "S", "Y"),
                               technical = c("lev6", "lev4", "NONE")) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(availStrat)) {
        stop("`availStrat` must be a data.frame.")
    }

    if (!all(c("Var1", "Var2", "Var3") %in% names(availStrat))) {
        stop("`availStrat` must contain columns `Var1`, `Var2`, and `Var3`.")
    }

    if (!is.data.frame(availSamp)) {
        stop("`availSamp` must be a data.frame.")
    }

    if (!is.data.frame(sampSizeRange)) {
        stop("`sampSizeRange` must be a data.frame.")
    }

    if (!inherits(CSobject, c("csData", "csDataVal"))) {
        stop("`CSobject` must be a 'csData' or 'csDataVal' object.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    ## ---------------------------------------------------------------------------
    ## Extract COST tables
    ## ---------------------------------------------------------------------------
    tr <- CSobject@tr
    hh <- CSobject@hh
    sl <- CSobject@sl
    hl <- CSobject@hl

    ## ---------------------------------------------------------------------------
    ## Initialise result container
    ## ---------------------------------------------------------------------------
    out_list <- vector("list", nrow(availStrat))

    ## ---------------------------------------------------------------------------
    ## Loop over strata
    ## ---------------------------------------------------------------------------
    for (indStr in seq_len(nrow(availStrat))) {

        tripCodeSelection <- selectSampStr(
            indStr = indStr,
            availStrat = availStrat,
            availSamp = availSamp,
            sampSizeRange = sampSizeRange,
            space = space,
            time = time,
            technical = technical
        )

        if (length(tripCodeSelection) > 0L &&
            !all(is.na(tripCodeSelection))) {

            tripCodeSelection <- tripCodeSelection[!is.na(tripCodeSelection)]

            ## generate replacement trip codes
            extention.list <- seq_along(tripCodeSelection)

            ## subset COST tables
            select.tr <- do.call(
                rbind,
                lapply(
                    seq_along(tripCodeSelection),
                    select.and.rename.tripcode2,
                    table = tr,
                    extention.list = extention.list,
                    selected.tripcode = tripCodeSelection
                )
            )

            select.hh <- do.call(
                rbind,
                lapply(
                    seq_along(tripCodeSelection),
                    select.and.rename.tripcode2,
                    table = hh,
                    extention.list = extention.list,
                    selected.tripcode = tripCodeSelection
                )
            )

            select.sl <- do.call(
                rbind,
                lapply(
                    seq_along(tripCodeSelection),
                    select.and.rename.tripcode2,
                    table = sl,
                    extention.list = extention.list,
                    selected.tripcode = tripCodeSelection
                )
            )

            select.hl <- do.call(
                rbind,
                lapply(
                    seq_along(tripCodeSelection),
                    select.and.rename.tripcode2,
                    table = hl,
                    extention.list = extention.list,
                    selected.tripcode = tripCodeSelection
                )
            )

            ## handle empty subsets safely
            if (is.null(select.hl) || nrow(select.hl) == 0L) {
                result <- NA_real_
            } else {
                subCS <- csData(
                    tr = select.tr,
                    hh = select.hh,
                    sl = select.sl,
                    hl = select.hl
                )

                result <- sum(as.numeric(as.character(subCS@hl$lenNum)), na.rm = TRUE)
            }

        } else {
            result <- NA_real_
        }

        out_list[[indStr]] <- data.frame(
            Var1 = as.character(availStrat$Var1[indStr]),
            Var2 = as.character(availStrat$Var2[indStr]),
            Var3 = as.character(availStrat$Var3[indStr]),
            nInd = result,
            stringsAsFactors = FALSE
        )
    }

    nMeas <- do.call(rbind, out_list)

    return(nMeas)
}
