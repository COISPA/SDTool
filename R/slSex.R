#' Add sex information from SL table to HL table
#'
#' Merges the \code{sex} field from the \code{sl} table into the \code{hl} table,
#' creating an auxiliary \code{lsex} field to preserve the original \code{hl}
#' sex information.
#'
#' This function is mainly used internally to ensure that the \code{sex} field
#' from the sampling level (\code{sl}) is available in the length level
#' (\code{hl}) and can be treated as part of the full key.
#'
#' @param slTab A data.frame corresponding to the COST \code{sl} table.
#' @param hlTab A data.frame corresponding to the COST \code{hl} table.
#'
#' @return A data.frame corresponding to the \code{hl} table enriched with:
#' \describe{
#'   \item{sex}{the sex field imported from the \code{sl} table}
#'   \item{lsex}{the original sex field from the \code{hl} table}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item checks consistency between \code{slTab} and \code{hlTab};
#'   \item selects the key columns from \code{slTab};
#'   \item copies \code{slTab$sex} into a new \code{lsex} field;
#'   \item merges \code{slTab} into \code{hlTab};
#'   \item swaps \code{sex} and \code{lsex} so that:
#'   \itemize{
#'     \item \code{sex} contains the value coming from \code{slTab},
#'     \item \code{lsex} stores the original \code{hlTab$sex}.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' hl_with_sex <- slSex(slTab, hlTab)
#' }
#'
#' @export
slSex <- function(slTab, hlTab) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(slTab) || !is.data.frame(hlTab)) {
        stop("`slTab` and `hlTab` must be data.frame objects.")
    }

    if (ncol(slTab) < 1 || ncol(hlTab) < 1) {
        stop("`slTab` and `hlTab` must have at least one column.")
    }

    if (names(slTab)[1] != names(hlTab)[1]) {
        stop("Tables must be consistent: first column names do not match.")
    }

    if (!"sex" %in% names(slTab)) {
        stop("`slTab` must contain a `sex` column.")
    }

    if (!"sex" %in% names(hlTab)) {
        stop("`hlTab` must contain a `sex` column.")
    }

    ## ---------------------------------------------------------------------------
    ## Identify key columns according to COST table structure
    ## ---------------------------------------------------------------------------
    ind <- NULL

    if (names(slTab)[1] == "sampType") {
        ind <- 1:14
    }

    if (names(slTab)[1] == "PSUid") {
        ind <- 1:15
    }

    if (is.null(ind)) {
        stop("Wrong input tables: unrecognized COST table structure.")
    }

    ## ---------------------------------------------------------------------------
    ## Prepare SL table for merge
    ## ---------------------------------------------------------------------------
    slTab <- slTab[, ind, drop = FALSE]
    slTab$lsex <- slTab$sex

    ## Keep original order of HL table
    hlTab$N <- seq_len(nrow(hlTab))

    ## ---------------------------------------------------------------------------
    ## Merge SL sex information into HL
    ## ---------------------------------------------------------------------------
    hlTab <- merge(hlTab, slTab, all.x = TRUE)

    ## Restore original order
    hlTab <- hlTab[order(hlTab$N), , drop = FALSE]
    hlTab <- hlTab[, -match("N", names(hlTab)), drop = FALSE]

    ## ---------------------------------------------------------------------------
    ## Swap sex fields:
    ## - sex  = value from SL
    ## - lsex = original value from HL
    ## ---------------------------------------------------------------------------
    original_hl_sex <- as.character(hlTab$sex)
    hlTab$sex <- as.character(hlTab$lsex)
    hlTab$lsex <- original_hl_sex

    rownames(hlTab) <- seq_len(nrow(hlTab))

    return(hlTab)
}
