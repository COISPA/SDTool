#' Select rows for one trip code and rename the trip code
#'
#' Extracts the rows of a table corresponding to one selected trip code and
#' replaces that trip code with the corresponding value from an extension list.
#'
#' @param i Integer index identifying the position of the trip code to be
#'   selected in \code{selected.tripcode} and the replacement code in
#'   \code{extention.list}.
#' @param table A data.frame containing at least the column \code{trpCode}.
#' @param extention.list A vector of replacement trip codes.
#' @param selected.tripcode A vector of trip codes to be selected from
#'   \code{table}.
#'
#' @return A data.frame containing the subset of \code{table} for the selected
#'   trip code, with \code{trpCode} renamed using \code{extention.list[i]}.
#'   If no rows are found, an empty data.frame is returned.
#'
#' @examples
#' \dontrun{
#' out <- select.and.rename.tripcode2(
#'   i = 1,
#'   table = hl_reduced,
#'   extention.list = extentionList,
#'   selected.tripcode = selected_samples
#' )
#' }
#'
#' @keywords internal
#' @export
select.and.rename.tripcode2 <- function(i,
                                        table,
                                        extention.list,
                                        selected.tripcode) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(table)) {
        stop("`table` must be a data.frame.")
    }

    if (!"trpCode" %in% names(table)) {
        stop("`table` must contain a `trpCode` column.")
    }

    if (length(i) != 1L || !is.numeric(i) || is.na(i)) {
        stop("`i` must be a single non-missing numeric index.")
    }

    i <- as.integer(i)

    if (i < 1L || i > length(selected.tripcode)) {
        stop("`i` is out of bounds for `selected.tripcode`.")
    }

    if (i < 1L || i > length(extention.list)) {
        stop("`i` is out of bounds for `extention.list`.")
    }

    ## ---------------------------------------------------------------------------
    ## Select and rename trip code
    ## ---------------------------------------------------------------------------
    tripcode <- selected.tripcode[i]

    sub.table <- table[table$trpCode == tripcode, , drop = FALSE]

    if (nrow(sub.table) > 0L) {
        sub.table$trpCode <- extention.list[i]
    }

    return(sub.table)
}
