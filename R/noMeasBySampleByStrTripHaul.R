#' Count measures by trip and haul within one stratum
#'
#' Subsets a COST catch sampling object to one stratification cell and computes
#' the number of measured individuals by trip and by trip-haul combination.
#'
#' @param i Integer index identifying the row of \code{tabStrat} corresponding
#'   to the target stratum.
#' @param tabStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, and \code{Var3}.
#' @param CSsp A \code{csData} or \code{csDataVal} object containing sampling
#'   data for one species.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{Var1}{space stratum}
#'   \item{Var2}{time stratum}
#'   \item{Var3}{technical stratum}
#'   \item{trpCode}{trip code}
#'   \item{commCat}{commercial category, fixed to \code{"ALL"}}
#'   \item{noMeasures}{number of measured individuals}
#'   \item{trpCodeHaul}{trip-haul identifier}
#' }
#'
#' If no valid records are found, the function returns an empty data.frame with
#' the same columns.
#'
#' @examples
#' \dontrun{
#' out <- noMeasBySampleByStrTripHaul(
#'   i = 1,
#'   tabStrat = tabStrat,
#'   CSsp = CSsp,
#'   space = "GSA",
#'   time = "Q",
#'   technical = "lev6"
#' )
#' }
#'
#' @keywords internal
noMeasBySampleByStrTripHaul <- function(i,
                                        tabStrat,
                                        CSsp,
                                        space = c("GSA", "Country", "ALL"),
                                        time = c("Q", "S", "Y"),
                                        technical = c("lev6", "lev4", "NONE")) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(tabStrat)) {
        stop("`tabStrat` must be a data.frame.")
    }

    if (!all(c("Var1", "Var2", "Var3") %in% names(tabStrat))) {
        stop("`tabStrat` must contain columns `Var1`, `Var2`, and `Var3`.")
    }

    if (!inherits(CSsp, c("csData", "csDataVal"))) {
        stop("`CSsp` must be a 'csData' or 'csDataVal' object.")
    }

    if (length(i) != 1L || !is.numeric(i) || is.na(i)) {
        stop("`i` must be a single non-missing numeric index.")
    }

    i <- as.integer(i)

    if (i < 1L || i > nrow(tabStrat)) {
        stop("`i` is out of bounds for `tabStrat`.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    ## ---------------------------------------------------------------------------
    ## Helpers
    ## ---------------------------------------------------------------------------
    to_date <- function(x) {
        if (inherits(x, "Date")) {
            return(x)
        }
        as.Date(x, format = "%d/%m/%Y")
    }

    get_quarter <- function(x) {
        x <- to_date(x)
        ((as.integer(format(x, "%m")) - 1L) %/% 3L) + 1L
    }

    get_semester <- function(x) {
        x <- to_date(x)
        ((as.integer(format(x, "%m")) - 1L) %/% 6L) + 1L
    }

    empty_result <- data.frame(
        Var1 = character(0),
        Var2 = character(0),
        Var3 = character(0),
        trpCode = character(0),
        commCat = character(0),
        noMeasures = numeric(0),
        trpCodeHaul = character(0),
        stringsAsFactors = FALSE
    )

    ## ---------------------------------------------------------------------------
    ## Target stratum
    ## ---------------------------------------------------------------------------
    case <- tabStrat[i, , drop = FALSE]

    ## ---------------------------------------------------------------------------
    ## Subset by space
    ## ---------------------------------------------------------------------------
    if (space == "GSA") {
        CSstr1 <- subsetCOST(CSsp, area %in% case$Var1, table = "hh")
    }

    if (space == "Country") {
        CSstr1 <- subsetCOST(CSsp, vslFlgCtry %in% case$Var1, table = "hh")
    }

    if (space == "ALL") {
        CSstr1 <- CSsp
    }

    ## ---------------------------------------------------------------------------
    ## Subset by time
    ## ---------------------------------------------------------------------------
    if (time == "Q") {
        CSstr2 <- subsetCOST(
            CSstr1,
            get_quarter(CSstr1@hh$date) %in% case$Var2,
            table = "hh"
        )
    }

    if (time == "S") {
        CSstr2 <- subsetCOST(
            CSstr1,
            get_semester(CSstr1@hh$date) %in% case$Var2,
            table = "hh"
        )
    }

    if (time == "Y") {
        CSstr2 <- CSstr1
    }

    ## ---------------------------------------------------------------------------
    ## Subset by technical variable
    ## ---------------------------------------------------------------------------
    if (technical == "lev6") {
        CSstr3 <- subsetCOST(CSstr2, foCatEu6 %in% case$Var3, table = "hh")
    }

    if (technical == "lev4") {
        CSstr3 <- subsetCOST(CSstr2, foCatEu5 %in% case$Var3, table = "hh")
    }

    if (technical == "NONE") {
        CSstr3 <- CSstr2
    }

    ## ---------------------------------------------------------------------------
    ## Extract HL table
    ## ---------------------------------------------------------------------------
    if (nrow(CSstr3@hl) == 0L) {
        return(empty_result)
    }

    HL_df <- data.frame(CSstr3@hl, stringsAsFactors = FALSE)

    if (!all(c("trpCode", "staNum", "lenNum") %in% names(HL_df))) {
        stop("`CSstr3@hl` must contain columns `trpCode`, `staNum`, and `lenNum`.")
    }

    HL_df$commCat <- "ALL"
    HL_df$trpCodeHaul <- paste(HL_df$trpCode, HL_df$staNum, sep = "_")
    HL_df$lenNum <- as.numeric(as.character(HL_df$lenNum))

    ## ---------------------------------------------------------------------------
    ## Summarise number of measures
    ## ---------------------------------------------------------------------------
    result_grouped <- dplyr::group_by(HL_df, trpCode, trpCodeHaul, commCat)

    result <- data.frame(
        dplyr::summarise(
            result_grouped,
            noMeasures = sum(lenNum, na.rm = TRUE)
        ),
        stringsAsFactors = FALSE
    )

    if (nrow(result) == 0L) {
        return(empty_result)
    }

    result <- data.frame(
        Var1 = rep(as.character(case$Var1[1]), nrow(result)),
        Var2 = rep(as.character(case$Var2[1]), nrow(result)),
        Var3 = rep(as.character(case$Var3[1]), nrow(result)),
        trpCode = result$trpCode,
        commCat = result$commCat,
        noMeasures = result$noMeasures,
        trpCodeHaul = result$trpCodeHaul,
        stringsAsFactors = FALSE
    )

    return(result)
}
