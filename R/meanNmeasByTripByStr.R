#' Compute mean number of measures by trip within one stratum
#'
#' Subsets a COST catch sampling object to one stratification cell and computes
#' the mean number of measured individuals per trip.
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
#' @return A one-row data.frame containing:
#' \describe{
#'   \item{Var1}{space stratum}
#'   \item{Var2}{time stratum}
#'   \item{Var3}{technical stratum}
#'   \item{meanNmeasByTrip}{mean number of measures per trip}
#' }
#'
#' If the selected stratum contains no valid trips, the function returns
#' \code{NA} for \code{meanNmeasByTrip}.
#'
#' @examples
#' \dontrun{
#' out <- meanNmeasByTripByStr(
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
#' @export
meanNmeasByTripByStr <- function(i,
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

    # ## ---------------------------------------------------------------------------
    # ## Helpers
    # ## ---------------------------------------------------------------------------
    # to_date <- function(x) {
    #     if (inherits(x, "Date")) {
    #         return(x)
    #     }
    #     as.Date(x, format = "%d-%m-%Y")
    # }
    #
    # get_quarter <- function(x) {
    #     x <- to_date(x)
    #     ((as.integer(format(x, "%m")) - 1L) %/% 3L) + 1L
    # }
    #
    # get_semester <- function(x) {
    #     ifelse(get_quarter(x) %in% c(1, 2), 1L, 2L)
    # }

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
    ## Compute mean number of measures by trip
    ## ---------------------------------------------------------------------------
    if (nrow(CSstr3@hl) == 0L || all(is.na(CSstr3@hl$trpCode))) {
        meanNmeasByTrip <- NA_real_
    } else {
        sampSize <- length(unique(CSstr3@hl$trpCode))
        nMeasTotal <- sum(as.numeric(as.character(CSstr3@hl$lenNum)), na.rm = TRUE)

        if (sampSize == 0) {
            meanNmeasByTrip <- NA_real_
        } else {
            meanNmeasByTrip <- nMeasTotal / sampSize
        }
    }

    result <- data.frame(
        Var1 = as.character(case$Var1[1]),
        Var2 = as.character(case$Var2[1]),
        Var3 = as.character(case$Var3[1]),
        meanNmeasByTrip = meanNmeasByTrip,
        stringsAsFactors = FALSE
    )

    return(result)
}
