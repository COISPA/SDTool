#' Create stratification variables Var1, Var2 and Var3
#'
#' Builds the stratification variables \code{Var1}, \code{Var2}, and
#' \code{Var3} in a sampling or landing data.frame, according to the selected
#' spatial, temporal, and technical stratification scheme.
#'
#' @param x A data.frame containing the variables needed for stratification.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#'
#' @return The input data.frame with three added columns:
#' \describe{
#'   \item{Var1}{space stratification variable}
#'   \item{Var2}{time stratification variable}
#'   \item{Var3}{technical stratification variable}
#' }
#'
#' @details
#' The function expects the following columns depending on the chosen options:
#' \itemize{
#'   \item \code{space = "GSA"}: column \code{area}
#'   \item \code{space = "Country"}: column \code{vslFlgCtry}
#'   \item \code{time = "Q"}: column \code{SampQuarter} or \code{quarter}
#'   \item \code{time = "S"}: column \code{SampSemester} or \code{semester}
#'   \item \code{technical = "lev6"}: column \code{foCatEu6}
#'   \item \code{technical = "lev4"}: column \code{foCatEu5}
#' }
#'
#' This helper is useful inside functions such as
#' \code{data_investigation()} when the same stratification logic must be
#' applied to different tables.
#'
#' @param is_landings Logical. If \code{FALSE}, the function looks for sampling
#'   time variables (\code{SampQuarter}, \code{SampSemester}). If \code{TRUE},
#'   it looks for landing time variables (\code{quarter}, \code{semester}).
#'
#' @examples
#' \dontrun{
#' x <- setVar123(
#'   x = x,
#'   space = "GSA",
#'   time = "Q",
#'   technical = "lev6"
#' )
#' }
#'
#' @keywords internal
setVar123 <- function(x,
                      space = c("GSA", "Country", "ALL"),
                      time = c("Q", "S", "Y"),
                      technical = c("lev6", "lev4", "NONE"),
                      is_landings = FALSE) {

    if (!is.data.frame(x)) {
        stop("`x` must be a data.frame.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    ## ---------------------------------------------------------------------------
    ## Var1: space
    ## ---------------------------------------------------------------------------
    if (space == "GSA") {
        if (!"area" %in% names(x)) {
            stop("Column `area` is required when `space = 'GSA'`.")
        }
        x$Var1 <- x$area
    }

    if (space == "Country") {
        if (!"vslFlgCtry" %in% names(x)) {
            stop("Column `vslFlgCtry` is required when `space = 'Country'`.")
        }
        x$Var1 <- x$vslFlgCtry
    }

    if (space == "ALL") {
        x$Var1 <- rep("ALL", nrow(x))
    }

    ## ---------------------------------------------------------------------------
    ## Var2: time
    ## ---------------------------------------------------------------------------
    if (time == "Q") {
        time_col <- if (is_landings) "quarter" else "SampQuarter"
        if (!time_col %in% names(x)) {
            stop("Column `", time_col, "` is required when `time = 'Q'`.")
        }
        x$Var2 <- x[[time_col]]
    }

    if (time == "S") {
        time_col <- if (is_landings) "semester" else "SampSemester"
        if (!time_col %in% names(x)) {
            stop("Column `", time_col, "` is required when `time = 'S'`.")
        }
        x$Var2 <- x[[time_col]]
    }

    if (time == "Y") {
        x$Var2 <- rep(1, nrow(x))
    }

    ## ---------------------------------------------------------------------------
    ## Var3: technical
    ## ---------------------------------------------------------------------------
    if (technical == "lev6") {
        if (!"foCatEu6" %in% names(x)) {
            stop("Column `foCatEu6` is required when `technical = 'lev6'`.")
        }
        x$Var3 <- x$foCatEu6
    }

    if (technical == "lev4") {
        if (!"foCatEu5" %in% names(x)) {
            stop("Column `foCatEu5` is required when `technical = 'lev4'`.")
        }
        x$Var3 <- x$foCatEu5
    }

    if (technical == "NONE") {
        x$Var3 <- rep("ALL", nrow(x))
    }

    x
}
