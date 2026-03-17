#' Find local maxima in a numeric vector
#'
#' Identifies the positions of local maxima in a numeric vector.
#'
#' A local maximum is detected as the end of an increasing run in the input
#' vector. This helper is mainly intended for internal use when analysing
#' density curves or other ordered numeric sequences.
#'
#' @param x A numeric vector.
#'
#' @return An integer vector containing the indices of local maxima in
#'   \code{x}. If no local maxima are found, an empty integer vector is
#'   returned.
#'
#' @details
#' The function works by detecting changes in the sign of successive
#' differences. It is designed for ordered numeric vectors, such as the
#' \code{y} values returned by \code{density()}.
#'
#' If the first two values are equal, the first candidate maximum is removed
#' to avoid a false detection at the beginning of the series.
#'
#' @examples
#' x <- c(1, 3, 2, 4, 4, 1, 5, 2)
#' localMaxima(x)
#'
#' @keywords internal
#' @export
localMaxima <- function(x) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.numeric(x)) {
        stop("`x` must be a numeric vector.")
    }

    if (length(x) < 2L) {
        return(integer(0))
    }

    ## ---------------------------------------------------------------------------
    ## Detect local maxima
    ## ---------------------------------------------------------------------------
    y <- diff(c(-Inf, x)) > 0
    y <- rle(y)$lengths
    y <- cumsum(y)
    y <- y[seq.int(1L, length(y), 2L)]

    if (length(x) >= 2L && x[[1]] == x[[2]] && length(y) > 0L) {
        y <- y[-1]
    }

    as.integer(y)
}
