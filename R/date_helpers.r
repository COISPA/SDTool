#' Convert to Date
#' @param x A date vector.
#' @return A Date vector.
#' @keywords internal
to_date <- function(x) {
    if (inherits(x, "Date")) {
        return(x)
    }
    as.Date(x, format = "%d/%m/%Y")
}

#' Get quarter from dates
#' @param x A date vector.
#' @return An integer vector.
#' @keywords internal
get_quarter <- function(x) {
    x <- to_date(x)
    ((as.integer(format(x, "%m")) - 1L) %/% 3L) + 1L
}

#' Get semester from dates
#' @param x A date vector.
#' @return An integer vector.
#' @keywords internal
get_semester <- function(x) {
    ifelse(get_quarter(x) %in% c(1, 2), 1L, 2L)
}
