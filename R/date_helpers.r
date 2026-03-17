#' Convert to Date
#' @param x A date vector.
#' @return A Date vector.
#' @keywords internal
#' @export
to_date <- function(x) {

    # se è già Date → ritorna subito
    if (inherits(x, "Date")) {
        return(x)
    }

    x <- as.character(x)

    # prova ISO: YYYY-MM-DD
    d <- as.Date(x, format = "%Y-%m-%d")

    # prova DD-MM-YYYY
    idx <- is.na(d)
    if (any(idx)) {
        d[idx] <- as.Date(x[idx], format = "%d-%m-%Y")
    }

    # prova DD/MM/YYYY
    idx <- is.na(d)
    if (any(idx)) {
        d[idx] <- as.Date(x[idx], format = "%d/%m/%Y")
    }

    # fallback generico (ultima spiaggia)
    idx <- is.na(d)
    if (any(idx)) {
        d[idx] <- as.Date(x[idx])
    }

    return(d)
}

#' Get quarter from dates
#' @param x A date vector.
#' @return An integer vector.
#' @keywords internal
#' @export
get_quarter <- function(x) {
    x <- to_date(x)
    ((as.integer(format(x, "%m")) - 1L) %/% 3L) + 1L
}

#' Get semester from dates
#' @param x A date vector.
#' @return An integer vector.
#' @keywords internal
#' @export
get_semester <- function(x) {
    ifelse(get_quarter(x) %in% c(1, 2), 1L, 2L)
}
