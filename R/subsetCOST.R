#' Subset a COST catch sampling object while preserving table links
#'
#' Subsets a \code{csData} or \code{csDataVal} object on one selected table
#' and propagates the subset to the linked tables in the COST hierarchy.
#'
#' The function preserves consistency among the \code{tr}, \code{hh},
#' \code{sl}, \code{hl}, and \code{ca} slots by subsetting upward and
#' downward across linked tables.
#'
#' @param x A \code{csData} or \code{csDataVal} object.
#' @param subset A logical expression indicating elements or rows to keep.
#'   The expression is evaluated within the selected table.
#' @param ... Additional arguments. Currently ignored.
#' @param table A character string indicating the table on which the subset
#'   must be applied. Allowed values are \code{"tr"}, \code{"hh"},
#'   \code{"sl"}, \code{"hl"}, and \code{"ca"}.
#'
#' @return A \code{csData} object, or a \code{csDataVal} object if the input
#'   \code{x} is of class \code{csDataVal}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item extracts the COST tables from the input object;
#'   \item applies the subset condition to the selected table;
#'   \item propagates the filtering upward and downward through linked tables;
#'   \item rebuilds a consistent COST object.
#' }
#'
#' Empty tables are replaced by the corresponding empty slots from
#' \code{csData()} in order to return a valid COST object.
#'
#' @examples
#' \dontrun{
#' x_sub <- subsetCOST(x, year == 2020, table = "hh")
#' x_sub2 <- subsetCOST(x, spp == "Merluccius merluccius", table = "hl")
#' }
#'
#' @export
subsetCOST <- function(x, subset, ..., table = "tr") {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!inherits(x, c("csData", "csDataVal"))) {
        stop("`x` must be a 'csData' or 'csDataVal' object.")
    }

    allowed_tables <- c("tr", "hh", "sl", "hl", "ca")
    if (!table %in% allowed_tables) {
        stop("`table` must be one of: ", paste(allowed_tables, collapse = ", "))
    }

    isVal <- inherits(x, "csDataVal")

    ## ---------------------------------------------------------------------------
    ## Extract tables
    ## ---------------------------------------------------------------------------
    tr <- tr(x)
    hh <- hh(x)
    sl <- sl(x)
    hl <- slSex(sl, hl(x))   # add artificial sex field for key consistency
    ca <- csData()@ca        # original code uses an empty ca template

    ## ---------------------------------------------------------------------------
    ## Helper: build primary / foreign keys
    ## ---------------------------------------------------------------------------
    fpKey <- function(tab, colIndex, sep = ":-:") {
        key <- tab[, colIndex, drop = FALSE]
        key <- apply(key, 1, paste, collapse = sep)
        key <- gsub("[[:space:]]", "", key)
        key
    }

    ## ---------------------------------------------------------------------------
    ## Identify tr records linked to ca versus hh
    ## ---------------------------------------------------------------------------
    indca <- fpKey(tr, 1:6) %in% fpKey(ca, 1:6)

    trca <- tr[indca, , drop = FALSE]
    trca$N <- which(indca)

    trhh <- tr[!indca, , drop = FALSE]
    trhh$N <- which(!indca)

    ## ---------------------------------------------------------------------------
    ## Apply subset to the selected table
    ## ---------------------------------------------------------------------------
    e <- substitute(subset)

    df0 <- switch(
        table,
        tr = tr,
        hh = hh,
        sl = sl,
        hl = hl,
        ca = ca
    )

    r <- eval(e, df0, parent.frame(n = 1))
    if (!is.logical(r)) {
        stop("`subset` must evaluate to a logical vector.")
    }

    df_sub <- df0[r, , drop = FALSE]

    if (table == "tr") tr <- df_sub
    if (table == "hh") hh <- df_sub
    if (table == "sl") sl <- df_sub
    if (table == "hl") hl <- df_sub
    if (table == "ca") ca <- df_sub

    ## ---------------------------------------------------------------------------
    ## Key indexes according to table hierarchy
    ## ---------------------------------------------------------------------------
    Up <- matrix(
        c("trhh", "hh", "sl", "trca",
          "1:6",  "1:7", "1:14", "1:6"),
        nrow = 2,
        byrow = TRUE
    )
    dimnames(Up) <- list(c("tab", "index"), c("hh", "sl", "hl", "ca"))

    Down <- matrix(
        c("hh",  "ca",  "sl",  "hl",
          "1:6", "1:6", "1:7", "1:14"),
        nrow = 2,
        byrow = TRUE
    )
    dimnames(Down) <- list(c("tab", "index"), c("tr", "tr", "hh", "sl"))

    ## ---------------------------------------------------------------------------
    ## Recursive propagation of subset through linked tables
    ## ---------------------------------------------------------------------------
    indSub <- FALSE

    subs <- function(tabName, tabKey) {
        if (tabName %in% dimnames(tabKey)[[2]]) {
            indSub <<- TRUE

            mat <- tabKey[, dimnames(tabKey)[[2]] %in% tabName, drop = FALSE]

            eval(
                parse(
                    text = paste(
                        mat["tab", ], " <<- ",
                        mat["tab", ], "[fpKey(",
                        mat["tab", ], ",", mat["index", ],
                        ") %in% fpKey(",
                        tabName, ",", mat["index", ],
                        "), ]",
                        sep = "",
                        collapse = ";"
                    )
                )
            )

            Recall(mat["tab", 1], tabKey)
        }
    }

    ## ---------------------------------------------------------------------------
    ## Apply propagation: upward first, then downward from tr
    ## ---------------------------------------------------------------------------
    subs(table, Up)

    if (indSub) {
        tr <- rbind.data.frame(trhh, trca)
        tr <- tr[order(tr$N), 1:(ncol(tr) - 1), drop = FALSE]
    }

    subs("tr", Down)

    ## ---------------------------------------------------------------------------
    ## Restore hl structure
    ## ---------------------------------------------------------------------------
    if (nrow(hl) > 0) {
        hl$sex <- hl$lsex
        hl <- hl[, -ncol(hl), drop = FALSE]
    }

    ## ---------------------------------------------------------------------------
    ## Rebuild valid COST object
    ## ---------------------------------------------------------------------------
    if (nrow(tr) < 1) tr <- csData()@tr
    if (nrow(hh) < 1) hh <- csData()@hh
    if (nrow(sl) < 1) sl <- csData()@sl
    if (nrow(hl) < 1) hl <- csData()@hl
    if (nrow(ca) < 1) ca <- csData()@ca

    res <- csData(
        tr = tr,
        hh = hh,
        sl = sl,
        hl = hl,
        ca = ca,
        desc = x@desc
    )

    if (isVal) {
        res <- csDataVal(res)
    }

    return(res)
}
