#' Select samples within one stratum
#'
#' Selects a set of trip codes (\code{trpCode}) belonging to a given
#' stratification cell, based on the available sampling pool and the minimum /
#' maximum sample size constraints stored in \code{sampSizeRange}.
#'
#' @param indStr Integer index identifying the row of \code{availStrat}
#'   corresponding to the stratum to be sampled.
#' @param availStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, and \code{Var3}.
#' @param availSamp A data.frame containing the available samples. It must
#'   include \code{trpCode} and the columns required by the selected
#'   stratification.
#' @param sampSizeRange A data.frame containing the minimum and maximum sample
#'   sizes by stratum. It must include the columns \code{Var1}, \code{Var2},
#'   \code{Var3}, \code{min_max}, and \code{solutions}.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#'
#' @return A character vector of selected \code{trpCode}s. If no valid sample
#'   size range is available for the selected stratum, the function returns
#'   \code{NA_character_}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item extracts the target stratum from \code{availStrat},
#'   \item identifies the pool of trip codes matching the spatial, temporal,
#'   and technical stratification variables,
#'   \item intersects the three pools,
#'   \item retrieves the allowed minimum and maximum sample size from
#'   \code{sampSizeRange},
#'   \item randomly selects a number of trip codes within that range.
#' }
#'
#' Sampling is performed with replacement, consistently with the original
#' implementation.
#'
#' @examples
#' \dontrun{
#' sel <- selectSampStr(
#'   indStr = 1,
#'   availStrat = availStrat,
#'   availSamp = availSamp,
#'   sampSizeRange = sampSizeRange,
#'   space = "GSA",
#'   time = "Q",
#'   technical = "lev6"
#' )
#' }
#'
#' @keywords internal
selectSampStr <- function(indStr,
                          availStrat,
                          availSamp,
                          sampSizeRange,
                          space = c("GSA", "Country", "ALL"),
                          time = c("Q", "S", "Y"),
                          technical = c("lev6", "lev4", "NONE")) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(availStrat)) {
        stop("`availStrat` must be a data.frame.")
    }

    if (!is.data.frame(availSamp)) {
        stop("`availSamp` must be a data.frame.")
    }

    if (!is.data.frame(sampSizeRange)) {
        stop("`sampSizeRange` must be a data.frame.")
    }

    if (!all(c("Var1", "Var2", "Var3") %in% names(availStrat))) {
        stop("`availStrat` must contain columns `Var1`, `Var2`, and `Var3`.")
    }

    if (!all(c("Var1", "Var2", "Var3", "min_max", "solutions") %in% names(sampSizeRange))) {
        stop("`sampSizeRange` must contain columns `Var1`, `Var2`, `Var3`, `min_max`, and `solutions`.")
    }

    if (!"trpCode" %in% names(availSamp)) {
        stop("`availSamp` must contain column `trpCode`.")
    }

    if (length(indStr) != 1L || !is.numeric(indStr) || indStr < 1 || indStr > nrow(availStrat)) {
        stop("`indStr` must be a valid row index of `availStrat`.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    ## ---------------------------------------------------------------------------
    ## Target stratum
    ## ---------------------------------------------------------------------------
    case <- availStrat[indStr, c("Var1", "Var2", "Var3"), drop = FALSE]

    ## ---------------------------------------------------------------------------
    ## Build pool by space
    ## ---------------------------------------------------------------------------
    if (space == "GSA") {
        if (!"area" %in% names(availSamp)) {
            stop("Column `area` is required when `space = 'GSA'`.")
        }
        sampPoolStr1 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$area == case$Var1]),
            group1 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (space == "Country") {
        if (!"vslFlgCtry" %in% names(availSamp)) {
            stop("Column `vslFlgCtry` is required when `space = 'Country'`.")
        }
        sampPoolStr1 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$vslFlgCtry == case$Var1]),
            group1 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (space == "ALL") {
        sampPoolStr1 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode),
            group1 = 1,
            stringsAsFactors = FALSE
        )
    }

    ## ---------------------------------------------------------------------------
    ## Build pool by time
    ## ---------------------------------------------------------------------------
    if (time == "Q") {
        if (!"SampQuarter" %in% names(availSamp)) {
            stop("Column `SampQuarter` is required when `time = 'Q'`.")
        }
        sampPoolStr2 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$SampQuarter == case$Var2]),
            group2 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (time == "S") {
        if (!"SampSemester" %in% names(availSamp)) {
            stop("Column `SampSemester` is required when `time = 'S'`.")
        }
        sampPoolStr2 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$SampSemester == case$Var2]),
            group2 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (time == "Y") {
        sampPoolStr2 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode),
            group2 = 1,
            stringsAsFactors = FALSE
        )
    }

    ## ---------------------------------------------------------------------------
    ## Build pool by technical variable
    ## ---------------------------------------------------------------------------
    if (technical == "lev6") {
        if (!"foCatEu6" %in% names(availSamp)) {
            stop("Column `foCatEu6` is required when `technical = 'lev6'`.")
        }
        sampPoolStr3 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$foCatEu6 == case$Var3]),
            group3 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (technical == "lev4") {
        if (!"foCatEu5" %in% names(availSamp)) {
            stop("Column `foCatEu5` is required when `technical = 'lev4'`.")
        }
        sampPoolStr3 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$foCatEu5 == case$Var3]),
            group3 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (technical == "NONE") {
        sampPoolStr3 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode),
            group3 = 1,
            stringsAsFactors = FALSE
        )
    }

    ## ---------------------------------------------------------------------------
    ## Intersect the three pools
    ## ---------------------------------------------------------------------------
    sampPoolStr <- merge(
        sampPoolStr1,
        merge(sampPoolStr2, sampPoolStr3, all = TRUE),
        all = TRUE
    )

    sampPoolStr$sum_groups <- rowSums(sampPoolStr[, c("group1", "group2", "group3")], na.rm = TRUE)

    sampPoolStr <- as.character(
        sampPoolStr$sampPoolStr[sampPoolStr$sum_groups == 3]
    )

    ## ---------------------------------------------------------------------------
    ## Retrieve min / max sample size for the selected stratum
    ## ---------------------------------------------------------------------------
    minSize <- sampSizeRange$solutions[
        as.character(sampSizeRange$Var1) == as.character(case$Var1[1]) &
            as.character(sampSizeRange$Var2) == as.character(case$Var2[1]) &
            as.character(sampSizeRange$Var3) == as.character(case$Var3[1]) &
            sampSizeRange$min_max == "min"
    ]

    maxSize <- sampSizeRange$solutions[
        as.character(sampSizeRange$Var1) == as.character(case$Var1[1]) &
            as.character(sampSizeRange$Var2) == as.character(case$Var2[1]) &
            as.character(sampSizeRange$Var3) == as.character(case$Var3[1]) &
            sampSizeRange$min_max == "max"
    ]

    ## ---------------------------------------------------------------------------
    ## Random selection
    ## ---------------------------------------------------------------------------
    if (length(minSize) == 0L || length(maxSize) == 0L || is.na(minSize[1]) || is.na(maxSize[1])) {
        return(NA_character_)
    }

    minSize <- as.integer(minSize[1])
    maxSize <- as.integer(maxSize[1])

    if (maxSize == minSize && minSize != 0) {
        maxSize <- maxSize + 1L
    }

    if (length(sampPoolStr) == 0L || maxSize < minSize) {
        return(NA_character_)
    }

    sampSize <- sample(minSize:maxSize, 1, replace = FALSE)
    selectSamp <- sample(sampPoolStr, sampSize, replace = TRUE)

    return(selectSamp)
}
