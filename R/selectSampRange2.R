#' Select acceptable sample-size ranges from CV curves
#'
#' Identifies acceptable sample-size ranges for a given stratum based on the
#' empirical relationship between sample size and coefficient of variation (CV).
#' The function filters the CV results for one stratum, removes missing values
#' and strong outliers, estimates a density of CV values, detects local maxima,
#' and derives candidate sample-size solutions.
#'
#' @param indStr Integer index identifying the row of \code{availStrat}
#'   corresponding to the target stratum.
#' @param availStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, and \code{Var3}.
#' @param cvFunS2 A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, \code{Var3}, \code{sampSize}, \code{cv},
#'   and optionally \code{recyclingRate}.
#' @param table_min_accepted_sample_size A data.frame containing the minimum
#'   accepted sample size by stratum. It must contain the columns
#'   \code{space}, \code{time}, \code{technical}, and \code{min_accepted_ss}.
#' @param output_dir Directory where the density plot will be saved.
#' @return A data.frame with the selected sample-size solutions and associated
#'   summary statistics. The returned object includes the columns:
#'   \code{Var1}, \code{Var2}, \code{Var3}, \code{solutions}, \code{maxCV},
#'   \code{minCV}, \code{maxRR}, \code{minRR}, \code{meanRR}, and
#'   \code{noIterations}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item extracts the CV results for the selected stratum;
#'   \item removes missing CV values;
#'   \item orders observations by sample size;
#'   \item removes strong outliers before the CV peak;
#'   \item estimates the density of CV values;
#'   \item finds local maxima of the density curve;
#'   \item derives candidate sample-size ranges and keeps the accepted ones.
#' }
#'
#' If no valid density can be estimated, or if all CV values are missing, the
#' function returns a one-row data.frame with \code{NA} solutions.
#'
#' @examples
#' \dontrun{
#' out <- selectSampRange2(
#'   indStr = 1,
#'   availStrat = availStrat,
#'   cvFunS2 = cvFunS2,
#'   table_min_accepted_sample_size = table_min_accepted_sample_size,
#'   output_dir = tempdir(),
#'   rs_shortcode = "HKE"
#' )
#' }
#'
#' @keywords internal
#' @export
selectSampRange2 <- function(indStr,
                             availStrat,
                             cvFunS2,
                             table_min_accepted_sample_size,
                             output_dir) {

    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(availStrat)) {
        stop("`availStrat` must be a data.frame.")
    }

    if (!is.data.frame(cvFunS2)) {
        stop("`cvFunS2` must be a data.frame.")
    }

    if (!is.data.frame(table_min_accepted_sample_size)) {
        stop("`table_min_accepted_sample_size` must be a data.frame.")
    }

    if (!all(c("Var1", "Var2", "Var3") %in% names(availStrat))) {
        stop("`availStrat` must contain columns `Var1`, `Var2`, and `Var3`.")
    }

    if (!all(c("Var1", "Var2", "Var3", "sampSize", "cv") %in% names(cvFunS2))) {
        stop("`cvFunS2` must contain columns `Var1`, `Var2`, `Var3`, `sampSize`, and `cv`.")
    }

    if (!all(c("space", "time", "technical", "min_accepted_ss") %in% names(table_min_accepted_sample_size))) {
        stop("`table_min_accepted_sample_size` must contain columns `space`, `time`, `technical`, and `min_accepted_ss`.")
    }

    if (length(indStr) != 1L || !is.numeric(indStr) || indStr < 1 || indStr > nrow(availStrat)) {
        stop("`indStr` must be a valid row index of `availStrat`.")
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    ## ---------------------------------------------------------------------------
    ## Helper: local maxima
    ## ---------------------------------------------------------------------------
    localMaxima <- function(x) {
        if (length(x) < 3) {
            return(integer(0))
        }
        which(diff(sign(diff(x))) == -2) + 1
    }

    ## ---------------------------------------------------------------------------
    ## Default empty result
    ## ---------------------------------------------------------------------------
    empty_result <- data.frame(
        Var1 = availStrat$Var1[indStr],
        Var2 = availStrat$Var2[indStr],
        Var3 = availStrat$Var3[indStr],
        solutions = NA,
        maxCV = NA,
        minCV = NA,
        maxRR = NA,
        minRR = NA,
        meanRR = NA,
        noIterations = NA,
        stringsAsFactors = FALSE
    )

    ## ---------------------------------------------------------------------------
    ## Extract CV results for selected stratum
    ## ---------------------------------------------------------------------------
    cvFun0 <- cvFunS2[
        as.character(cvFunS2$Var1) == as.character(availStrat$Var1[indStr]) &
            as.character(cvFunS2$Var2) == as.character(availStrat$Var2[indStr]) &
            as.character(cvFunS2$Var3) == as.character(availStrat$Var3[indStr]),
        ,
        drop = FALSE
    ]

    if (nrow(cvFun0) == 0L || all(is.na(cvFun0$cv))) {
        return(empty_result)
    }

    ## ---------------------------------------------------------------------------
    ## Remove missing CV values and order by sample size
    ## ---------------------------------------------------------------------------
    cvFun <- cvFun0[!is.na(cvFun0$cv), , drop = FALSE]
    cvFunOrder <- cvFun[order(cvFun$sampSize), , drop = FALSE]

    if (nrow(cvFunOrder) == 0L) {
        return(empty_result)
    }

    ## ---------------------------------------------------------------------------
    ## Remove outliers before CV maximum
    ## ---------------------------------------------------------------------------
    cvMAXi <- which.max(cvFunOrder$cv)

    if (length(cvMAXi) == 1L && cvMAXi != 1L) {
        indOutliers <- unlist(lapply(seq_len(cvMAXi - 1L), function(i) {
            diffCV <- abs((cvFunOrder$cv[i] - cvFunOrder$cv[cvMAXi]) / cvFunOrder$cv[cvMAXi])
            if (diffCV > 0.5) i else NULL
        }))

        if (length(indOutliers) > 0L) {
            cvFunOrder <- cvFunOrder[-indOutliers, , drop = FALSE]
        }
    }

    if (nrow(cvFunOrder) == 0L) {
        return(empty_result)
    }

    ## ---------------------------------------------------------------------------
    ## Prepare data for density estimation
    ## ---------------------------------------------------------------------------
    dat <- data.frame(
        x = cvFunOrder$sampSize,
        y = cvFunOrder$cv,
        stringsAsFactors = FALSE
    )

    d <- try(stats::density(dat$y, bw = "sj"), silent = TRUE)

    min_accepted_sample_size <- as.numeric(
        as.character(
            table_min_accepted_sample_size$min_accepted_ss[
                as.character(table_min_accepted_sample_size$space) == as.character(availStrat$Var1[indStr]) &
                    as.character(table_min_accepted_sample_size$time) == as.character(availStrat$Var2[indStr]) &
                    as.character(table_min_accepted_sample_size$technical) == as.character(availStrat$Var3[indStr])
            ]
        )
    )

    if (inherits(d, "try-error") || length(min_accepted_sample_size) == 0L || is.na(min_accepted_sample_size[1])) {
        return(empty_result)
    }

    min_accepted_sample_size <- min_accepted_sample_size[1]

    ## ---------------------------------------------------------------------------
    ## Identify candidate solutions from local maxima
    ## ---------------------------------------------------------------------------
    loc_max_idx <- localMaxima(d$y)

    if (length(loc_max_idx) == 0L) {
        return(empty_result)
    }

    loc.max <- d$x[loc_max_idx] * max(dat$x)

    SOLs <- data.frame(
        Var1 = NA,
        Var2 = NA,
        Var3 = NA,
        solutions = c(round(loc.max), round(loc.max)[length(round(loc.max))]),
        maxCV = NA,
        minCV = NA,
        maxRR = NA,
        minRR = NA,
        meanRR = NA,
        noIterations = NA,
        stringsAsFactors = FALSE
    )

    for (rig in seq_len(nrow(SOLs))) {

        if (rig == 1) {
            dat_this_rig <- dat[dat$x <= as.numeric(SOLs$solutions[rig]), , drop = FALSE]
            cvFun0_this_rig <- cvFun0[cvFun0$sampSize <= as.numeric(SOLs$solutions[rig]), , drop = FALSE]
        } else if (rig == nrow(SOLs)) {
            dat_this_rig <- dat[dat$x > as.numeric(SOLs$solutions[rig - 1]), , drop = FALSE]
            cvFun0_this_rig <- cvFun0[cvFun0$sampSize > as.numeric(SOLs$solutions[rig - 1]), , drop = FALSE]
        } else {
            dat_this_rig <- dat[
                dat$x > as.numeric(SOLs$solutions[rig - 1]) &
                    dat$x <= as.numeric(SOLs$solutions[rig]),
                ,
                drop = FALSE
            ]
            cvFun0_this_rig <- cvFun0[
                cvFun0$sampSize > as.numeric(SOLs$solutions[rig - 1]) &
                    cvFun0$sampSize <= as.numeric(SOLs$solutions[rig]),
                ,
                drop = FALSE
            ]
        }

        if (nrow(dat_this_rig) > 0L) {
            SOLs$minCV[rig] <- round(min(dat_this_rig$y, na.rm = TRUE), 2)
            SOLs$maxCV[rig] <- round(max(dat_this_rig$y, na.rm = TRUE), 2)
        }

        if ("recyclingRate" %in% names(cvFun0_this_rig) && nrow(cvFun0_this_rig) > 0L) {
            SOLs$minRR[rig] <- min(cvFun0_this_rig$recyclingRate, na.rm = TRUE)
            SOLs$maxRR[rig] <- max(cvFun0_this_rig$recyclingRate, na.rm = TRUE)
            SOLs$meanRR[rig] <- round(mean(cvFun0_this_rig$recyclingRate, na.rm = TRUE), 2)
        }

        SOLs$noIterations[rig] <- nrow(cvFun0_this_rig)
    }

    accepted_solutions <- which(SOLs$solutions >= min_accepted_sample_size)

    ## ---------------------------------------------------------------------------
    ## Save density plot
    ## ---------------------------------------------------------------------------
    file_name <- file.path(
        output_dir,
        paste0(
"Optimization/","CV Density function ",
            cvFun0$Var1[1], "-", cvFun0$Var2[1], "-", cvFun0$Var3[1],
            ".png"
        )
    )

    grDevices::png(file = file_name, width = 21, height = 21, bg = "white", units = "cm", res = 200)

    plot(d, main = "Density plot of the CV values")
    graphics::abline(v = d$x[loc_max_idx])

    if (length(accepted_solutions) > 0L) {
        for (i in accepted_solutions[seq_len(min(length(accepted_solutions), 4))]) {
            try(graphics::abline(v = d$x[loc_max_idx][i], col = "red"), silent = TRUE)
        }
    }

    graphics::text(d$x[loc_max_idx], 0, round(loc.max))
    graphics::text(
        c(d$x[loc_max_idx], max(d$x)),
        c(d$y[loc_max_idx] + 0.1, min(d$y) + 0.1),
        paste(SOLs$maxCV, "-", SOLs$minCV),
        cex = 0.7
    )

    grDevices::dev.off()

    ## ---------------------------------------------------------------------------
    ## Final output
    ## ---------------------------------------------------------------------------
    SOLs$solutions[nrow(SOLs)] <- paste(">", SOLs$solutions[nrow(SOLs)], sep = "")

    if (length(accepted_solutions) == 0L) {
        result <- empty_result
    } else if (nrow(SOLs[accepted_solutions, , drop = FALSE]) <= 4) {
        result <- SOLs[accepted_solutions, , drop = FALSE]
    } else {
        result <- SOLs[accepted_solutions[1:4], , drop = FALSE]
    }

    result$Var1 <- availStrat$Var1[indStr]
    result$Var2 <- availStrat$Var2[indStr]
    result$Var3 <- availStrat$Var3[indStr]

    return(result)
}
