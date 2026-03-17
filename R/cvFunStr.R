#' Compute bootstrap CV results for one stratum
#'
#' Performs a bootstrap procedure within one stratum to estimate the
#' coefficient of variation (CV) of raised length frequencies for a target
#' species.
#'
#' The function:
#' \enumerate{
#'   \item filters the landings object and the available sample pool to the
#'   selected stratum;
#'   \item repeatedly samples trip codes with replacement;
#'   \item rebuilds a reduced COST catch sampling object;
#'   \item randomly selects one landing year within the filtered landings;
#'   \item computes the raised length distribution and associated CV;
#'   \item stores iteration-level results.
#' }
#'
#' @param caseInd Integer index identifying the row of \code{availStrat}
#'   corresponding to the target stratum.
#' @param stratification A \code{strIni()} object describing the stratification
#'   used in COST raising.
#' @param availStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, \code{Var3}, and \code{maxSampSize}.
#' @param availSamp A data.frame containing the available samples. It must
#'   include at least \code{trpCode}, \code{date}, and the columns required
#'   by the selected stratification.
#' @param CSobject A \code{csData} or \code{csDataVal} object containing the
#'   catch sampling tables.
#' @param CLobjectList A COST landing object containing the landing table.
#' @param sppName Character string with the target species name.
#' @param nIter Number of bootstrap iterations.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#' @param lengthStr_dataframe Optional data.frame used to accumulate length
#'   structure results across iterations. Default is \code{NULL}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{res_all_iter}{A data.frame with one row per iteration and columns
#'     \code{cv}, \code{nMeas}, \code{sampSize}, and \code{recyclingRate}.}
#'   \item{LengthStr_dataframe}{A data.frame containing the sampled length
#'     structures, if available.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- cvFunStr(
#'   caseInd = 1,
#'   stratification = COSTcore::strIni(),
#'   availStrat = availStrat,
#'   availSamp = availSamp,
#'   CSobject = CS_boot,
#'   CLobjectList = CL_boot,
#'   sppName = "Merluccius merluccius",
#'   nIter = 100,
#'   space = "GSA",
#'   time = "Q",
#'   technical = "lev6"
#' )
#' }
#'
#' @keywords internal
#' @export
#' @import COSTcore COSTdbe

cvFunStr <- function(caseInd,
                     stratification,
                     availStrat,
                     availSamp,
                     CSobject,
                     CLobjectList,
                     sppName,
                     nIter,
                     space = c("GSA", "Country", "ALL"),
                     time = c("Q", "S", "Y"),
                     technical = c("lev6", "lev4", "NONE"),
                     lengthStr_dataframe = NULL) {


    if (FALSE){
        caseInd=nr
        stratification=COSTcore::strIni()
        availStrat
        availSamp
        CSobject=CS_boot
        CLobjectList=CL_boot
        target_SPECIES
        nIter=5
        space
        time
        technical
        lengthStr_dataframe =LengthStr_dataframe
    }
    ## ---------------------------------------------------------------------------
    ## Checks
    ## ---------------------------------------------------------------------------
    if (!is.data.frame(availStrat)) {
        stop("`availStrat` must be a data.frame.")
    }

    if (!all(c("Var1", "Var2", "Var3", "maxSampSize") %in% names(availStrat))) {
        stop("`availStrat` must contain columns `Var1`, `Var2`, `Var3`, and `maxSampSize`.")
    }

    if (!is.data.frame(availSamp)) {
        stop("`availSamp` must be a data.frame.")
    }

    if (!inherits(CSobject, c("csData", "csDataVal"))) {
        stop("`CSobject` must be a 'csData' or 'csDataVal' object.")
    }

    if (!inherits(CLobjectList, c("clData", "clDataVal", "csData", "csDataVal"))) {
        stop("`CLobjectList` must be a valid COST landing object.")
    }

    if (length(caseInd) != 1L || !is.numeric(caseInd) || is.na(caseInd)) {
        stop("`caseInd` must be a single non-missing numeric index.")
    }

    caseInd <- as.integer(caseInd)

    if (caseInd < 1L || caseInd > nrow(availStrat)) {
        stop("`caseInd` is out of bounds for `availStrat`.")
    }

    if (!is.numeric(nIter) || length(nIter) != 1L || is.na(nIter) || nIter < 1) {
        stop("`nIter` must be a positive integer.")
    }

    nIter <- as.integer(nIter)

    if (!is.character(sppName) || length(sppName) != 1L || !nzchar(sppName)) {
        stop("`sppName` must be a single non-empty character string.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    if (is.null(lengthStr_dataframe)) {
        lengthStr_dataframe <- data.frame()
    }

    ## ---------------------------------------------------------------------------
    ## Helpers
    ## ---------------------------------------------------------------------------
    # to_date <- function(x) {
    #     if (inherits(x, "Date")) {
    #         return(x)
    #     }
    #     as.Date(x, format = "%d/%m/%Y")
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
    ## Selected stratum
    ## ---------------------------------------------------------------------------
    case <- availStrat[caseInd, c("Var1", "Var2", "Var3"), drop = FALSE]

    print("______________________________________________________________________", quote = FALSE)
    print(
        paste(
            "Stratum:",
            toupper(as.character(case$Var1[1])),
            toupper(as.character(case$Var2[1])),
            toupper(as.character(case$Var3[1]))
        ),
        quote = FALSE
    )
    print("______________________________________________________________________", quote = FALSE)

    ## ---------------------------------------------------------------------------
    ## Subset CL object to the selected stratum
    ## ---------------------------------------------------------------------------
    CLsub <- CLobjectList

    if (space == "GSA") {
        CLsub <- COSTcore::subset(CLsub, area == as.character(case$Var1[1]), link = TRUE)
    }

    if (space == "Country") {
        CLsub <-  COSTcore::subset(CLsub, vslFlgCtry == as.character(case$Var1[1]), link = TRUE)
    }

    if (time == "Q") {
        CLsub <-  COSTcore::subset(CLsub, quarter == as.character(case$Var2[1]))
    }

    if (time == "S") {
        if (as.character(case$Var2[1]) == "1") {
            CLsub <-  COSTcore::subset(CLsub, quarter == 1 | quarter == 2)
        } else {
            CLsub <- subset(CLsub, quarter == 3 | quarter == 4)
        }
    }

    if (technical == "lev6") {
        CLsub <-  COSTcore::subset(CLsub, foCatEu6 == as.character(case$Var3[1]))
    }

    if (technical == "lev4") {
        CLsub <-  COSTcore::subset(CLsub, foCatEu5 == as.character(case$Var3[1]))
    }

    ## ---------------------------------------------------------------------------
    ## Build pool of available sample trip codes in the selected stratum
    ## ---------------------------------------------------------------------------
    if (space == "GSA") {
        sampPoolStr1 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$area == case$Var1]),
            group1 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (space == "Country") {
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

    if (time == "Q") {
        sampPoolStr2 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[get_quarter(availSamp$date) %in% case$Var2]),
            group2 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (time == "S") {
        sampPoolStr2 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[get_semester(availSamp$date) %in% case$Var2]),
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

    if (technical == "lev6") {
        sampPoolStr3 <- data.frame(
            sampPoolStr = unique(availSamp$trpCode[availSamp$foCatEu6 == case$Var3]),
            group3 = 1,
            stringsAsFactors = FALSE
        )
    }

    if (technical == "lev4") {
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

    sampPoolStr <- merge(
        sampPoolStr1,
        merge(sampPoolStr2, sampPoolStr3, all = TRUE),
        all = TRUE
    )

    sampPoolStr$sum_groups <- rowSums(sampPoolStr[, c("group1", "group2", "group3")], na.rm = TRUE)
    sampPoolStr <- sampPoolStr[!is.na(sampPoolStr$sum_groups), , drop = FALSE]
    sampPoolStr <- as.character(sampPoolStr$sampPoolStr[sampPoolStr$sum_groups == 3])

    print("**********************************************************************", quote = FALSE)
    print(paste("total number of trips available for the bootstrap:", length(sampPoolStr)), quote = FALSE)
    print("**********************************************************************", quote = FALSE)

    ## ---------------------------------------------------------------------------
    ## Bootstrap settings
    ## ---------------------------------------------------------------------------
    max.sample.size <- availStrat$maxSampSize[caseInd]
    min.sample.size <- 2

    res_all_iter <- data.frame(
        cv = numeric(0),
        nMeas = numeric(0),
        sampSize = numeric(0),
        recyclingRate = numeric(0)
    )

    LANDINGS_EXISTS <- TRUE

    ## ---------------------------------------------------------------------------
    ## Iterations
    ## ---------------------------------------------------------------------------
    for (n_it in seq_len(nIter)) {

        if (!LANDINGS_EXISTS) {
            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = NA_real_,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        if (length(sampPoolStr) == 0L || is.na(max.sample.size) || max.sample.size < min.sample.size) {
            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = NA_real_,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        sampSize <- sample(min.sample.size:max.sample.size, 1, replace = FALSE)
        sampSize <- round(sampSize, 0)

        selectSamp <- sample(sampPoolStr, sampSize, replace = TRUE)
        rate <- (length(selectSamp) - length(unique(selectSamp))) / length(selectSamp)

        print(paste("No. of samples selected for bootstrap iteration n°", n_it), quote = FALSE)
        print(length(selectSamp), quote = FALSE)

        extentionList <- seq_along(selectSamp)

        tr <- CSobject@tr
        hh <- CSobject@hh
        sl <- CSobject@sl
        hl <- CSobject@hl

        selectTR <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = tr, extention.list = extentionList,
                                            selected.tripcode = selectSamp))
        selectHH <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = hh, extention.list = extentionList,
                                            selected.tripcode = selectSamp))
        selectSL <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = sl, extention.list = extentionList,
                                            selected.tripcode = selectSamp))
        selectHL <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = hl, extention.list = extentionList,
                                            selected.tripcode = selectSamp))

        selectSL <- selectSL[selectSL$catchCat != "DIS", , drop = FALSE]
        selectHL <- selectHL[selectHL$catchCat != "DIS", , drop = FALSE]

        no_trip_positive_to_species <- length(unique(selectSL$trpCode))

        if (no_trip_positive_to_species <= 1) {
            if (no_trip_positive_to_species == 0) {
                print("Impossible to calculate CV due to only samples of discard in the stratum!", quote = FALSE)
            } else {
                print("Impossible to calculate CV due to the number of trips positive to the species in the stratum (equal to 1)!", quote = FALSE)
            }

            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = sampSize,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        subCS <- COSTcore::csData(
            tr = selectTR,
            hh = selectHH,
            sl = selectSL,
            hl = selectHL
        )

        years_for_selection <- unique(CLsub@cl$year)

        if (length(years_for_selection) == 0L) {
            LANDINGS_EXISTS <- FALSE
            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = sampSize,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        CL_ye <- if (length(years_for_selection) == 1L) {
            years_for_selection
        } else {
            sample(years_for_selection, 1)
        }

        CLobject_ye <- COSTcore::subset(CLsub, year == CL_ye, link = TRUE)

        if (nrow(CLobject_ye@cl) == 0L) {
            LANDINGS_EXISTS <- FALSE
            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = sampSize,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        CLobject_ye@cl$year <- 2999

        nMeas <- sum(subCS@hl$lenNum, na.rm = TRUE)

        CSval <- COSTcore::csDataVal(subCS)
        CLval <- COSTcore::clDataVal(CLobject_ye)

        CScons <- COSTcore::csDataCons(CSval, stratification)
        CLcons <- COSTcore::clDataCons(CLval, stratification)

        dbeOutput <- COSTdbe::dbeObject(
            desc = sppName,
            species = sppName,
            catchCat = "Lan",
            strataDesc = stratification,
            methodDesc = "analytical"
        )

        raising <- try(
            COSTdbe::RaiseLgth(
                dbeOutput,
                CScons,
                CLcons,
                spp = sppName,
                taxon = sppName,
                strDesc = stratification
            ),
            silent = TRUE
        )

        if (inherits(raising, "try-error")) {
            print("Impossible to calculate CV due a problem in COST function!", quote = FALSE)

            result <- data.frame(
                cv = NA_real_,
                nMeas = NA_real_,
                sampSize = sampSize,
                recyclingRate = NA_real_
            )
            res_all_iter <- rbind(res_all_iter, result)
            next
        }

        cv <- raising@lenNum$DCRcvIndicator
        rr <- round(rate, 2)

        if (length(raising@lenStruc) > 0 && nrow(raising@lenStruc[[1]]) > 0) {
            lenstructure <- raising@lenStruc[[1]]

            result2 <- data.frame(
                iteration = rep(n_it, nrow(lenstructure)),
                Var1 = rep(as.character(case$Var1[1]), nrow(lenstructure)),
                Var2 = rep(as.character(case$Var2[1]), nrow(lenstructure)),
                Var3 = rep(as.character(case$Var3[1]), nrow(lenstructure)),
                lenstructure[, 4:5],
                stringsAsFactors = FALSE
            )

            lengthStr_dataframe <- rbind(lengthStr_dataframe, result2)
        }

        result <- data.frame(
            cv = cv,
            nMeas = nMeas,
            sampSize = sampSize,
            recyclingRate = rr
        )

        res_all_iter <- rbind(res_all_iter, result)
    }

    if (!LANDINGS_EXISTS) {
        print("Impossible to calculate CV on the raised LFD for all the iterations due to MISSING LANDINGS in the stratum!", quote = FALSE)
    } else {
        print(paste(nIter, "iterations done!"), quote = FALSE)
    }

    return(list(
        res_all_iter = res_all_iter,
        LengthStr_dataframe = lengthStr_dataframe
    ))
}
