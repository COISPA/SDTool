#' Compute bootstrap CV with fixed number of trips per sharing level
#'
#' Performs a bootstrap procedure within one stratum using a fixed number of
#' trips per sharing level. For each iteration, the function:
#' \enumerate{
#'   \item filters landing and sampling data to the target stratum;
#'   \item determines the pool of trips available in the stratum;
#'   \item samples trips according to the requested sharing structure;
#'   \item rebuilds a reduced COST catch sampling object;
#'   \item randomly selects one landing year within the filtered landings;
#'   \item computes the raised length distribution and its CV.
#' }
#'
#' @param caseInd Integer index identifying the row of \code{availStrat}
#'   corresponding to the target stratum.
#' @param stratification A \code{strIni()} object describing the stratification
#'   used in COST raising.
#' @param availStrat A data.frame containing at least the columns
#'   \code{Var1}, \code{Var2}, \code{Var3}.
#' @param availSamp A data.frame containing the available samples for the target
#'   species.
#' @param CSobject A \code{csData} or \code{csDataVal} object containing the
#'   catch sampling tables for the target species.
#' @param CLobjectList A COST landing object containing the landing table.
#' @param sppName Character string with the target species name.
#' @param nIter Number of bootstrap iterations.
#' @param loca_SH_levels A data.frame describing the sharing levels and number
#'   of trips to be sampled.
#' @param table_NO_TRIPS_4_BOOTSTRAP A data.frame containing the number of trips
#'   to be used for bootstrap by stratum. It must contain the columns
#'   \code{space}, \code{time}, \code{technical}, and \code{NO_TRIPS_4_BOOTSTRAP}.
#' @param cs_all A \code{csData} or \code{csDataVal} object containing all
#'   species, used to identify available trips in the stratum.
#' @param space Spatial stratification. Allowed values are \code{"GSA"},
#'   \code{"Country"}, and \code{"ALL"}.
#' @param time Temporal stratification. Allowed values are \code{"Q"},
#'   \code{"S"}, and \code{"Y"}.
#' @param technical Technical stratification. Allowed values are
#'   \code{"lev6"}, \code{"lev4"}, and \code{"NONE"}.
#' @param sharing_level Optional character string defining the sharing level.
#'   Allowed values are \code{"space"}, \code{"time"}, \code{"technical"}, or
#'   \code{NULL}.
#' @param SD_var1 Spatial variable used in the sampling design.
#' @param SD_var2 Temporal variable used in the sampling design.
#' @param SD_var3 Technical variable used in the sampling design.
#' @param current_species Integer index of the current species in
#'   \code{target_SPECIES}.
#' @param target_SPECIES Character vector of target species names.
#' @param trips_bootstrap_first_species Optional nested list used to store
#'   trip selections for the first species and reuse them for the others.
#' @param LengthStr_dataframe Optional data.frame used to accumulate raised
#'   length structures across iterations.
#' @param CVbyLength_dataframe Optional data.frame used to accumulate CV by
#'   length class across iterations.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{res_all_iter}{A data.frame with one row per iteration and columns
#'   \code{cv}, \code{nMeas}, and \code{sampSize}.}
#'   \item{LengthStr_dataframe}{Updated length-structure data.frame.}
#'   \item{CVbyLength_dataframe}{Updated CV-by-length data.frame.}
#'   \item{trips_bootstrap_first_species}{Updated list of bootstrap trip
#'   selections.}
#' }
#'
#' @keywords internal
#' @export
cvFunStr_FixedNoTrips <- function(caseInd,
                                  stratification,
                                  availStrat,
                                  availSamp,
                                  CSobject,
                                  CLobjectList,
                                  sppName,
                                  nIter,
                                  loca_SH_levels,
                                  table_NO_TRIPS_4_BOOTSTRAP,
                                  cs_all,
                                  space = c("GSA", "Country", "ALL"),
                                  time = c("Q", "S", "Y"),
                                  technical = c("lev6", "lev4", "NONE"),
                                  sharing_level = NULL,
                                  SD_var1 = NULL,
                                  SD_var2 = NULL,
                                  SD_var3 = NULL,
                                  current_species = 1L,
                                  target_SPECIES = sppName,
                                  trips_bootstrap_first_species = list(),
                                  LengthStr_dataframe = data.frame(),
                                  CVbyLength_dataframe = data.frame()) {

    if (!is.data.frame(availStrat)) {
        stop("`availStrat` must be a data.frame.")
    }
    if (!all(c("Var1", "Var2", "Var3") %in% names(availStrat))) {
        stop("`availStrat` must contain columns `Var1`, `Var2`, and `Var3`.")
    }
    if (!is.data.frame(availSamp)) {
        stop("`availSamp` must be a data.frame.")
    }
    if (!is.data.frame(loca_SH_levels)) {
        stop("`loca_SH_levels` must be a data.frame.")
    }
    if (!is.data.frame(table_NO_TRIPS_4_BOOTSTRAP)) {
        stop("`table_NO_TRIPS_4_BOOTSTRAP` must be a data.frame.")
    }
    if (!inherits(CSobject, c("csData", "csDataVal"))) {
        stop("`CSobject` must be a 'csData' or 'csDataVal' object.")
    }
    if (!inherits(cs_all, c("csData", "csDataVal"))) {
        stop("`cs_all` must be a 'csData' or 'csDataVal' object.")
    }
    if (length(caseInd) != 1L || caseInd < 1L || caseInd > nrow(availStrat)) {
        stop("`caseInd` must be a valid row index of `availStrat`.")
    }

    space <- match.arg(space)
    time <- match.arg(time)
    technical <- match.arg(technical)

    to_date <- function(x) {
        if (inherits(x, "Date")) return(x)
        as.Date(x, format = "%d/%m/%Y")
    }

    get_quarter <- function(x) {
        x <- to_date(x)
        ((as.integer(format(x, "%m")) - 1L) %/% 3L) + 1L
    }

    get_semester <- function(x) {
        ifelse(get_quarter(x) %in% c(1, 2), 1L, 2L)
    }

    case <- availStrat[caseInd, 1:3, drop = FALSE]

    max_no_trips_boot <- table_NO_TRIPS_4_BOOTSTRAP$NO_TRIPS_4_BOOTSTRAP[
        as.character(table_NO_TRIPS_4_BOOTSTRAP$space) == as.character(case$Var1[1]) &
            as.character(table_NO_TRIPS_4_BOOTSTRAP$time) == as.character(case$Var2[1]) &
            as.character(table_NO_TRIPS_4_BOOTSTRAP$technical) == as.character(case$Var3[1])
    ]

    if (length(max_no_trips_boot) == 0L || is.na(max_no_trips_boot[1])) {
        max_no_trips_boot <- 0
    } else {
        max_no_trips_boot <- as.numeric(max_no_trips_boot[1])
    }

    print(
        paste(
            "MAX_NO_TRIPS to be used in the BOOTSTRAP for sharing in strata",
            case$Var1[1], "-", case$Var2[1], "-", case$Var3[1], ":", max_no_trips_boot
        ),
        quote = FALSE
    )

    CLsub <- CLobjectList

    if (space == "GSA") {
        CLsub <- subset(CLsub, area == as.character(case$Var1[1]), link = TRUE)
    }
    if (space == "Country") {
        CLsub <- subset(CLsub, vslFlgCtry == as.character(case$Var1[1]), link = TRUE)
    }

    if (time == "Q") {
        CLsub <- subset(CLsub, quarter == as.character(case$Var2[1]))
    }
    if (time == "S") {
        if (as.character(case$Var2[1]) == "1") {
            CLsub <- subset(CLsub, quarter == 1 | quarter == 2)
        } else {
            CLsub <- subset(CLsub, quarter == 3 | quarter == 4)
        }
    }

    if (technical == "lev6") {
        CLsub <- subset(CLsub, foCatEu6 == as.character(case$Var3[1]))
    }
    if (technical == "lev4") {
        CLsub <- subset(CLsub, foCatEu5 == as.character(case$Var3[1]))
    }

    availSamp_all_species <- merge(cs_all@hh, cs_all@hl)
    availSamp_all_species$SampQuarter <- get_quarter(availSamp_all_species$date)
    availSamp_all_species$SampSemester <- get_semester(availSamp_all_species$date)

    availSamp_all_species <- setVar123(
        x = availSamp_all_species,
        space = space,
        time = time,
        technical = technical,
        is_landings = FALSE
    )

    availSamp_all_species <- availSamp_all_species[
        as.character(availSamp_all_species$Var1) == as.character(case$Var1[1]) &
            as.character(availSamp_all_species$Var2) == as.character(case$Var2[1]) &
            as.character(availSamp_all_species$Var3) == as.character(case$Var3[1]),
        ,
        drop = FALSE
    ]

    sampPoolStr <- unique(availSamp_all_species$trpCode)

    print("**********************************************************************", quote = FALSE)
    print(paste("total number of trips available for the bootstrap:", length(sampPoolStr)), quote = FALSE)
    print("**********************************************************************", quote = FALSE)

    if (!is.null(sharing_level) && length(sharing_level) > 0L) {
        sharing_level <- match.arg(sharing_level, c("space", "time", "technical"))

        if (sharing_level == "space") {
            SH_levels_Str <- loca_SH_levels[
                as.character(loca_SH_levels[, 2]) == as.character(case$Var1[1]),
                ,
                drop = FALSE
            ]
        } else if (sharing_level == "time") {
            SH_levels_Str <- loca_SH_levels[
                as.character(loca_SH_levels[, 1]) == as.character(case$Var1[1]) &
                    as.character(loca_SH_levels[, 4]) == as.character(case$Var3[1]),
                ,
                drop = FALSE
            ]
        } else {
            if (as.character(case$Var3[1]) != "ALL") {
                SH_levels_Str <- loca_SH_levels[
                    as.character(loca_SH_levels[, 1]) == as.character(case$Var1[1]) &
                        as.character(loca_SH_levels[, 4]) == as.character(case$Var3[1]),
                    ,
                    drop = FALSE
                ]
            } else {
                SH_levels_Str <- loca_SH_levels
            }
        }
    } else {
        SH_levels_Str <- loca_SH_levels[
            as.character(loca_SH_levels$space_original) == as.character(case$Var1[1]) &
                as.character(loca_SH_levels$time_original) == as.character(case$Var2[1]) &
                as.character(loca_SH_levels$technical_original) == as.character(case$Var3[1]),
            ,
            drop = FALSE
        ]
    }

    hl_this_stratum <- CSobject@hl[as.character(CSobject@hl$trpCode) %in% sampPoolStr, , drop = FALSE]

    if (nrow(hl_this_stratum) > 0L) {
        hl2 <- data.frame(
            Var1 = rep(as.character(case$Var1[1]), nrow(hl_this_stratum)),
            Var2 = rep(as.character(case$Var2[1]), nrow(hl_this_stratum)),
            Var3 = rep(as.character(case$Var3[1]), nrow(hl_this_stratum)),
            hl_this_stratum,
            stringsAsFactors = FALSE
        )
        LengthStr_dataframe <- rbind(LengthStr_dataframe, hl2)
    }

    res_all_iter <- data.frame(cv = numeric(0), nMeas = numeric(0), sampSize = numeric(0))
    LANDINGS_EXISTS <- TRUE

    key_case <- paste(case$Var1[1], case$Var2[1], case$Var3[1], sep = "_")
    if (is.null(trips_bootstrap_first_species[[key_case]])) {
        trips_bootstrap_first_species[[key_case]] <- vector("list", nIter)
    }

    for (n_it in seq_len(nIter)) {

        if (!LANDINGS_EXISTS) {
            res_all_iter <- rbind(
                res_all_iter,
                data.frame(cv = NA_real_, nMeas = NA_real_, sampSize = NA_real_)
            )
            next
        }

        loca_SH_levels_Str <- SH_levels_Str

        if (!is.null(sharing_level) && length(sharing_level) > 0L) {
            if (sharing_level == "space") {
                if (SD_var1 == "GSA") {
                    trips_by_level <- availSamp_all_species[, c("trpCode", "area"), drop = FALSE]
                } else if (SD_var1 == "Country") {
                    trips_by_level <- availSamp_all_species[, c("vslFlgCtry", "trpCode"), drop = FALSE]
                } else {
                    trips_by_level <- availSamp_all_species[, c("area", "trpCode"), drop = FALSE]
                    trips_by_level[, 1] <- "ALL"
                }
            } else if (sharing_level == "time") {
                if (SD_var2 == "Q") {
                    trips_by_level <- availSamp_all_species[, c("trpCode", "SampQuarter"), drop = FALSE]
                } else if (SD_var2 == "S") {
                    trips_by_level <- availSamp_all_species[, c("trpCode", "SampSemester"), drop = FALSE]
                } else {
                    trips_by_level <- availSamp_all_species[, c("year", "trpCode"), drop = FALSE]
                }
            } else {
                if (SD_var3 == "lev4") {
                    trips_by_level <- availSamp_all_species[, c("trpCode", "foCatEu5"), drop = FALSE]
                } else {
                    trips_by_level <- availSamp_all_species[, c("trpCode", "foCatEu6"), drop = FALSE]
                }
            }
        } else {
            trips_by_level <- availSamp_all_species[, c("year", "trpCode"), drop = FALSE]
            loca_SH_levels_Str$time_original <- 2999
            loca_SH_levels_Str$time <- 2999
        }

        trips_by_level_selection <- unique(trips_by_level[trips_by_level$trpCode %in% sampPoolStr, , drop = FALSE])
        selected_samples <- character(0)

        if (current_species == 1L || is.null(trips_bootstrap_first_species[[key_case]][[n_it]])) {

            lengths_sharing <- numeric(0)

            for (ll in seq_len(nrow(loca_SH_levels_Str))) {

                if (nrow(loca_SH_levels_Str) > 1L) {
                    if (!is.null(sharing_level) && sharing_level == "space") {
                        sampPoolStr_this_level <- trips_by_level_selection[
                            as.character(trips_by_level_selection[, 2]) == as.character(loca_SH_levels_Str[ll, 1]),
                            1
                        ]
                    } else if (!is.null(sharing_level) && sharing_level == "technical") {
                        sampPoolStr_this_level <- trips_by_level_selection[
                            as.character(trips_by_level_selection[, 2]) == as.character(loca_SH_levels_Str[ll, 3]),
                            1
                        ]
                    } else {
                        sampPoolStr_this_level <- trips_by_level_selection[
                            as.character(trips_by_level_selection[, 2]) == as.character(loca_SH_levels_Str[ll, 2]),
                            1
                        ]
                    }
                } else {
                    sampPoolStr_this_level <- trips_by_level_selection[, 1]
                }

                if (length(sampPoolStr_this_level) > 0L) {
                    if (nrow(loca_SH_levels_Str) > 1L && max_no_trips_boot != 0) {
                        trips_to_be_used <- sample(sampPoolStr_this_level, max_no_trips_boot, replace = TRUE)
                    } else {
                        trips_to_be_used <- sampPoolStr_this_level
                    }

                    SH_level_sel <- sample(trips_to_be_used, loca_SH_levels_Str$no_of_trips[ll], replace = TRUE)
                    selected_samples <- c(selected_samples, SH_level_sel)
                    lengths_sharing <- c(lengths_sharing, length(SH_level_sel))
                } else {
                    print("ATTENTION: 0 selected!")
                    lengths_sharing <- c(lengths_sharing, 0)
                }
            }

            print(paste("No. of samples selected for bootstrap iteration n°", n_it), quote = FALSE)
            print(length(selected_samples), quote = FALSE)
            cat(lengths_sharing, sep = " - ")
            cat("\n")

            trips_bootstrap_first_species[[key_case]][[n_it]] <- selected_samples
        } else {
            selected_samples <- trips_bootstrap_first_species[[key_case]][[n_it]]
            print(paste("No. of samples selected for bootstrap iteration n°", n_it), quote = FALSE)
            print(length(selected_samples), quote = FALSE)
        }

        sampSize <- length(selected_samples)
        extentionList <- seq_len(sampSize)

        selectTR <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = CSobject@tr, extention.list = extentionList, selected.tripcode = selected_samples
        ))
        selectHH <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = CSobject@hh, extention.list = extentionList, selected.tripcode = selected_samples
        ))
        selectSL <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = CSobject@sl, extention.list = extentionList, selected.tripcode = selected_samples
        ))
        selectHL <- do.call("rbind", lapply(seq_len(sampSize), select.and.rename.tripcode2,
                                            table = CSobject@hl, extention.list = extentionList, selected.tripcode = selected_samples
        ))

        selectSL <- selectSL[selectSL$catchCat != "DIS", , drop = FALSE]
        selectHL <- selectHL[selectHL$catchCat != "DIS", , drop = FALSE]

        no_trip_positive_to_species <- length(unique(selectSL$trpCode))

        if (no_trip_positive_to_species <= 1L) {
            if (no_trip_positive_to_species == 0L) {
                print("Impossible to calculate CV due to only samples of discard in the stratum!", quote = FALSE)
            } else {
                print("Impossible to calculate CV due to the number of trips positive to the species in the stratum (equal to 1)!", quote = FALSE)
            }

            res_all_iter <- rbind(
                res_all_iter,
                data.frame(cv = NA_real_, nMeas = NA_real_, sampSize = sampSize)
            )
            next
        }

        subCS <- csData(tr = selectTR, hh = selectHH, sl = selectSL, hl = selectHL)

        years_for_selection <- unique(CLsub@cl$year)
        if (length(years_for_selection) == 0L) {
            LANDINGS_EXISTS <- FALSE
            res_all_iter <- rbind(
                res_all_iter,
                data.frame(cv = NA_real_, nMeas = NA_real_, sampSize = sampSize)
            )
            next
        }

        CL_ye <- if (length(years_for_selection) == 1L) years_for_selection else sample(years_for_selection, 1)
        CLobject_ye <- subset(CLsub, year == CL_ye, link = TRUE)

        if (nrow(CLobject_ye@cl) == 0L) {
            LANDINGS_EXISTS <- FALSE
            res_all_iter <- rbind(
                res_all_iter,
                data.frame(cv = NA_real_, nMeas = NA_real_, sampSize = sampSize)
            )
            next
        }

        CLobject_ye@cl$year <- 2999

        nMeas <- sum(subCS@hl$lenNum, na.rm = TRUE)

        CSval <- csDataVal(subCS)
        CLval <- clDataVal(CLobject_ye)
        CScons <- csDataCons(CSval, stratification)
        CLcons <- clDataCons(CLval, stratification)

        dbeOutput <- dbeObject(
            desc = sppName,
            species = sppName,
            catchCat = "Lan",
            strataDesc = stratification,
            methodDesc = "analytical"
        )

        raising <- try(
            RaiseLgth(dbeOutput, CScons, CLcons, spp = sppName, taxon = sppName, strDesc = stratification),
            silent = TRUE
        )

        if (inherits(raising, "try-error")) {
            print("Impossible to calculate CV due a problem in COST function!", quote = FALSE)
            res_all_iter <- rbind(
                res_all_iter,
                data.frame(cv = NA_real_, nMeas = NA_real_, sampSize = sampSize)
            )
            next
        }

        cv <- raising@lenNum$DCRcvIndicator

        if (length(raising@lenStruc) > 0L && nrow(raising@lenStruc[[1]]) > 0L) {
            lenstructure <- raising@lenStruc[[1]]
            result2 <- data.frame(
                species = rep(target_SPECIES[current_species], nrow(lenstructure)),
                iteration = rep(n_it, nrow(lenstructure)),
                Var1 = rep(as.character(case$Var1[1]), nrow(lenstructure)),
                Var2 = rep(as.character(case$Var2[1]), nrow(lenstructure)),
                Var3 = rep(as.character(case$Var3[1]), nrow(lenstructure)),
                lenstructure[, 4:5],
                stringsAsFactors = FALSE
            )
            colnames(result2) <- c("species", "iteration", "Var1", "Var2", "Var3", "length", "value")
            LengthStr_dataframe <- rbind(LengthStr_dataframe, result2)
        }

        if (!is.null(raising@lenNum$cv) && nrow(raising@lenNum$cv) > 0L) {
            lenCV <- raising@lenNum$cv
            result3 <- data.frame(
                species = rep(target_SPECIES[current_species], nrow(lenCV)),
                iteration = rep(n_it, nrow(lenCV)),
                Var1 = rep(as.character(case$Var1[1]), nrow(lenCV)),
                Var2 = rep(as.character(case$Var2[1]), nrow(lenCV)),
                Var3 = rep(as.character(case$Var3[1]), nrow(lenCV)),
                lenCV[, 4:5],
                stringsAsFactors = FALSE
            )
            colnames(result3) <- c("species", "iteration", "Var1", "Var2", "Var3", "length", "CV")
            CVbyLength_dataframe <- rbind(CVbyLength_dataframe, result3)
        }

        res_all_iter <- rbind(
            res_all_iter,
            data.frame(cv = cv, nMeas = nMeas, sampSize = sampSize)
        )
    }

    if (!LANDINGS_EXISTS) {
        print("Impossible to calculate CV on the raised LFD for all the iterations due to MISSING LANDINGS in the stratum!", quote = FALSE)
    } else {
        print(paste(nIter, "iterations done!"), quote = FALSE)
    }

    list(
        res_all_iter = res_all_iter,
        LengthStr_dataframe = LengthStr_dataframe,
        CVbyLength_dataframe = CVbyLength_dataframe,
        trips_bootstrap_first_species = trips_bootstrap_first_species
    )
}
