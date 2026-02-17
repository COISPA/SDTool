#' Dummy CL dataset for GSA99 (HKE)
#'
#' A small dummy dataset derived from CL-GSA99 2014-2024_HKE.csv.
#' Used for examples and testing of `convert_cl_to_costCL()`.
#'
#' @format A data.frame with 100 rows and the following variables:
#' \describe{
#'   \item{vslFlgCtry}{Vessel flag country}
#'   \item{year}{Year}
#'   \item{quarter}{Quarter}
#'   \item{month}{Month}
#'   \item{area}{GSA area}
#'   \item{taxon}{Species}
#'   \item{foCatNat}{National fishing activity category}
#'   \item{foCatEu6}{EU fishing activity category level 6}
#'   \item{harbour}{Harbour name}
#'   \item{landWt}{Landing weight}
#'   \item{landValue}{Landing value}
#' }
#'
#' @source Simulated subset of CL-GSA99 2014-2024_HKE.csv
"cl_gsa99_dummy"
