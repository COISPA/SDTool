#' Data investigation on COST tables (CS and CL)
#'

#' @param CS A data.frame containing the RCG Catch Sampling table.
#' @param CL A data.frame containing the RCG landing table.
#' @param metier_shortcode data frame reporting the short codes for the metier considered.
#' @param space	space stratification (allowed values: GSA, Country,ALL)
#' @param time	time stratification (allowed values: Y, Q,S)
#' @param technical	gear/metier stratification (allowed values: lev4, lev6,NONE)
#' @param target_SPECIES target species on which apply the sampling optimization algorithm
#' @return graphs and tables
#' @export
#'
#' @examples

investigateData <- function(
        CS, CL,metier_shortcode=metier_shortcode_default, space="ALL",time="Y", technical="lev4",target_SPECIES){

}
