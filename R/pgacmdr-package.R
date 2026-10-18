#' @keywords internal
#'
#' @import shiny
#' @import shinydashboard
#' @importFrom Rcpp evalCpp
#' @importFrom magrittr %>%
#' @importFrom waiter useWaiter
#' @importFrom htmltools tags HTML tagList
#' @importFrom utils combn
#' @importFrom stats aggregate median sd
#' @importFrom utils head write.table
#' @useDynLib pgacmdr, .registration = TRUE
"_PACKAGE"

# Suppress R CMD check NOTEs for symbols that cannot be resolved statically.
utils::globalVariables(c(
  # data-masking column names used in dplyr verbs
  "Population", "Hs", "ID", "Observed_FIS", "Observed_FIT", "P_value",
  "Observed_FST", "Observed_HT", "Observed_HS",
  # remaining names flagged by R CMD check (column references / withTags)
  "setNames", "pop", "Allele", "Frequency", "Marker", "Na",
  "Significant", "CI_L", "CI_U", "median",
  "thead", "tr", "th", "table",
  "Ho", "He", "Locus", "N_Samples"
))
