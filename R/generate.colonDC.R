

#' Colon cancer data
#'
#' A dataset containing the individual baseline and follow-up data on more than 15,000
#' colon cancer patient. The data is a data cleaned version of the \code{rstpm2::colon} data set.
#'
#' \itemize{
#'   \item sex. Either female or male
#'   \item age. Diagnostic age
#'   \item stage. Clinical stage at diagnosis (either unknown, localised, regional, or distant)
#'   \item statusDC. Alive indicator with cause of death
#'   \item subsite. Anatomical subsite of tumour (either coecum and ascending, transverse,
#'   descending and sigmoid, or other and NOS)
#'   \item dx. Date of diagnosis
#'   \item exit. Date of study exit
#'   \item status. Alive indicator (0 = alive, 1 = dead)
#'   \item FU. Follow-up time measured in days
#'   \item FUyear. Follow-up time measured in years
#'   \item agedays. Diagnostic age in days.
#' }
#'
#' @name colonDC
#' @usage data(colonDC)
#' @format A data frame with 15564 rows and 11 variables
"colonDC"
