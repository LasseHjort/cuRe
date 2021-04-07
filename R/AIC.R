#' Akaike information criterion
#'
#' This function computes the AIC for a fitted parametric cure model.
#'
#' @param object Object of class cm or gfcm.
#' @param k Number to control if either AIC or BIC is to be computed (default is 2 equal to AIC).
#' @export
#' @method AIC cm
AIC.cm <- function(object, k = 2){
  npar <- length(object$optim$par)
  lnL <- object$ML
  aic <- as.vector(-2*lnL + k * npar)

  return(aic)
}

#' Akaike information criterion
#'
#' This function computes the AIC for a fitted parametric cure model.
#'
#' @param object Object of class cm or gfcm.
#' @param k Number to control if either AIC or BIC is to be computed (default is 2 equal to AIC).
#' @export
#' @method AIC gfcm
AIC.gfcm <- function(object, k = 2){
  npar <- length(object$coefs) + length(object$coefs.spline)
  lnL <- -object$NegMaxLik
  aic <- as.vector(-2*lnL + k * npar)

  return(aic)
}





#' Bayesian information criterion
#'
#' This function computes the BIC for a fitted parametric cure model.
#'
#' @param object Object of class cm or gfcm.
#' @export
#' @method BIC cm
BIC.cm <- function(object, ...)
{
  bic <- AIC(object, k = log(nobs(object)))

  return(bic)
}

#' Bayesian information criterion
#'
#' This function computes the BIC for a fitted parametric cure model.
#'
#' @param object Object of class cm or gfcm.
#' @export
#' @method BIC gfcm
BIC.gfcm <- function(object, ...)
{
  bic <- AIC(object, k = log(nobs(object)))

  return(bic)
}

nobs <- function(object, ...){
  val <- nrow(object$data)
  return(val)
}




# @method nobs cm
# nobs.cm <- function(object, ...){
#   val <- nrow(object$data)
#   return(val)
# }
#
# # @method nobs gfcm
# nobs.gfcm <- function(object, ...){
#   val <- nrow(object$data)
#   return(val)
# }
#
