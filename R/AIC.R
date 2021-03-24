#' Akaike information criterion
#'
#' This function computes the AIC for a fitted parametric cure model
#'
#' @param object Object of class cm or gfcm.
#' @param k Number to control if either AIC or BIC is to be computed (default is 2 equal to AIC)
#' @export
#' @method AIC cm
AIC.cm <- function(object, k = 2){
  npar <- length(object$optim$par)
  lnL <- object$ML
  aic <- as.vector(-2*lnL + k * npar)

  return(aic)
}

#' @export
#' @method BIC cm
BIC.cm <- function(object, ...)
{
  bic <- AIC(object, k = log(nobs(object)))

  return(bic)
}

#' @method nobs cm
nobs.cm <- function(object, ...){
  val <- nrow(object$data)
  return(val)
}

