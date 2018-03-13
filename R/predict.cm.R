
#' Predict function for cure models
#'
#' This function is used to make predictions of the cure models.
#'
#' @param fit Object of class \code{cm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, \code{probcure} for the
#' conditional probability of being cured, and \code{survuncured} for the disease-specific survival of the uncured.
#' @param ci Logical. If \code{TRUE}, confidence intervals are computed.
#' @param pars Numerical vector containing the parameters values of the model.
#' @return An object of class \code{matrix} including the predictions.
#' @export
#'
predict.cm <- function(fit, newdata = NULL, type = c("surv", "curerate", "probcure", "survuncured", "hazarduncured",
                                                     "cumhazuncured", "densityuncured", "failuncured", "oddsuncured",
                                                     "loghazarduncured", "hazard", "density", "fail",
                                                     "loghazard", "odds", "cumhaz"),
                       time = NULL, var.type = c("ci", "se", "n"), pars = NULL, link = NULL, keep.attributes = F){
  type <- match.arg(type)
  if(!is.null(pars)){
    groups <- factor(rep(1:length(fit$coefs), fit$n.param.formula), 1:length(fit$coefs))
    fit$coefs <- split(pars, f = groups)
  }
  is_null_newdata <- is.null(newdata)

  #Check if covariates are included in the model in cases where newdata is not provided
  if(is_null_newdata){
    classes <- sapply(lapply(fit$all.formulas, rstpm2:::rhs), class)
    if(any(classes != "numeric")){
      stop("'newdata' must be specified for model including covariates")
    }
    newdata <- data.frame(x = 1)
    #colnames(newdata) <- "(Intercept)"
  }

  if(is.null(time)) time <- 0

  all.formulas <- lapply(fit$all.formulas, function(x){
    rstpm2:::lhs(x) <- NULL
    x
  }
  )

  X.all <- lapply(all.formulas, get_design, data = newdata)

  if(fit$ci){
    var.type <- match.arg(var.type)
  } else {
    var.type <- "n"
  }

  pred <- if (!var.type %in% c("ci", "se")) {
    lapply(1:nrow(newdata), function(i){
      data.frame(Estimate = local.cm(fit, type = type, X.all = lapply(X.all, function(X) X[i,, drop = F]),
                                     time = time))
    })
  } else {
    gd <- NULL
    if (is.null(link)){
      if(!fit$excess){
        link <- switch(type, linkS = "I", linkpi = "I", curerate = "cloglog",
                       probcure = "cloglog", survuncured = "cloglog",
                       hazarduncured = "log", cumhazuncured = "log",
                       densityuncured = "log", failuncured = "cloglog",
                       oddsuncured = "cloglog", loghazarduncured = "I",
                       surv = "cloglog", hazard = "log", density = "log", fail = "cloglog",
                       loghazard = "I", odds = "cloglog", cumhaz = "log")
      } else {
        link <- switch(type, linkS = "I", linkpi = "I", curerate = "cloglog",
                       probcure = "cloglog", survuncured = "log",
                       hazarduncured = "I", cumhazuncured = "I",
                       densityuncured = "I", failuncured = "log2",
                       oddsuncured = "cloglog", loghazarduncured = "I",
                       surv = "log", hazard = "I", density = "I", fail = "log2",
                       loghazard = "I", odds = "cloglog", cumhaz = "I")
      }
    }

    var.link <- switch(link, I = function(x) x, log = function(x) log(x),
                       cloglog = function(x) log(-log(x)), log2 = function(x) -log(x))
    var.link.inv <- switch(link, I = function(x) x, log = function(x) exp(x),
                           cloglog = function(x) exp(-exp(x)), log2 = function(x) exp(-x))


    lapply(1:nrow(newdata), function(i){
      res <- predictnl.default.cm(fit, local.cm, var.link = var.link, type = type,
                                  X.all = lapply(X.all, function(X) X[i,, drop = F]), time = time)
      if(var.type == "ci"){
        lower <- var.link.inv(res$Estimate - res$SE * qnorm(0.975))
        upper <- var.link.inv(res$Estimate + res$SE * qnorm(0.975))
        res$lower <- pmin(lower, upper)
        res$upper <- pmax(lower, upper)
        res <- subset(res, select = -SE)
      }

      res$Estimate <- var.link.inv(res$Estimate)
      res
    })
  }
  if (keep.attributes)
    attr(pred, "newdata") <- newdata
  return(pred)
}


predictnl.default.cm <- function (fit, fun, ...)
{
  localf <- function(coef, ...) {
    fit$coefs <- split(coef, rep(1:length(fit$coefs), fit$n.param.formula))
    fun(fit, ...)
  }
  numDeltaMethod.cm(fit, localf, ...)
}


numDeltaMethod.cm <- function (fit, fun, ...)
{
  coef <- unlist(fit$coefs)
  est <- fun(coef, ...)
  Sigma <- fit$covariance
  gd <- rstpm2:::grad(fun, coef, ...)
  se.est <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}


local.cm <- function(fit, type = "surv", var.link = function(x) x,
                     X.all, time) {
  lps <- lapply(1:length(fit$coefs), function(i) as.vector(X.all[[i]] %*% fit$coefs[[i]]))
  pi <- as.vector(get.link(fit$link)(lps[[1]]))
  Su <- fit$surv.fun(x = time, lps = lps)
  fu <- fit$dens.fun(x = time, lps = lps)
  Hu <- -log(Su)
  S <- if(fit$type == "mixture") pi + (1 - pi) * Su else pi ^ (1 - Su)
  f <- if(fit$type == "mixture") (1 - pi) * fu else - log(pi) * fu * S
  H <- -log(S)
  est <- switch(type, linkS = lps, linkpi = lps[[1]], curerate = pi,
                probcure = pi / S, survuncured = Su,
                hazarduncured = Su * fu, cumhazuncured = Hu,
                densityuncured = fu, failuncured = 1 - Su,
                oddsuncured = (1 - Su)/Su, loghazarduncured = log(Su * fu),
                surv = S, hazard = S * f, density = f, fail = 1 - S,
                loghazard = log(S * f), odds = (1 - S) / S, cumhaz = H)

  est <- var.link(est)
  return(est)
}

