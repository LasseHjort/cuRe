#' Predict function for flexible mixture cure model
#'
#' Function for doing predictions for class \code{gfcm}.
#'
#' @param object Object of class \code{gfcm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the data which
#' the model was fitted on.
#' @param type Prediction type (see details). The default is \code{surv}.
#' @param time Optional time points at which to compute predictions.
#' This argument is not used if type is \code{curerate}.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param pars Numerical vector containing the parameters values of the model.
#' In general, this argument can be ignored by the user.
#' @param link Character, indicating the link function for the variance calculations.
#' Possible values are "\code{log}", "\code{cloglog}" for \eqn{log(-log(x))} , "\code{mlog}" for -log(x),
#' and "\code{I}" for the indentity.
#' If \code{NULL} (default), the function will determine \code{link} from \code{type}.
#' @param keep.attributes Logical. If \code{TRUE}, \code{newdata} will be added to the attributes of the output.
#' @param indi Logical. If \code{TRUE}, each line in \code{newdata} is treated as an individual observations. This
#' functionality allows predictions for each observation at more than one time point.
#' @param ... Additional arguments. Currently not used.
#' @return A list containing the predictions of each individual in \code{newdata}.
#' @details
#' Possible values for argument \code{type} are:\cr
#' \code{surv}: Survival function\cr
#' \code{curerate}: The cure fraction\cr
#' \code{probcure}: The conditional probability of being cured\cr
#' \code{survuncured}: The survival of the uncured\cr
#' \code{hazarduncured}: The hazard function of the uncured\cr
#' \code{cumhazuncured}: The cumulative hazard of the uncured\cr
#' \code{densityuncured}: The density function of the uncured\cr
#' \code{failuncured}: The distribution function of the uncured, i.e., 1 - \code{survuncured}\cr
#' \code{oddsuncured}: Odds of the uncured, i.e., (1 - \code{survuncured}) / \code{survuncured}\cr
#' \code{loghazarduncured}: The log-hazard of the uncured\cr
#' \code{hazard}: The hazard function\cr
#' \code{density}: The density function\cr
#' \code{fail}: The distribution function\cr
#' \code{loghazard}: The log-hazard function\cr
#' \code{odds}: The odds, i.e., (1 - \code{surv}) / \code{surv}\cr
#' \code{cumhaz}: The cumulative hazard function
#' @export

predict.gfcm <- function (object, newdata = NULL,
                          type = c("surv", "curerate", "probcure", "survuncured", "hazarduncured",
                                   "cumhazuncured", "densityuncured", "failuncured", "oddsuncured",
                                   "loghazarduncured", "hazard", "density", "fail",
                                   "loghazard", "odds", "cumhaz"),
                          indi = TRUE, time = NULL, var.type = c("ci", "se", "n"), pars = NULL,
                          link = NULL, keep.attributes = FALSE, ...)
{
  use.gr = FALSE
  type <- match.arg(type)
  args <- object$args

  if(!is.null(pars)){
    object$coefs <- pars[1:length(object$coefs)]
    object$coefs.spline <- pars[(length(object$coefs) + 1):length(pars)]
  }

  calcX <- !is.null(newdata)
  if (is.null(newdata)) {
    if(indi){
      vars <- c(all.vars(formula), all.vars(object$cr.formula))
      vars <- vars[!vars %in% c(as.character(object$timeExpr), as.character(object$eventExpr))]
      if(length(vars) != 0){
        stop("'newdata' needs to be specified with option 'indi = TRUE' when covariates are present")
      }
      newdata <- data.frame(x = 1)
      colnames(newdata) <- "(Intercept)"
    } else {
      X <- object$args$X
      XD <- object$args$XD
      X.cr <- object$args$X.cr
      y <- object$args$time
      time <- y
      newdata <- as.data.frame(object$data)
    }
  }

  lpfunc <- function(delta, fit, data, var) {
    data[[var]] <- data[[var]] + delta
    lpmatrix.lm(fit, data)
  }

  if (is.null(time)) {
    if(indi){
      dtimes <- object$data[[object$timeVar]][object$event]
      time <- seq(min(dtimes), max(dtimes), length.out = 300)[-1]
    }else{
      time <- eval(object$timeExpr, newdata, parent.frame())
    }
    if(type == "curerate"){
      time <- 1
    }
  }

  if(indi){
    newdata.list <- split(newdata, f = 1:nrow(newdata))
    for(i in 1:length(newdata.list)){
      newdata.list[[i]] <- newdata.list[[i]][rep(1, length(time)),, drop = F]
      newdata.list[[i]][, object$timeVar] <- time
    }

    X <- lapply(1:length(newdata.list), function(i){
      object$transX(lpmatrix.lm(object$lm.obj, newdata.list[[i]]), newdata.list[[i]])
    })

    XD <- lapply(1:length(newdata.list), function(i){
      XD.tmp <- grad(lpfunc, 0, object$lm.obj, newdata.list[[i]], object$timeVar)
      object$transXD(matrix(XD.tmp, nrow = nrow(X[[i]])))
    })

    X.cr <- lapply(1:length(newdata.list), function(i){
      lpmatrix.lm(object$lm.obj.cr, newdata.list[[i]])
      #model.matrix(object$cr.formula, data = newdata.list[[i]])
    })

  } else {
    if(calcX){
      X <- object$transX(lpmatrix.lm(object$lm.obj, newdata),
                         newdata)
      XD <- grad(lpfunc, 0, object$lm.obj, newdata, object$timeVar)
      XD <- object$transXD(matrix(XD, nrow = nrow(X)))
      X.cr <- model.matrix(object$cr.formula, data = newdata)
    }
  }


  var.type <- match.arg(var.type)
  pred <- if (!var.type %in% c("ci", "se")) {
    if(is.list(X)){
      lapply(1:length(X), function(i){
        data.frame(Estimate = local(object, newdata, type, var.link = function(x) x,
                                    X = X[[i]], XD = XD[[i]], X.cr = X.cr[[i]]))
      })
    } else {
      data.frame(Estimate = local(object, newdata, type, var.link = function(x) x,
                                  X = X, XD = XD, X.cr = X.cr))
    }
  } else {
    gd <- NULL
    #beta <- object$coefs.spline
    if (is.null(link)){
      if(!object$excess){
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
                       densityuncured = "I", failuncured = "mlog",
                       oddsuncured = "cloglog", loghazarduncured = "I",
                       surv = "log", hazard = "I", density = "I", fail = "mlog",
                       loghazard = "I", odds = "cloglog", cumhaz = "I")
      }
    }

    var.link <- switch(link, I = function(x) x, log = function(x) log(x),
                       cloglog = function(x) log(-log(x)), mlog = function(x) -log(x))
    var.link.inv <- switch(link, I = function(x) x, log = function(x) exp(x),
                           cloglog = function(x) exp(-exp(x)), mlog = function(x) exp(-x))

    if (use.gr) {
      if (type == "hazard" && link %in% c("I", "log")) {
        betastar <- beta
        gd <- switch(link, I = t(object$link.surv$gradh(X %*%
                                                          betastar, XD %*% betastar, list(X = X, XD = XD))),
                     log = t(object$link.surv$gradh(X %*% betastar, XD %*%
                                                      betastar, list(X = X, XD = XD))/object$link.surv$h(X %*%
                                                                                                           betastar, XD %*% betastar)))
      }
    }

    lapply(1:length(X), function(i){
      res <- predictnl.default(object, local, var.link = var.link, newdata = newdata[i,, drop = F],
                               type = type, gd = if (use.gr) gd else NULL,
                               X = X[[i]], XD = XD[[i]], X.cr = X.cr[[i]])
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



predictnl.default <- function (object, fun, newdata = NULL, gd = NULL, ...)
{
  if (is.null(newdata) && !is.null(object$data))
    newdata <- object$data
  localf <- function(coef, ...) {
    object$coefs <- coef[1:length(object$coefs)]
    object$coefs.spline <- coef[(length(object$coefs) + 1):length(coef)]
    fun(object, ...)
  }
  numDeltaMethod(object, localf, newdata = newdata, gd = gd,
                 ...)
}


numDeltaMethod <- function (object, fun, gd = NULL, ...)
{
  coef <- c(object$coefs, object$coefs.spline)
  est <- fun(coef, ...)
  Sigma <- object$covariance
  if (is.null(gd))
    gd <- grad(fun, coef, ...)
  se.est <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}

local <- function(object, newdata, type = "surv", var.link = function(x) x,
                  X = X, XD = XD, X.cr = X.cr) {
  gamma <- object$coefs
  beta <- object$coefs.spline
  #tt <- object$lm.obj$terms
  link.surv <- object$link.surv
  link.type.cr <- object$link.type.cr
  eta_pi <- as.vector(X.cr %*% gamma)
  eta <- as.vector(X %*% beta)
  etaD <- as.vector(XD %*% beta)

  pi <- get.link(link.type.cr)(eta_pi)
  Su <- link.surv$ilink(eta)
  hazu <- link.surv$h(eta, etaD)
  Hu = link.surv$H(eta)
  S <- object$cure.type$surv(pi, Su)
  dSu <- link.surv$gradS(eta, etaD)
  haz <- object$cure.type$haz(pi, dSu, S)
  if (!object$excess && any(haz < 0))
    warning(sprintf("Predicted hazards less than zero (n=%i).",
                    sum(haz < 0)))
  H <- -log(S)
  Sigma = object$covariance
  est <- switch(type, linkS = eta, linkpi = eta_pi, curerate = pi,
                probcure = pi / S, survuncured = Su,
                hazarduncured = hazu, cumhazuncured = Hu,
                densityuncured = Su * hazu, failuncured = 1 - Su,
                oddsuncured = (1 - Su)/Su, loghazarduncured = log(hazu),
                surv = S, hazard = haz, density = S * haz, fail = 1 - S,
                loghazard = log(haz), odds = (1 - S) / S, cumhaz = H)

  est <- var.link(est)
  return(est)
}

expit <- function(x) {
  ifelse(x==-Inf, 0, ifelse(x==Inf, 1, 1/(1+exp(-x))))
}
logit <- function(p) {
  ifelse(p==0, -Inf, ifelse(p==1, Inf, log(p/(1-p))))
} # numerical safety for large values?

