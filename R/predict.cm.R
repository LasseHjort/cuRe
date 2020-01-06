#' Predict function for flexible mixture cure model
#'
#' Function for doing predictions for class \code{cm}.
#' @usage \method{predict}{cm}(object, newdata = NULL,
#'         type = c("surv", "curerate", "probcure", "survuncured", "hazarduncured",
#'         "cumhazuncured","densityuncured", "failuncured", "oddsuncured",
#'         "loghazarduncured","hazard", "density", "fail", "loghazard",
#'         "odds", "cumhaz"), time = NULL, var.type = c("ci", "se", "n"),
#'         pars = NULL, link = NULL, keep.attributes = F, \dots)
#' @param object Object of class \code{cm} to do predictions from.
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
#' @method predict cm


predict.cm <- function(object, newdata = NULL,
                       type = c("surv", "curerate", "probcure", "survuncured", "hazarduncured",
                                "cumhazuncured", "densityuncured", "failuncured", "oddsuncured",
                                "loghazarduncured", "hazard", "density", "fail",
                                "loghazard", "odds", "cumhaz"),
                       time = NULL, var.type = c("ci", "se", "n"), pars = NULL,
                       link = NULL, keep.attributes = F, ...){
  type <- match.arg(type)
  if(!is.null(pars)){
    groups <- factor(rep(1:length(object$coefs), object$n.param.formula), 1:length(object$coefs))
    object$coefs <- split(pars, f = groups)
  }
  is_null_newdata <- is.null(newdata)

  #Check if covariates are included in the model in cases where newdata is not provided
  if(is_null_newdata){
    classes <- sapply(lapply(object$all.formulas, rhs), class)
    if(any(classes != "numeric")){
      stop("'newdata' must be specified for model including covariates")
    }
    newdata <- data.frame(x = 1)
    #colnames(newdata) <- "(Intercept)"
  }

  #if(is.null(time)) time <- 0
  if (is.null(time)) {
    dtimes <- object$data[[object$timeVar]][object$event]
    time <- seq(min(dtimes), max(dtimes), length.out = 300)[-1]
  }
  if(type == "curerate"){
    time <- 1
  }

  all.formulas <- lapply(object$all.formulas, function(x){
    lhs(x) <- NULL
    x
  }
  )

  X.all <- lapply(all.formulas, get_design, data = newdata)

  if(object$ci){
    var.type <- match.arg(var.type)
  } else {
    var.type <- "n"
  }

  pred <- if (!var.type %in% c("ci", "se")) {
    lapply(1:nrow(newdata), function(i){
      data.frame(Estimate = local.cm(object, type = type, X.all = lapply(X.all, function(X) X[i,, drop = F]),
                                     time = time))
    })
  } else {
    gd <- NULL
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


    lapply(1:nrow(newdata), function(i){
      res <- predictnl.default.cm(object, local.cm, var.link = var.link, type = type,
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


