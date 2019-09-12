#' Long term survival predictions
#'
#' Function for computing survival estimates using a relative survival model and the
#' expected general population survival.
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{gfcm}, \code{stpm2},
#' \code{pstpm2}, and \code{cm}.
#' @param type Prediction type (see details). The default is \code{surv}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and the maximum follow-up time is selected.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}.
#' @param exp.fun Object of class \code{list} containing functions for the expected survival
#' of each row in \code{newdata}. If not specified, the function computes the expected
#' survival based on \code{newdata} using the \code{survival::survexp} function. If \code{newdata} is not provided,
#' the expected survival is based on the data which the model was fitted on.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package if \code{exp.fun = NULL}.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @param scale Numeric. Passed to the \code{survival::survexp} function and defaults to 365.24.
#' That is, the time scale is assumed to be in years.
#' @param smooth.exp Logical. If \code{TRUE}, the general population survival function is smoothed by the function
#' \code{smooth.spline} using the the argument \code{all.knots = TRUE}.
#' @param link Character, indicating the link function for the variance calculations.
#' Possible values are "\code{log}", "\code{cloglog}" for \eqn{log(-log(x))} , "\code{mlog}" for -log(x),
#' and "\code{I}" for the indentity.
#' @return An object of class \code{lts} containing the predictions of each individual in \code{newdata}.
#' @details
#' Possible values for argument \code{type} are:\cr
#' \code{surv}: Survival function computed by \eqn{S(t) = R(t)S^*(t)}\cr
#' \code{hazard}: Hazard function computed by \eqn{h(t) = \lambda(t) + h^*(t)}\cr
#' \code{cumhaz}: The cumulative hazard function computed by \eqn{H(t) = \Lambda(t) + H^*(t)}\cr
#' \code{loghazard}: The log-hazard function computed by \eqn{\log(\lambda(t) + h^*(t))}\cr
#' \code{fail}: The distribution function computed by \eqn{1 - R(t)S^*(t)}
#' @export
#' @example inst/predict.lts.ex.R


lts <- function(fit, type = c("surv", "hazard", "cumhaz", "loghaz", "fail"),
                newdata = NULL, time = NULL, var.type = c("ci", "se", "n"),
                exp.fun = NULL, ratetable = cuRe::survexp.dk, rmap, scale = 365.24,
                smooth.exp = FALSE, link = NULL){

  var.type <- match.arg(var.type)
  type <- match.arg(type)

  if(is.null(time)){
    if(any(class(fit) %in% c("stpm2", "pstpm2"))){
      time <- seq(1e-05, max(fit@data[[fit@timeVar]]), length.out = 100)
    } else {
      time <- seq(1e-05, max(fit$time), length.out = 100)
    }
  }

  is_null_newdata <- is.null(newdata)
  if(is_null_newdata){
    if(any(class(fit) %in% c("stpm2", "pstpm2"))){
      data <- fit@data
      newdata <- data.frame(arbritary_var = 0)
    }else{
      data <- fit$data
    }
  }


  if(is.null(exp.fun)){
    #The time points for the expected survival
    times <- seq(0, max(time) + 1, by = 0.1)

    #Extract expected survival function
    if(is_null_newdata){
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = scale, times = times * scale)))
    }else{
      expected <- vector("list", nrow(newdata))
      for(i in 1:length(expected)){
        expected[[i]] <- do.call("survexp",
                                 list(formula = ~ 1, rmap = substitute(rmap),
                                      data = newdata[i, ], ratetable = ratetable,
                                      scale = scale, times = times * scale))
      }
    }
    if(smooth.exp){
      exp.fun <- lapply(1:length(expected), function(i){
        smooth.obj <- smooth.spline(x = expected[[i]]$time, y = expected[[i]]$surv, all.knots = T)
        function(time) predict(smooth.obj, x = time)$y
      })
    } else {
      exp.fun <- lapply(1:length(expected), function(i){
        function(time){
          s <- summary(expected[[i]], time)
          names(s$surv) <- s$time
          survs <- s$surv[as.character(time)]
          names(survs) <- NULL
          survs
        }
      })
    }
  }

  if(any(class(fit) %in% c("stpm2", "pstpm2"))){
    if(class(fit) == "stpm2"){
      response_name <- all.vars(fit@call.formula)[1]
    }else{
      response_name <- all.vars(fit@fullformula)[1]
    }

    fit_tmp <- fit
    rel_surv <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        #res <- rep(NA, length(t))
        fit_tmp@fullcoef <- pars
        #wh <- which(t != 0)
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        #res[wh] <-
        as.numeric(predict(fit_tmp, newdata = newdata_tmp, type = "surv"))
        #res[-wh] <- 1
        #res
      }
    })

    excess_haz <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        fit_tmp@fullcoef <- pars
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        as.numeric(predict(fit_tmp, newdata = newdata_tmp, type = "hazard"))
      }
    })
    model.params <- fit@fullcoef
    cov <- fit@vcov
  }else{
    if ("cuRe" %in% class(fit)) {
      rel_surv <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, var.type = "n")[[1]]$Estimate
      })

      excess_haz <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, type = "hazard",
                                  var.type = "n")[[1]]$Estimate
      })
    }
    model.params <- c(unlist(fit$coefs), fit$coefs.spline)
    cov <- fit$covariance
  }

  expected_haz <- lapply(1:length(exp.fun), function(i){
    cum_haz_smooth <- function(t) -log(exp.fun[[i]](t))
    function(t, pars) numDeriv::grad(func = cum_haz_smooth, t)
  })

  if(is.null(link)){
    link <- switch(type, surv = "cloglog", hazard = "log", cumhaz = "log", loghazard = "I", fail = "mlog")
  }

  var.link <- switch(link, I = function(x) x, log = function(x) log(x),
                     cloglog = function(x) log(-log(x)), mlog = function(x) -log(x))
  var.link.inv <- switch(link, I = function(x) x, log = function(x) exp(x),
                         cloglog = function(x) exp(-exp(x)), mlog = function(x) exp(-x))

  res <- vector("list", length(exp.fun))
  for(i in 1:length(exp.fun)){
    est <- lts.local(exp.fun = exp.fun[[i]], rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                     expected_haz = expected_haz[[i]], time = time,
                     model.params = model.params, var.link = var.link, type = type)
    D <- data.frame(Estimate = est)
    if(var.type != "n"){
      gr <- numDeriv::jacobian(lts.local, x = model.params, time = time, exp.fun = exp.fun[[i]],
                     rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                     expected_haz = expected_haz[[i]], var.link = var.link, type = type)
      D$SE <- sqrt(apply(gr, 1, function(x) x %*% cov %*% x))
      if(var.type == "ci"){
        lower <- var.link.inv(D$Estimate - D$SE * qnorm(0.975))
        upper <- var.link.inv(D$Estimate + D$SE * qnorm(0.975))
        D$lower <- pmin(lower, upper)
        D$upper <- pmax(lower, upper)
        D <- subset(D, select = -SE)
      }
    }
    D$Estimate <- var.link.inv(D$Estimate)
    res[[i]] <- D
  }
  attributes(res) <- list(time = time, type = type, var.type = var.type)
  class(res) <- "lts"
  res
}



lts.local <- function(exp.fun, rel_surv, excess_haz, expected_haz, time, model.params, var.link, type){
  Est <- if(type == "surv"){
    exp.fun(time) * rel_surv(time, model.params)
  } else if(type == "hazard"){
    expected_haz(time) + excess_haz(time, model.params)
  } else if(type == "cumhaz"){
    -log(exp.fun(time) * rel_surv(time, model.params))
  } else if(type == "loghazard"){
    log(expected_haz(time) + excess_haz(time, model.params))
  } else if(type == "fail"){
    1 - (exp.fun(time) * rel_surv(time, model.params))
  }
  return(var.link(Est))
}
