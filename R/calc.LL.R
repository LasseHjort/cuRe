#' Loss of lifetime estimation
#'
#' Function for computing mean residual lifetime and loss of lifetime estimates based on relative survival models.
#'
#' @param object Fitted model to do predictions from. Possible classes are
#' \code{gfcm}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{tau} is selected.
#' @param type Type of life expectation estimate.
#' Possible values are \code{ll} (default) which gives the loss of lifetime, and \code{mrl}
#' which gives the mean residual lifetime.
#' @param tau The upper limit of the integral (see details). Default is 100.
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
#' @param n Number of knots used for the Gauss-Legendre quadrature.
#' @param scale Numeric. Passed to the \code{survival::survexp} function and defaults to 365.24.
#' That is, the time scale is assumed to be in years.
#' @param smooth.exp Logical. If \code{TRUE}, the general population survival function is smoothed by the function
#' \code{smooth.spline} using the the argument \code{all.knots = TRUE}.
#' @param pars A vector of parameter values for the model given in \code{object}. Currently not used.
#' @return An object of class \code{le} containing the life expectancy estimates
#' of each individual in \code{newdata}.
#' @details The mean residual lifetime function and loss of lifetime function are based on numerical
#' integration of the area under the observed and expected conditional survival functions.
#' If \code{type = "ll"}, the function computes
#' \deqn{\frac{\int_t^\infty S^*(u)}{S^*(t)} - \frac{\int_t^\infty S(u)}{S(t)}.}
#' If \code{type = "mrl"}, the function computes
#' \deqn{\frac{\int_t^\infty S(u)}{S(t)},}
#' for a given t. The function \eqn{S^*(t)} is the general population survival function and \eqn{S(t)}
#' is the observed survival function. Integration to infinity is not required in studies of human mortality,
#' so an upper limit, \code{tau}, is chosen instead. As most humans die before they 100 years, this is
#' the default setting of the function. The integral is computed by Gauss-Legendre quadrature
#' and the point wise variance is estimated using the delta method and numerical differentiation.
#' @export
#' @example inst/calc.LL.ex.R



calc.LL <- function(object, newdata = NULL, type = c("ll", "mrl"), time = NULL,
                    tau = 100, var.type = c("ci", "se", "n"), exp.fun = NULL, ratetable = survexp.dk,
                    rmap, smooth.exp = FALSE, scale = ayear, pars = NULL, n = 100){
  type <- match.arg(type)
  var.type <- match.arg(var.type)

  #Replace coefficients if new ones are provided
  # if(!is.null(pars)){
  #   if(any(class(object) %in% c("stpm2", "pstpm2"))){
  #     object@fullcoef <- pars
  #   } else {
  #     object$coefs <- pars[1:length(object$coefs)]
  #     object$coefs.spline <- pars[(length(object$coefs) + 1):length(pars)]
  #   }
  # }

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }

  is_null_newdata <- is.null(newdata)
  if(is_null_newdata){
    if(any(class(object) %in% c("stpm2", "pstpm2"))){
      data <- object@data
      newdata <- data.frame(arbritary_var = 0)
    }else{
      data <- object$data
    }
  }

  if(is.null(exp.fun)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.1)

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


  #Extract relative survival function
  if(any(class(object) %in% c("stpm2", "pstpm2"))){
    if(class(object) == "stpm2"){
      response_name <- all.vars(object@call.formula)[1]
    }else{
      response_name <- all.vars(object@fullformula)[1]
    }

    object_tmp <- object
    rel_surv <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        res <- rep(NA, length(t))
        object_tmp@fullcoef <- pars
        wh <- t != 0
        if(any(wh)){
          suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t[wh]))
          names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
          res[wh] <- as.numeric(predict(object_tmp, newdata = newdata_tmp, type = "surv"))
        }
        res[!wh] <- 1
        res
      }
    })
    model.params <- object@fullcoef
    cov <- object@vcov
  } else {
    if ("cuRe" %in% class(object)) {
      rel_surv <- lapply(1:length(exp.fun), function(i){
        function(t, pars){
          res <- rep(NA, length(t))
          wh <- t != 0
          if(any(wh)){
            res[wh] <- predict(object, newdata = newdata[i,, drop = F],
                               time = t[wh], pars = pars, var.type = "n")[[1]]$Estimate
          }
          res[!wh] <- 1
          res
        }
      })
    }
    model.params <- c(unlist(object$coefs), object$coefs.spline)
    cov <- object$covariance
  }

  gaussxw <- statmod::gauss.quad(n)

  Ests <- lapply(1:length(exp.fun), function(i){
    #Calculate loss of lifetime
    Est <- calcMean(rel_surv[[i]], exp.fun = exp.fun[[i]], time = time,
                    tau = tau, pars = model.params,
                    gaussxw = gaussxw)

    if(type == "ll"){
      Est_exp <- calcExpMean(exp.fun = exp.fun[[i]], time = time,
                             tau = tau, gaussxw = gaussxw)
      Est <- Est_exp - Est
    }

    res <- data.frame(Estimate = Est)
    if(var.type %in% c("ci", "se")){
      #Calculate variances numerically by the delta method
      J <- numDeriv::jacobian(calcMean, x = model.params, rel_surv = rel_surv[[i]],
                              exp.fun = exp.fun[[i]], time = time, tau = tau,
                              gaussxw = gaussxw)
      res$SE <- sqrt(apply(J, MARGIN = 1, function(x) x %*% cov %*% x))
      if(var.type == "ci"){
        res$lower.ci <- res$Estimate - res$SE * qnorm(0.975)
        res$upper.ci <- res$Estimate + res$SE * qnorm(0.975)
        res <- subset(res, select = -SE)
      }
    }
    res
  })
  attributes(Ests) <- list(time = time, type = type, var.type = var.type)
  class(Ests) <- "le"
  Ests
}


calcMean <- function(rel_surv, exp.fun, time, tau, pars, gaussxw){
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_rel_t <- rel_surv(time, pars)
  eval_gen_t <- exp.fun(time)
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * gaussxw$nodes + scale2[i]
    eval_gen <- exp.fun(points)
    eval_rel <- rel_surv(points, pars)
    eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel / eval_rel_t[i]))
  }
  scale * eval / eval_gen_t
}

calcExpMean <- function(exp.fun, time, tau, gaussxw){
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_gen_t <- exp.fun(time)
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * gaussxw$nodes + scale2[i]
    eval_gen <- exp.fun(points)
    eval[i] <- sum(gaussxw$weights * (eval_gen))
  }
  scale * eval / eval_gen_t
}

