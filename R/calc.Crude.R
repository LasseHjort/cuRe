#' Crude event probabilities
#'
#' Function for computing crude event probabilties based on relative survival models.
#'
#' @param object Fitted model to do predictions from. Possible classes are
#' \code{gfcm}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Probability to compute. Possible values are \code{disease} (default),
#' \code{other}, and \code{condother} (see details).
#' @param time Time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{tau} is selected.
#' @param tau Upper bound of the integral used to compute the probability of disease-related death (see details).
#' The argument is only used for \code{type = condother}. Default is 100.
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
#' @param reverse Logical. If \code{TRUE}, 1 - probability is provided (default is \code{FALSE}).
#' Only applicable for \code{type = condother}.
#' @param scale Numeric. Passed to the \code{survival::survexp} function and defaults to 365.24.
#' That is, the time scale is assumed to be in years.
#' @param link Link function for computing variance in order to restrict confidence intervals to [0, 1].
#' Default is \code{loglog}.
#' @param n Number of knots used for the Gauss-Legendre quadrature.
#' @param smooth.exp Logical. If \code{TRUE}, the general population survival function is smoothed by the function
#' \code{smooth.spline} using the the argument \code{all.knots = TRUE}.
#' @param pars A vector of parameter values for the model given in \code{object}. Currently not used.
#' @return A list containing the crude probability estimates
#' of each individual in \code{newdata}.
#' @details The function estimates crude probabilities by using the relative survival, expected survival,
#' and the cause-specific hazard function.
#' The crude cumulative incidence of disease-related death (\code{type = "disease"}) is
#' \deqn{P(T \leq t, D = disease) = \int_0^t S^*(u) R(u) \lambda(u)du.}
#' The crude cumulative incidence of death from other causes (\code{type = "other"}) is
#' \deqn{P(T \leq t, D = other) = \int_0^t S^*(u) R(u) h^*(u)du.}
#' The conditional probability of eventually dying from other causes than disease (\code{type = "condother"}) is
#' \deqn{P(D = other| T > t) = \frac{P(D = disease) - P(T \leq t, D = disease)}{P(T > t)}.}
#' The probability of disease-related death, P(D = disease),
#' can be computed by using \code{type = "disease"} and choosing a sufficiently large time point.
#' For P(D = other| T>t), the argument \code{tau} controls this time point (default is 100).
#' @references Eloranta, S., et al. (2014) The application of cure models in the presence of competing risks: a tool
#' for improved risk communication in population-based disease patient survival. \emph{Epidemiology}, 12:86.
#' @export
#' @import date
#' @importFrom statmod gauss.quad
#' @example inst/calc.Crude.ex.R

calc.Crude <- function(object, newdata = NULL, type = c("disease", "other", "condother"),
                       time = NULL, tau = 100, reverse = FALSE, var.type = c("ci", "se", "n"),
                       exp.fun = NULL, ratetable = cuRe::survexp.dk, rmap, scale = ayear,
                       smooth.exp = FALSE, pars = NULL, link = "loglog", n = 100){

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
        #res <- rep(NA, length(t))
        object_tmp@fullcoef <- pars
        #wh <- which(t != 0)
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        #res[wh] <-
        as.numeric(predict(object_tmp, newdata = newdata_tmp, type = "surv"))
        #res[-wh] <- 1
        #res
      }
    })

    excess_haz <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        object_tmp@fullcoef <- pars
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        as.numeric(predict(object_tmp, newdata = newdata_tmp, type = "hazard"))
      }
    })
    model.params <- object@fullcoef
    cov <- object@vcov
  }else{
    if ("cuRe" %in% class(object)) {
      rel_surv <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(object, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, var.type = "n")[[1]]$Estimate
      })

      excess_haz <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(object, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, type = "hazard",
                                  var.type = "n")[[1]]$Estimate
      })
    }
    model.params <- c(unlist(object$coefs), object$coefs.spline)
    cov <- object$covariance
  }

  expected_haz <- lapply(1:length(exp.fun), function(i){
    cum_haz_smooth <- function(t) -log(exp.fun[[i]](t))
    function(t, pars) numDeriv::grad(func = cum_haz_smooth, t)
  })


  probfun <- switch(type,
                    disease = prob_cuminc,
                    other = prob_cuminc,
                    condother = cprob_time)

  cs_haz <- switch(type,
                   disease = excess_haz,
                   other = expected_haz,
                   condother = excess_haz)

  gaussxw <- statmod::gauss.quad(n)

  probs <- lapply(1:length(exp.fun), function(i){
    prob <- probfun(time = time, rel_surv = rel_surv[[i]], cs_haz = cs_haz[[i]],
                    exp.fun = exp.fun[[i]], reverse = reverse,
                    pars = model.params, tau = tau, link = link,
                    gaussxw = gaussxw)
    res <- data.frame(Estimate = prob)
    if(var.type %in% c("ci", "se")){
      prob_gr <- numDeriv::jacobian(probfun, x = model.params, time = time,
                                    rel_surv = rel_surv[[i]], cs_haz = cs_haz[[i]],
                                    exp.fun =  exp.fun[[i]], tau = tau,
                                    link = link, reverse = reverse, gaussxw = gaussxw)
      res$SE <- sqrt(apply(prob_gr, 1, function(x) x %*% cov %*% x))
      if(var.type == "ci"){
        ci1 <- get.link(link)(res$Estimate - qnorm(0.975) * res$SE)
        ci2 <- get.link(link)(res$Estimate + qnorm(0.975) * res$SE)
        res <- subset(res, select = -SE)
        res$lower.ci <- pmin(ci1, ci2)
        res$upper.ci <- pmax(ci1, ci2)
      }
    }
    res$Estimate <- get.link(link)(res$Estimate)
    if(type %in% c("disease", "other")){
      res[time == 0,] <- 0
    }
    res
  })

  attributes(probs) <- list(time = time, type = type, reverse = reverse, ci = var.type == "ci")
  class(probs) <- "crude"
  probs
}


prob_cuminc <- function(time, rel_surv, cs_haz, exp.fun, pars,
                        tau, link, reverse, gaussxw){
  scale <- time / 2
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    if(time[i] == 0){
      eval[i] <- 0
    } else {
      points <- scale[i] * (gaussxw$nodes + 1)
      eval_gen <- exp.fun(points)
      eval_rel <- rel_surv(points, pars)
      eval_haz <- cs_haz(points, pars)
      eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel * eval_haz))
    }
  }
  prob <- scale * eval
  prob[time == 0] <- 0
  get.inv.link(link)(prob)
}


cprob_time <- function(time, rel_surv, cs_haz, exp.fun, pars, tau,
                       link, reverse, n, gaussxw){
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  zs <- gaussxw$nodes
  wt <- gaussxw$weights
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * zs + scale2[i]
    eval_gen <- exp.fun(points)
    eval_rel <- rel_surv(points, pars)
    eval_haz <- cs_haz(points, pars)
    eval[i] <- sum(wt * (eval_gen * eval_rel * eval_haz))
  }
  eval_surv_t <- rep(NA, length(time))
  eval_surv_t[time == 0] <- 1
  if(any(time != 0)){
    eval_surv_t[time != 0] <- rel_surv(time[time != 0], pars) * exp.fun(time[time != 0])
  }
  prob <- scale * eval / eval_surv_t
  if(!reverse) prob <- 1 - prob
  get.inv.link(link)(prob)
}
