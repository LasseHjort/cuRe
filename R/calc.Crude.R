#' Crude event probabilities
#'
#' Function for computing crude event probabilties based on relative survival models.
#'
#' @param object Fitted model to do predictions from. Possible classes are
#' \code{gfcm}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Probability to compute. Possible values are \code{cancer} (default),
#' \code{other}, and \code{condother} (see details).
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{tau} is selected.
#' @param tau Upper bound of the cancer related death integral (see details).
#' The argument is only used for \code{type = condother}. Default is 100.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}.
#' @param exp.fun Object of class \code{list} containing functions for the expected survival
#' of each row in \code{newdata}. If not specified, the function computes the expected
#' survival using the \code{survival::survexp} function and smoothing by \code{smooth.spline}.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package if \code{exp.fun = NULL}.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @param reverse Logical. If \code{TRUE}, 1 - probability is provided (default is \code{FALSE}).
#' Only applicable for \code{type = condother}.
#' @param Link Link function for computing variance in order to bound confidence intervals. Default is \code{loglog}.
#' @param n Number of knots used for the Gauss-Legendre quadrature.
#' @return A list containing the crude probability estimates
#' of each individual in \code{newdata}.
#' @details The function estimates crude probabilities by using the relative survival, expected survival,
#' and the cause-specific hazard function.
#' The crude cumulative incidence of cancer related death (\code{type = "cancer"}) is
#' \deqn{P(T \leq t, D = cancer) = \int_0^t S^*(u) R(u) \lambda(u)du.}
#' The crude cumulative incidence of death from other causes (\code{type = "other"}) is
#' \deqn{P(T \leq t, D = other) = \int_0^t S^*(u) R(u) h^*(u)du.}
#' The conditional probability of eventually dying from other causes than cancer (\code{type = "condother"}) is
#' \deqn{P(D = other| T > t) = \frac{P(D = cancer) - P(T \leq t, D = cancer)}{P(T > t)}.}
#' The proportion of patients bound to die from the disease (P(D = cancer))
#' can be computed by using \code{type = "cancer"} and choosing a sufficiently large time point (e.g., 100 years).
#' @references Eloranta, S., et al. (2014) The application of cure models in the presence of competing risks: a tool
#' for improved risk communication in population-based cancer patient survival. \emph{Epidemiology}, 12:86.
#' @export
#' @example inst/calc.Crude.ex.R
#' @import statmod

calc.Crude <- function(object, newdata = NULL, type = c("cancer", "other", "condother"),
                       time = NULL, tau = 100, reverse = FALSE, scale = ayear,
                       var.type = c("ci", "se", "n"), exp.fun = NULL, ratetable = survexp.dk, rmap,
                       link = "loglog", n = 100){

  type <- match.arg(type)
  var.type <- match.arg(var.type)

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }
  if(is.null(exp.fun)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.1)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(object) %in% c("stpm2", "pstpm2"))){
        data <- object@data
        #if(class(data) == "list") data <- do.call(cbind, data)
        newdata <- data.frame(arbritary_var = 0)
      }else{
        data <- object$data
      }
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
    exp.fun <- lapply(1:length(expected), function(i){
      smooth.obj <- smooth.spline(x = expected[[i]]$time, y = expected[[i]]$surv, all.knots = T)
      function(time) predict(smooth.obj, x = time)$y
    })
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
                    cancer = prob_cuminc,
                    other = prob_cuminc,
                    condother = cprob_time)

  cs_haz <- switch(type,
                   cancer = excess_haz,
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
    if(type %in% c("cancer", "other")){
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
