#' Crude event probabilities
#'
#' Function for computing crude event probabilties based on relative survival models.
#'
#' @param object Fitted model to do predictions from. Possible classes are
#' \code{fmc}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Probability to compute. Possible values are \code{cancer} (default),
#' \code{other}, and \code{othertime} (see details).
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{last.point} is selected.
#' @param last.point Upper bound of the cancer related death integral (see details).
#' The argument is only used for \code{type = othertime}. Default is 100.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are computed.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}.
#' @param expected Object of class \code{list} containing objects of class \code{survexp},
#' with the expected survival of each row in \code{newdata}. If not specified, the function computes the expected
#' survival.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package if \code{expected = NULL}.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @param reverse Logical. If \code{TRUE}, 1 - probability is provided (default is \code{FALSE}).
#' Only applicable for \code{type = othertime}.
#' @param Link Link function for computing variance in order to bound confidence intervals. Default is \code{loglog}.
#' @return An object of class \code{crude} containing the crude probability estimates
#' of each individual in \code{newdata}.
#' @details The function estimates crude probabilities by using the relative survival, expected survival,
#' and the cause-specific hazard functions.
#' The crude cumulative incidence of cancer related death (\code{cancer}) is
#' \deqn{P(T \leq t, D = cancer) = \int_0^t S^*(u) R(u) \lambda(u)du.}
#' The crude cumulative incidence of death from other causes (\code{other})is
#' \deqn{P(T \leq t, D = other) = \int_0^t S^*(u) R(u) h^*(u)du.}
#' The conditional probability of eventually dying from other causes than cancer (\code{othertime})
#' \deqn{P(D = other| T > t) = \frac{P(D = cancer) - P(T \leq t, D = cancer)}{P(T > t)}.}
#' The proportion of patients bound to die from the disease (P(D = cancer))
#' can be computed by using \code{cancer} and choosing a sufficiently large time point (default is 100).
#' @references Eloranta, S., et al. (2014) The application of cure models in the presence of competing risks: a tool
#' for improved risk communication in population-based cancer patient survival. \emph{Epidemiology}, 12:86.
#' @export
#' @example inst/calc.Crude.ex.R
#' @import statmod

calc.Crude <- function(object, newdata = NULL, type = c("cancer", "other", "othertime"),
                       time = NULL, last.point = 100, reverse = FALSE,
                       ci = T, expected = NULL, ratetable = survexp.dk, rmap, link = "loglog", n = 100){

  type <- match.arg(type)

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, last.point, length.out = 100)
  }
  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, last.point + 1, by = 0.05)

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
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- vector("list", nrow(newdata))
      for(i in 1:length(expected)){
        expected[[i]] <- do.call("survexp",
                                 list(formula = ~ 1, rmap = substitute(rmap),
                                      data = newdata[i, ], ratetable = ratetable,
                                      scale = ayear, times = times * ayear))
      }
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
    rel_surv <- lapply(1:length(expected), function(i){
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

    excess_haz <- lapply(1:length(expected), function(i){
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
      rel_surv <- lapply(1:length(expected), function(i){
        function(t, pars) predict(object, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, var.type = "n")[[1]]$Estimate
      })

      excess_haz <- lapply(1:length(expected), function(i){
        function(t, pars) predict(object, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, type = "hazard",
                                  var.type = "n")[[1]]$Estimate
      })
    }
    model.params <- c(unlist(object$coefs), object$coefs.spline)
    cov <- object$covariance
  }

  expected_haz <- lapply(1:length(expected), function(i){
    D <- data.frame(Cum_haz = c(0, -log(summary(expected[[i]])$surv)), Time = c(-0.1, expected[[i]]$time))
    sm_fit <- loess(Cum_haz ~ Time, data = D, span = 0.1)
    cum_haz_smooth <- function(t) predict(sm_fit, newdata = data.frame(Time = t))
    function(t, pars) numDeriv::grad(func = cum_haz_smooth, t)
  })


  probfun <- switch(type,
                    cancer = prob_cuminc,
                    other = prob_cuminc,
                    othertime = cprob_time)

  cs_haz <- switch(type,
                   cancer = excess_haz,
                   other = expected_haz,
                   othertime = excess_haz)

  gaussxw <- statmod::gauss.quad(n)

  probs <- lapply(1:length(expected), function(i){
    prob <- probfun(time = time, rel_surv = rel_surv[[i]], cs_haz = cs_haz[[i]],
                    expected =  expected[[i]], reverse = reverse,
                    pars = model.params, last.point = last.point, link = link,
                    gaussxw = gaussxw)
    res <- data.frame(prob = prob)
    if(ci){
      prob_gr <- numDeriv::jacobian(probfun, x = model.params, time = time,
                                    rel_surv = rel_surv[[i]], cs_haz = cs_haz[[i]],
                                    expected =  expected[[i]], last.point = last.point,
                                    link = link, reverse = reverse, gaussxw = gaussxw)
      res$var <- apply(prob_gr, 1, function(x) x %*% cov %*% x)
      ci1 <- get.link(link)(res$prob - qnorm(0.975) * sqrt(res$var))
      ci2 <- get.link(link)(res$prob + qnorm(0.975) * sqrt(res$var))
      res$lower.ci <- pmin(ci1, ci2)
      res$upper.ci <- pmax(ci1, ci2)
    }
    res$prob <- get.link(link)(res$prob)
    if(type %in% c("cancer", "other")){
      res[time == 0,] <- 1
    }
    res
  })

  probs <- list(prob = probs, time = time, ci = ci, type = type, reverse = reverse)
  class(probs) <- "crude"
  probs
}


# #Integration function using the trapez method
# int.trapez <- function(func, time, pars, n = 10000){
#   eps <- .Machine$double.eps
#   t_new <- sort(unique(c(seq(eps, max(time), length.out = n), time)))
#   t_new <- t_new[t_new != 0]
#   df_time <- diff(t_new)
#   surv_eval <- func(t_new, pars)
#   surv_diff <- diff(surv_eval)
#   inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
#   vals_pop <- cumsum(c(0, inner * df_time))
#   if(any(time == 0)) t_new[1] <- 0
#   vals_pop[t_new %in% time]
# }

# #Integration function using the rectangular method
# int.square <- function(func, time, pars, n = 100000){
#   t_new <- unique(sort(c(seq(0, max(time), length.out = n), time)))
#   df_time <- diff(t_new)
#   mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
#   vals_pop <- c(0, cumsum(func(mid_points, pars) * df_time))
#   vals_pop[t_new %in% time]
# }
#
# prob_cancer <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link, reverse, n, gaussxw){
#   if(all(time == 0)){
#     return(rep(0, length(time)))
#   }else{
#     dens <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
#     prob <- int.square(dens, time = time, pars = pars, n = n)
#     if(reverse) prob <- 1 - prob
#     get.inv.link(link)(prob)
#   }
# }
#
# prob_other <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link, reverse, n, gaussxw){
#   if(all(time == 0)){
#     return(rep(0, length(time)))
#   }else{
#     dens <- function(t, pars) rel_surv(t, pars) * expected_haz(t) * exp_function(t, expected)
#     prob <- int.square(dens, time = time, pars = pars, n = n)
#     if(reverse) prob <- 1 - prob
#     get.inv.link(link)(prob)
#   }
# }
#
# prob_other_time <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link, reverse, n, gaussxw){
#   dens1 <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
#   int_1 <- int.square(dens1, time = time, pars = pars, n = n)
#   btd <- int.square(dens1, time = last.point, pars = pars, n = n)
#   wh <- time == 0
#   prob_t <- rep(NA, length(time))
#   prob_t[!wh] <- rel_surv(time[!wh], pars) * exp_function(time[!wh], expected)
#   prob_t[wh] <- 1
#   prob <- 1 - (btd - int_1) / prob_t
#   if(reverse) prob <- 1 - prob
#   get.inv.link(link)(prob)
# }

prob_cuminc <- function(time, rel_surv, cs_haz, expected, pars,
                        last.point, link, reverse, gaussxw){
  scale <- time / 2
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    if(time[i] == 0){
      eval[i] <- 0
    } else {
      points <- scale[i] * (gaussxw$nodes + 1)
      eval_gen <- exp_function(points, expected)
      eval_rel <- rel_surv(points, pars)
      eval_haz <- cs_haz(points, pars)
      eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel * eval_haz))
    }
  }
  prob <- scale * eval
  prob[time == 0] <- 0
  get.inv.link(link)(prob)
}


cprob_time <- function(time, rel_surv, cs_haz, expected, pars, last.point,
                       link, reverse, n, gaussxw){
  scale <- (last.point - time) / 2
  scale2 <- (last.point + time) / 2
  zs <- gaussxw$nodes
  wt <- gaussxw$weights
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * zs + scale2[i]
    eval_gen <- exp_function(points, expected)
    eval_rel <- rel_surv(points, pars)
    eval_haz <- cs_haz(points, pars)
    eval[i] <- sum(wt * (eval_gen * eval_rel * eval_haz))
  }
  eval_surv_t <- rel_surv(time, pars) * exp_function(time, expected)
  prob <- scale * eval / eval_surv_t
  if(!reverse) prob <- 1 - prob
  get.inv.link(link)(prob)
}
