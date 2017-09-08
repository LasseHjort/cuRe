#Integration function using the trapez method
int.trapez <- function(func, time, pars, n = 10000){
  eps <- .Machine$double.eps
  t_new <- sort(unique(c(seq(eps, max(time), length.out = n), time)))
  t_new <- t_new[t_new != 0]
  df_time <- diff(t_new)
  surv_eval <- func(t_new, pars)
  surv_diff <- diff(surv_eval)
  inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
  vals_pop <- cumsum(c(0, inner * df_time))
  if(any(time == 0)) t_new[1] <- 0
  vals_pop[t_new %in% time]
}

#Integration function using the rectangular method
int.square <- function(func, time, pars){
  t_new <- unique(sort(c(seq(0, max(time), length.out = 100000), time)))
  df_time <- diff(t_new)
  mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  vals_pop <- c(0, cumsum(func(mid_points, pars) * df_time))
  vals_pop[t_new %in% time]
}

prob_cancer <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link){
  if(all(time == 0)){
    return(rep(0, length(time)))
  }else{
    dens <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
    get.inv.link(int.square(dens, time = time, pars = pars), type = link)
  }
}

prob_other <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link){
  if(all(time == 0)){
    return(rep(0, length(time)))
  }else{
    dens <- function(t, pars) rel_surv(t, pars) * expected_haz(t) * exp_function(t, expected)
    get.inv.link(int.square(dens, time = time, pars = pars), type = link)
  }
}

prob_other_time <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point, link){
    dens1 <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
    int_1 <- int.square(dens1, time = time, pars = pars)
    btd <- int.square(dens1, time = last.point, pars = pars)
    dens2 <- function(t, pars) rel_surv(t, pars) * expected_haz(t) * exp_function(t, expected)
    int_2 <- int.square(dens2, time = time, pars = pars)
    get.inv.link(1 - (btd - int_1) / (1 - int_1 - int_2), type = link)
}


#' Crude event probabilities
#'
#' Function for computing crude event probabilties from relative survival models.
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fmcm}, \code{stpm2}, \code{pstpm2}, and \code{CureModel}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Character indicating which measure is desired. Possible values are \code{cancer} (default),
#' \code{other}, and \code{othertime} (see details).
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{last.point} is selected.
#' @param last.point Constant at which the bound to tie probability is calculated. Default is 100.
#' @param ci Logical indicating whether confidence intervals should be computed.
#' @param ratetable Object of class \code{ratetable} to compute the general population survival from. Default is survexp.dk.
#' @param Link Character indicating the used link function for computing confidence intervals.
#' @return An object of class \code{crude} containing the crude probability estimates
#' of each individual in \code{newdata}.
#' @details The types of crude probabilities is the typical measures in relative survival, namely
#' a crude disease-specific probability (\code{cancer}), the probability of dying from other causes than the
#' disease (\code{other}), and the conditional probability of eventually dying from other causes than the
#' disease given survival until time t (\code{othertime}). The proportion of patients bound to die from the
#' disease can be computed by using \code{cancer} and choosing a sufficiently large time point.
#' @references Eloranta, S., et al. (2014) The application of cure models in the presence of competing risks: a tool
#' for improved risk communication in population-based cancer patient survival. \emph{Epidemiology}, 12:86.
#' @export
calc.Crude <- function(fit, newdata = NULL, type = "cancer", time = NULL, last.point = 100,
                       ci = T, ratetable = survexp.dk, link = "crudeprob"){
  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, last.point, length.out = 100)
  }

  #Extract expected survival function
  times <- seq(0, last.point + 1, by = 0.1)
  if(is.null(newdata)){
    expected <- list(survexp(~ 1, rmap = list(age = age, sex = sex, year = diag_date), data = fit$data,
                             ratetable = ratetable, scale = year,
                             times = times * year))
  }else{
    expected <- lapply(1:nrow(newdata), function(x)survexp(~ 1,
                                                           rmap = list(age = age, sex = sex, year = diag_date),
                                                           data = newdata[x,],
                                                           ratetable = ratetable, scale = year,
                                                           times = times * year))
  }

  rel_surv <- lapply(1:length(expected),
                     function(i) function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                                           time = t, pars = pars, ci = F)$res[[1]]$Est)
  excess_haz <- lapply(1:length(expected),
                       function(i) function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                                             time = t, pars = pars,
                                                             type = "ehaz", ci = F)$res[[1]]$Est)

  expected_haz <- lapply(1:length(expected), function(i){
    D <- data.frame(Cum_haz = c(0, -log(summary(expected[[i]])$surv)), Time = c(-0.1, times))
    sm_fit <- loess(Cum_haz ~ Time, data = D, span = 0.1)
    cum_haz_smooth <- function(t) predict(sm_fit, newdata = data.frame(Time = t))
    function(t) numDeriv::grad(func = cum_haz_smooth, t)
    # f <- function(t) numDeriv::grad(func = cum_haz_smooth, t)
    # curve(f, from = 0, to = 100)
    # plot(Cum_haz ~ Time, data = D, type = "l")
    # curve(cum_haz_smooth, from = 0, to = 100, add = T, col = 2)
  })


  if(type == "cancer"){
    probfun <- prob_cancer
  }else if(type == "other"){
    probfun <- prob_other
  }else if(type == "othertime"){
    probfun <- prob_other_time
  }

  probs <- lapply(1:length(expected), function(i){
    prob <- probfun(time = time, rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                    expected_haz = expected_haz[[i]], expected =  expected[[i]],
                    pars = c(fit$coefs, fit$coefs.spline), last.point = last.point, link = link)
    res <- data.frame(prob = prob)
    if(ci){
      prob_gr <- jacobian(probfun, x = c(fit$coefs, fit$coefs.spline), time = time,
                          rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                          expected_haz = expected_haz[[i]], expected =  expected[[i]],
                          last.point = last.point, link = link)
      res$var <- apply(prob_gr, 1, function(x) x %*% fit$covariance %*% x)
      res$lower.ci <- get.link(res$prob - sqrt(res$var) * qnorm(0.975), type = link)
      res$upper.ci <- get.link(res$prob + sqrt(res$var) * qnorm(0.975), type = link)
    }
    res$prob <- get.link(res$prob, type = link)
    if(type %in% c("cancer", "other")){
      res[time == 0,] <- 0
    }
    res
  })

  probs <- list(prob = probs, time = time, ci = ci, type = type)
  class(probs) <- "crude"
  probs
}


