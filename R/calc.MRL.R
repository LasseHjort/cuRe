#' Mean residual lifetime
#'
#' Function for computing mean (restricted) residual lifetime from an estimated relative survival model
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fmcm}, \code{stpm2}, \code{pstpm2}, \code{CureModel}, and \code{flexsurvreg}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{tau} is selected.
#' @param tau The upper limit of the integral. Default is 100.
#' @param ci Logical indicating whether confidence intervals should be computed
#' @param ratetable Object of class \code{ratetable} to compute the general population survival from.
#' @return A object of class \code{lol} containing the loss of lifetime estiamtes
#' of each individual in \code{newdata}.
#' @export
#'
calc.MRL <- function(fit, time = NULL, ci = T, tau = 100, newdata = NULL, ratetable = survexp.dk){
  #The time points for the expected survival
  times <- seq(0, tau + 1, by = 0.1)

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }

  #Extract expected survival function
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

  #Extract relative survival function
  rel_surv <- lapply(1:length(expected), function(i) function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                                                               time = t, pars = pars, ci = F)$res[[1]]$Est)

  MRL <- lapply(1:length(expected), function(i){
    #Calculate loss of lifetime
    EL <- .calcArea_MRL(rel_surv[[i]], exp_function, time = time, tau = tau, pars = c(fit$coefs, fit$coefs.spline),
                        expected[[i]])
    res <- data.frame(EL = EL)
    if(ci){
      #Calculate variances numerically by the delta method
      J <- jacobian(.calcArea_MRL, x = c(fit$coefs, fit$coefs.spline), rel_surv = rel_surv[[i]],
                    exp_function = exp_function, time = time, tau = tau, expected = expected[[i]])
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      res$lower.ci <- res$EL - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res$EL + sqrt(res$Var) * qnorm(0.975)
    }
    res
  })
  all_res <- list(MRL = MRL, time = time)
  class(all_res) <- "lol"
  all_res
}


.calcArea_MRL <- function(rel_surv, exp_function, time, tau, pars, expected){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
  df_time <- -diff(t_new)
  mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
  vals_pop <- rev(vals_pop[t_new %in% time])
  vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
}
