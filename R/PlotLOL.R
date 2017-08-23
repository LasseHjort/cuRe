
# m_0 <- function(rel_surv, exp_function, time, tau, pars, expected){
#   t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
#   df_time <- -diff(t_new)
#   mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
#   vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
#   vals_pop <- rev(vals_pop[t_new %in% time])
#   vals_exp <- c(0, cumsum(exp_function(mid_points, expected) * df_time))
#   vals_exp <- rev(vals_exp[t_new %in% time])
#   vals_exp / exp_function(time, expected) - vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
# }

m_0 <- function(rel_surv, exp_function, time, tau, pars, expected){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 1000))), decreasing = T)
  df_time <- -diff(t_new)
  exp_eval <- exp_function(t_new, expected)
  surv_eval <- rel_surv(t_new, pars) * exp_eval
  surv_diff <- diff(surv_eval)
  inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
  vals_pop <- cumsum(c(0, inner * df_time))
  surv_diff <- diff(exp_eval)
  inner <- abs(surv_diff) / 2 + pmin(exp_eval[-length(exp_eval)], exp_eval[-1])
  vals_exp <- cumsum(c(0, inner * df_time))
  these <- t_new %in% time
  rev(vals_exp[these]) / rev(exp_eval[these]) - rev(vals_pop[these]) / rev(surv_eval[these])
}

exp_function <- function(t, expected){
  s <- summary(expected, t)
  names(s$surv) <- s$time
  a <- s$surv[as.character(t)]
  names(a) <- NULL
  a
}

#Calculation of the loss of lifetime function
calc.LL <- function(fit, time = NULL, ci = T, tau = 100, newdata = NULL, ratetable = survexp.dk){
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

  LOL <- lapply(1:length(expected), function(i){
    #Calculate loss of lifetime
    LL <- m_0(rel_surv[[i]], exp_function, time = time, tau = tau, pars = c(fit$coefs, fit$coefs.spline),
              expected[[i]])
    res <- data.frame(LL = LL)
    if(ci){
      #Calculate variances numerically by the delta method
      J <- jacobian(m_0, x = c(fit$coefs, fit$coefs.spline), rel_surv = rel_surv[[i]],
                    exp_function = exp_function, time = time, tau = tau, expected = expected[[i]])
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      res$lower.ci <- res$LL - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res$LL + sqrt(res$Var) * qnorm(0.975)
    }
    res
  })
  all_res <- list(LOL = LOL, time = time)
  class(all_res) <- "LOL"
  all_res
}

#Plot function for the loss of lifetime function
plot.LOL <- function(obj, ylim = NULL, xlim = NULL, ci = T, col = 1,
                     ylab = NULL, xlab = "Time", add = F){
  if(is.null(ylab)){
    ylab <- "Loss of lifetime"
  }
  if(length(col) == 1){
    col <- rep(col, length(obj$LOL))
  }
  if(is.null(ylim)){
    if(ci){
      ylim <- range(unlist(lapply(obj$LOL, function(x) x[, c("lower.ci", "upper.ci")])))
    }else{
      ylim <- range(unlist(lapply(obj$LOL, function(x) x[, "LL"])))
    }
  }
  for(i in 1:length(obj$LOL)){
    if(i == 1 & !add){
      plot(LL ~ obj$time, data = obj$LOL[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], xlab = xlab, ylab = ylab)
    }else{
      lines(LL ~ obj$time, data = obj$LOL[[i]], col = col[i])
    }
    if(ci){
      lines(lower.ci ~ obj$time, data = obj$LOL[[i]], lty = 2, col = col[i])
      lines(upper.ci ~ obj$time, data = obj$LOL[[i]], lty = 2, col = col[i])
    }
  }
}


m_0_MRL <- function(rel_surv, exp_function, time, tau, pars, expected){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
  df_time <- -diff(t_new)
  mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
  vals_pop <- rev(vals_pop[t_new %in% time])
  vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
}

#Calculation of the loss of lifetime function
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
    EL <- m_0_MRL(rel_surv[[i]], exp_function, time = time, tau = tau, pars = c(fit$coefs, fit$coefs.spline),
              expected[[i]])
    res <- data.frame(EL = EL)
    if(ci){
      #Calculate variances numerically by the delta method
      J <- jacobian(m_0_MRL, x = c(fit$coefs, fit$coefs.spline), rel_surv = rel_surv[[i]],
                    exp_function = exp_function, time = time, tau = tau, expected = expected[[i]])
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      res$lower.ci <- res$EL - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res$EL + sqrt(res$Var) * qnorm(0.975)
    }
    res
  })
  all_res <- list(MRL = MRL, time = time)
  class(all_res) <- "LOL"
  all_res
}
