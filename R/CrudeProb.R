#Integration function using the trapez method
int.trapez <- function(func, time, pars, n = 1000){
  t_new <- unique(sort(c(seq(0, max(time), length.out = n), time)))
  df_time <- diff(t_new)
  surv_eval <- func(t_new, pars)
  surv_diff <- diff(surv_eval)
  inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
  vals_pop <- cumsum(c(0, inner * df_time))
  vals_pop[t_new %in% time]
}

#Integration function using the rectangular method
int.square <- function(func, time, pars){
  t_new <- unique(sort(c(seq(0, max(time), length.out = 5000), time)))
  df_time <- diff(t_new)
  mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  vals_pop <- c(0, cumsum(func(mid_points, pars) * df_time))
  vals_pop[t_new %in% time]
}

prob_cancer <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point){
  if(all(time == 0)){
    return(rep(0, length(time)))
  }else{
    dens <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
    get.inv.link(int.square(dens, time = time, pars = pars), type = "crudeprob")
  }
}

prob_other <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point){
  if(all(time == 0)){
    return(rep(0, length(time)))
  }else{
    dens <- function(t, pars) rel_surv(t, pars) * expected_haz(t) * exp_function(t, expected)
    get.inv.link(int.square(dens, time = time, pars = pars), type = "crudeprob")
  }
}

prob_other_time <- function(time, rel_surv, excess_haz, expected_haz, expected, pars, last.point){
    dens1 <- function(t, pars) rel_surv(t, pars) * excess_haz(t, pars) * exp_function(t, expected)
    int_1 <- int.square(dens1, time = time, pars = pars)
    btd <- int.square(dens1, time = last.point, pars = pars)
    dens2 <- function(t, pars) rel_surv(t, pars) * expected_haz(t) * exp_function(t, expected)
    int_2 <- int.square(dens2, time = time, pars = pars)
    get.inv.link(1 - (btd - int_1) / (1 - int_1 - int_2), type = "crudeprob")
}



calc.Crude <- function(fit, time, newdata = NULL, last.point = 100, type = "cancer",
                       ci = T, ratetable = survexp.dk){
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
  }else if(type == "probother"){
    probfun <- prob_other_time
  }

  probs <- lapply(1:length(expected), function(i){
    prob <- probfun(time = time, rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                    expected_haz = expected_haz[[i]], expected =  expected[[i]],
                    pars = c(fit$coefs, fit$coefs.spline), last.point = last.point)
    res <- data.frame(prob = prob)
    if(ci){
      prob_gr <- jacobian(prob_cancer, x = c(fit$coefs, fit$coefs.spline), time = time,
                          rel_surv = rel_surv[[i]], excess_haz = excess_haz[[i]],
                          expected_haz = expected_haz[[i]], expected =  expected[[i]],
                          last.point = last.point)
      res$var <- apply(prob_gr, 1, function(x) x %*% fit$covariance %*% x)
      res$lower.ci <- get.link(res$prob - sqrt(res$var) * qnorm(0.975), type = "crudeprob")
      res$upper.ci <- get.link(res$prob + sqrt(res$var) * qnorm(0.975), type = "crudeprob")
    }
    res$prob <- get.link(res$prob, type = "crudeprob")
    if(type %in% c("cancer", "other")){
      res[time == 0,] <- 0
    }
    res
  })

  probs <- list(prob = probs, time = time, ci = ci, type = type)
  class(probs) <- "CrudeProb"
  probs
}

plot.CrudeProb <- function(obj, ylim = c(0, 1), xlim = NULL, ci = T, col = 1, ylab = NULL, xlab = "Time"){
  if(is.null(ylab)){
    ylab <- switch(obj$type, cancer = "Cumulative incidence of cancer related death",
                   other = "Cumulative incidence of death from other causes than cancer",
                   probother = "Probability of eventually dying from other causes than cancer")
  }
  if(is.null(xlim)) xlim <- range(obj$time)

  if(length(col) == 1){
    col <- rep(col, length(obj$prob))
  }
  if(ci & !obj$ci){
    ci <- FALSE
  }

  for(i in 1:length(obj$prob)){
    if(i == 1){
      plot(prob ~ obj$time, data = obj$prob[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], ylab = ylab, xlab = xlab)
    }else{
      lines(prob ~ obj$time, data = obj$prob[[i]], col = col[i])
    }
    if(ci){
      lines(lower.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i])
      lines(upper.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i])
    }
  }
}
