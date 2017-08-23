m2 <- function(survival_function, exp_function, tau, time, extra = FALSE){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 500))), decreasing = T)
  df_time <- -diff(t_new)
  vals_pop <- cumsum(c(0, survival_function(t_new[-length(t_new)]) * df_time))
  vals_pop <- rev(vals_pop[t_new %in% time])
  vals_exp <- cumsum(c(0, exp_function(t_new[-length(t_new)]) * df_time))
  vals_exp <- rev(vals_exp[t_new %in% time])
  rmrl <- vals_exp / exp_function(time) - vals_pop / survival_function(time)
  if(extra){
    return(data.frame(LossOfLifetime = rmrl, Time = time))
  }else{
    return(data.frame(LossOfLifetime = rmrl / (tau - time), Time = time))
  }
}


m <- function(survival_function, exp_function, time, tau, extra = FALSE){
  od <- order(time, decreasing = T)
  time <- time[od]
  rmrl_pop <- rep(0, length(time))
  rmrl_ref <- rep(0, length(time))
  for(i in 1:length(time)){
    if(i == 1){
      if(time[i] == tau){
        rmrl_pop[i] <- 0
        rmrl_ref[i] <- 0
      }else{
        rmrl_pop[i] <- integrate(survival_function, lower = time[i], upper = tau, subdivisions = 2000,
                                 rel.tol = .Machine$double.eps ^ 0.25)$value
        rmrl_ref[i] <- integrate(exp_function, lower = time[i], upper = tau, subdivisions = 2000,
                                 rel.tol = .Machine$double.eps ^ 0.25)$value
      }
    }else{
      rmrl_pop[i] <- integrate(survival_function, lower = time[i], upper = time[i - 1])$value
      rmrl_ref[i] <- integrate(exp_function, lower = time[i], upper = time[i - 1])$value
    }
  }
  rmrl_pop <- cumsum(rmrl_pop)[od] / survival_function(time[od])
  rmrl_ref <- cumsum(rmrl_ref)[od] / exp_function(time[od])
  if(extra){
    res <- rmrl_ref - rmrl_pop
  }else{
    res <- (rmrl_ref - rmrl_pop) / (tau - time[od]) * year
  }
  data.frame(LossOfLifetime = res, Time = time[od])
}


calc_stat_cure_flex <- function(formula, data, tau, ClinRelThres){
  var <- as.character(formula[3])
  if(var != "1"){
    if(class(data[, var]) != "factor"){
      data[, var] <- factor(data[, var])
    }
  }
  sfit <- flexsurvspline(formula, data = data, k = 3)
  formula2 <- reformulate(as.character(formula[3]), response = "FU")
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = EOT),
                      data = data, ratetable = survexp.dk, scale = year, method = "conditional")

  if(var != "1"){
    uni <- lapply(levels(data[, var]), function(x){
      cat("Calculating for", var, "=", x, "\n")
      exp_function <- function(t) summary(expected, t)$surv[,paste0(var, "=", x)]
      X <- rep(0, nlevels(data[, var]) - 1)
      wh <- which(levels(data[, var]) == x)
      X[wh - 1] <- 1
      surv_function <- function(t){
        (1 - psurvspline(t,
                         gamma = sfit$coefficients[grepl("gamma", names(sfit$coefficients))],
                         knots = sfit$knots)) ^ exp(X %*% sfit$coefficients[grepl(var, names(sfit$coefficients))])
      }
      obj_func <- function(t) m(surv_function, exp_function, tau = tau, time = t)$LossOfLifetime - ClinRelThres
      uniroot.all(obj_func, interval = c(0, tau))
    })
    names(uni) <- levels(data[, var])
  }else{
    exp_function <- function(t) summary(expected, t)$surv
    surv_function <- function(t){
      1 - psurvspline(t,
                      gamma = sfit$coefficients[grepl("gamma", names(sfit$coefficients))],
                      knots = sfit$knots)
    }
    obj_func <- function(t) m(surv_function, exp_function, tau = tau, time = t)$LossOfLifetime - ClinRelThres
    uni <- uniroot.all(obj_func, interval = c(0, tau))
  }
  uni
}


calc_stat_cure_flex_extra <- function(fit, ClinRelThres){
  times <- seq(0, 100, by = 0.2)
  expected <- survexp( ~ 1, rmap = list(age = age, sex = sex, year = diag_date),
                      data = fit$Data, ratetable = survexp.dk, scale = year, times = times * year)
  exp_function <- function(t){
    s <- summary(expected, t)
    names(s$surv) <- s$time
    a <- s$surv[as.character(t)]
    names(a) <- NULL
    a
  }
  surv_function <- function(t) exp_function(t) * (1 - psurvspline(t, gamma = fit$Lambda, knots = fit$knots))
  f <- function(t) exp_function(t) - 0.05
  tau <- uniroot(f, interval = c(0, 100))$root
  uni <- sapply(ClinRelThres, function(x){
    obj_func <- function(t) m2(surv_function, exp_function, tau = tau, time = t, extra = T)$LossOfLifetime - x
    uniroot.all(obj_func, interval = c(0, tau))
  })
  uni
}


calc_stat_cure_flex_extra_cov <- function(fit, ClinRelThres, z, newdata){
  times <- seq(0, 100, by = 0.1)
  formula2 <- reformulate(fit$var)
  expected <- survexp(~ 1, rmap = list(age = age, sex = sex, year = year),
                      data = newdata, ratetable = survexp.dk, scale = year,
                      times = times * year)
  exp_function <- function(t){
    s <- summary(expected, t)
    names(s$surv) <- s$time
    a <- s$surv[as.character(t)]
    names(a) <- NULL
    a
  }
  surv_function <- function(t) exp_function(t) * flex_spline_surv(t, lambda1 = fit$Lambda1,
                                                                  lambda2 = fit$Lambda2, z = z, knots = fit$knots)
  f <- function(t) exp_function(t) - 0.05
  tau <- uniroot(f, interval = c(0, 100))$root
  uni <- sapply(ClinRelThres, function(x){
    obj_func <- function(t) m2(surv_function, exp_function, tau = tau, time = t, extra = T)$LossOfLifetime - x
    uniroot.all(obj_func, interval = c(0, tau))
  })
  ifelse(length(unlist(uni)) == 0, 0, uni)
}


