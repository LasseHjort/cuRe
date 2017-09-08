#' Loss of lifetime estimation
#'
#' Function for computing loss of lifetime estimates from an estimated relative survival model
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fmcm}, \code{stpm2}, \code{pstpm2}, and \code{CureModel}.
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

calc.LL <- function(fit, newdata = NULL, time = NULL, tau = 100, ci = T, ratetable = survexp.dk){
  #The time points for the expected survival
  times <- seq(0, tau + 1, by = 0.1)

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }

  #Extract expected survival function
  if(is.null(newdata)){
    if(class(fit) == "stpm2"){
      data <- fit@data
      newdata <- data.frame(arbritary_var = 0)
    }else{
      data <- fit$data
    }
    expected <- list(survexp(~ 1, rmap = list(age = age, sex = sex, year = diag_date), data = data,
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
  if(class(fit) == "stpm2"){
    fit_tmp <- fit
    rel_surv <- lapply(1:length(expected), function(i){
      function(t, pars){
        res <- rep(NA, length(t))
        fit_tmp@fullcoef <- pars
        response_name <- as.character(fit@call.formula[[2]])[2]
        wh <- which(t != 0)
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t[wh]))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        res[wh] <- as.numeric(predict(fit_tmp, newdata = newdata_tmp))
        res[-wh] <- 1
        res
      }
      })
    model.params <- fit@fullcoef
    cov <- fit@vcov
  }else{
    rel_surv <- lapply(1:length(expected), function(i) function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                                                                 time = t, pars = pars, ci = F)$res[[1]]$Est)
    model.params <- c(fit$coefs, fit$coefs.spline)
    cov <- fit$covariance
  }

  LOL <- lapply(1:length(expected), function(i){
    #Calculate loss of lifetime
    LL <- .calcArea(rel_surv[[i]], exp_function, time = time, tau = tau, pars = model.params,
              expected[[i]])
    res <- data.frame(LL = LL)
    if(ci){
      #Calculate variances numerically by the delta method
      J <- jacobian(.calcArea, x = model.params, rel_surv = rel_surv[[i]],
                    exp_function = exp_function, time = time, tau = tau, expected = expected[[i]])
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% cov %*% x)
      res$lower.ci <- res$LL - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res$LL + sqrt(res$Var) * qnorm(0.975)
    }
    res
  })
  all_res <- list(LOL = LOL, time = time)
  class(all_res) <- "lol"
  all_res
}


plot.lol <- function(obj, ylim = NULL, xlim = NULL, ci = T, col = 1,
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



# .calcArea <- function(rel_surv, exp_function, time, tau, pars, expected){
#   t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
#   df_time <- -diff(t_new)
#   mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
#   vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
#   vals_pop <- rev(vals_pop[t_new %in% time])
#   vals_exp <- c(0, cumsum(exp_function(mid_points, expected) * df_time))
#   vals_exp <- rev(vals_exp[t_new %in% time])
#   vals_exp / exp_function(time, expected) - vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
# }

.calcArea <- function(rel_surv, exp_function, time, tau, pars, expected){
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
