#' Loss of lifetime estimation
#'
#' Function for computing loss of lifetime estimates from an estimated relative survival model
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fmcm}, \code{stpm2}, \code{pstpm2}, and \code{CureModel}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and \code{tau} is selected.
#' @param type Character indicating the type of life expectation estimate.
#' Possible values are \code{ll} (default) which gives the loss of lifetime and \code{mrl},
#' which gives the mean residual lifetime
#' @param tau The upper limit of the integral. Default is 100.
#' @param ci Logical indicating whether confidence intervals should be computed
#' @param ratetable Object of class \code{ratetable} to compute the general population survival from.
#' @param expected Object of class \class{list} containing objects of class \code{survexp}
#' denoting the expected survival of each row in newdata. If not specified, the function computes the expected
#' survival.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @return A object of class \code{lol} containing the loss of lifetime estiamtes
#' of each individual in \code{newdata}.
#' @export
#' @examples
#' D$bhaz <- extract_general(time = "FU", age = "agedays", sex = "sex",
#'                           date = "dx", data = D, ratetable = survexp.dk)
#' fit <- stpm2(Surv(FUyear, status2) ~ 1, data = D, df = 2, bhazard = D$bhaz)
#' res <- calc.LL(fit, time = seq(0, 20, length.out = 100),
#'                rmap = list(age = agedays, sex = sex, year = dx))
#' plot(res)



calc.LL <- function(fit, newdata = NULL, time = NULL, type = "ll",
                    tau = 100, ci = T, expected = NULL, ratetable = survexp.dk,
                    rmap = rmap){

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        #if(class(data) == "list") data <- do.call(cbind, data)
        newdata <- data.frame(arbritary_var = 0)
      }else{
        data <- fit$data
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
  if(any(class(fit) %in% c("stpm2", "pstpm2"))){
    if(class(fit) == "stpm2"){
      response_name <- all.vars(fit@call.formula)[1]
    }else{
      response_name <- all.vars(fit@fullformula)[1]
    }

    fit_tmp <- fit
    rel_surv <- lapply(1:length(expected), function(i){
      function(t, pars){
        res <- rep(NA, length(t))
        fit_tmp@fullcoef <- pars
        wh <- which(t != 0)
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t[wh]))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        res[wh] <- as.numeric(predict(fit_tmp, newdata = newdata_tmp, type = "surv"))
        res[-wh] <- 1
        res
      }
      })
    model.params <- fit@fullcoef
    cov <- fit@vcov
  }else{
    rel_surv <- lapply(1:length(expected), function(i){
      function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                time = t, pars = pars, ci = F)$res[[1]]$Est
      })
    model.params <- c(unlist(fit$coefs), fit$coefs.spline)
    cov <- fit$covariance
  }

  .calcArea <- switch(type,
                      ll = .calcArea.LL,
                      mrl = .calcArea.MRL)

  Ests <- lapply(1:length(expected), function(i){
    #Calculate loss of lifetime
    Est <- .calcArea(rel_surv[[i]], exp_function, time = time,
                     tau = tau, pars = model.params, expected[[i]])
    res <- data.frame(Est)
    names(res) <- type
    if(ci){
      #Calculate variances numerically by the delta method
      J <- pracma::jacobian(.calcArea, x = model.params, rel_surv = rel_surv[[i]],
                    exp_function = exp_function, time = time, tau = tau,
                    expected = expected[[i]])
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% cov %*% x)
      res$lower.ci <- res[, type] - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res[, type] + sqrt(res$Var) * qnorm(0.975)
    }
    res
  })
  all_res <- list(Ests = Ests, time = time, type = type)
  class(all_res) <- "le"
  all_res
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

.calcArea.LL <- function(rel_surv, exp_function, time, tau, pars, expected){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
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

.calcArea.MRL <- function(rel_surv, exp_function, time, tau, pars, expected){
  t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
  df_time <- -diff(t_new)
  mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
  vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
  vals_pop <- rev(vals_pop[t_new %in% time])
  vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
}

exp_function <- function(t, expected){
  s <- summary(expected, t)
  names(s$surv) <- s$time
  a <- s$surv[as.character(t)]
  names(a) <- NULL
  a
}
