
#' Loss of lifetime estimation
#'
#' Function for computing loss of lifetime function based on relative survival models.
#'
#' @param object Fitted model to do predictions from. Possible classes are
#' \code{fmc}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions.
#' If empty, a grid of 100 time points between 0 and \code{tau} is selected.
#' @param type Type of life expectation estimate.
#' Possible values are \code{ll} (default) which gives the loss of lifetime, and \code{mrl}
#' which gives the mean residual lifetime.
#' @param tau The upper limit of the integral. Default is 100.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are computed.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}
#' @param expected Object of class \code{list} containing objects of class \code{survexp},
#' with the expected survival of each row in newdata. If not specified, the function computes the expected
#' survival.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @return An object of class \code{le} containing the life expectancy estimates
#' of each individual in \code{newdata}.
#' @details The function computes
#' \deqn{L(t) = \frac{\int_0^\infty S^*(u)}{S^*(t)} - \frac{\int_0^\infty S(u)}{S(t)}}
#' for a given t, where \eqn{S^*(t)} and \eqn{S(t)} is the survival function in the general
#' population and the patient population, respectively. The integral is computed by Gauss-Legendre quadrature
#' and the point wise variance is estimated using the delta method and numerical differentiation.
#' If \code{type = mrl}, only \eqn{\int_0^\infty S^*(u) / S^*(t)} is computed.
#' @export
#' @example inst/calc.LL.ex.R



calc.LL <- function(object, newdata = NULL, time = NULL, type = c("ll", "mrl"),
                    tau = 100, ci = T, expected = NULL, ratetable = survexp.dk,
                    rmap, pars = NULL, n = 100){

  type <- match.arg(type)
  if(!type %in% c("ll", "mrl"))
    stop("Argument 'type' is wrongly specified, must be either 'll' and 'mrl'")

  #Replace coefficients if new ones are provided
  if(!is.null(pars)){
    if(any(class(object) %in% c("stpm2", "pstpm2"))){
      object@fullcoef <- pars
    } else {
      object$coefs <- pars[1:length(object$coefs)]
      object$coefs.spline <- pars[(length(object$coefs) + 1):length(pars)]
    }
  }

  #Time points at which to evaluate integral
  if(is.null(time)){
    time <- seq(0, tau, length.out = 100)
  }

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.05)

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
        res <- rep(NA, length(t))
        object_tmp@fullcoef <- pars
        wh <- t != 0
        if(any(wh)){
          suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t[wh]))
          names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
          res[wh] <- as.numeric(predict(object_tmp, newdata = newdata_tmp, type = "surv")) 
        }
        res[!wh] <- 1
        res
      }
    })
    model.params <- object@fullcoef
    cov <- object@vcov
  } else {
    if ("gfcm" %in% class(object)) {
      rel_surv <- lapply(1:length(expected), function(i){
        function(t, pars){
          res <- rep(NA, length(t))
          wh <- t != 0
          if(any(wh)){
            res[wh] <- predict(object, newdata = newdata[i,, drop = F],
                               time = t[wh], pars = pars, ci = F)[[1]]$Estimate
          }
          res[!wh] <- 1
          res
        }
      })
    } else {
      rel_surv <- lapply(1:length(expected), function(i){
        function(t, pars){
          res <- rep(NA, length(t))
          wh <- t != 0
          if(any(wh)){
            res[wh] <- predict(object, newdata = newdata[i,, drop = F],
                               time = t[wh], pars = pars, ci = F)$res[[1]]$Est
          }
          res[!wh] <- 1
          res
        }
      })
    }
    model.params <- c(unlist(object$coefs), object$coefs.spline)
    cov <- object$covariance
  }

  gaussxw <- statmod::gauss.quad(n)

  Ests <- lapply(1:length(expected), function(i){
    #Calculate loss of lifetime
    Est <- calcMean(rel_surv[[i]], exp_function, time = time,
                    tau = tau, pars = model.params, expected = expected[[i]],
                    gaussxw = gaussxw)

    if(type == "ll"){
      Est_exp <- calcExpMean(exp_function, time = time,
                             tau = tau, expected = expected[[i]],
                             gaussxw = gaussxw)
      Est <- Est_exp - Est
    }

    res <- data.frame(Est)
    names(res) <- type
    if(ci){
      #Calculate variances numerically by the delta method
      J <- numDeriv::jacobian(calcMean, x = model.params, rel_surv = rel_surv[[i]],
                              exp_function = exp_function, time = time, tau = tau,
                              expected = expected[[i]], gaussxw = gaussxw)
      res$Var <- apply(J, MARGIN = 1, function(x) x %*% cov %*% x)
      res$lower.ci <- res[, type] - sqrt(res$Var) * qnorm(0.975)
      res$upper.ci <- res[, type] + sqrt(res$Var) * qnorm(0.975)
    }
    if(type == "ll2"){
      names(res)[names(res) == type] <- "ll"
    }
    res
  })
  all_res <- list(Ests = Ests, time = time, type = ifelse(type == "ll2", "ll", type), ci = ci)
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

# .calcArea.LL <- function(rel_surv, exp_function, time, tau, pars, expected, gaussxw){
#   t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
#   df_time <- -diff(t_new)
#   exp_eval <- exp_function(t_new, expected)
#   surv_eval <- rel_surv(t_new, pars) * exp_eval
#   surv_diff <- diff(surv_eval)
#   inner <- abs(surv_diff) / 2 + pmin(surv_eval[-length(surv_eval)], surv_eval[-1])
#   vals_pop <- cumsum(c(0, inner * df_time))
#   surv_diff <- diff(exp_eval)
#   inner <- abs(surv_diff) / 2 + pmin(exp_eval[-length(exp_eval)], exp_eval[-1])
#   vals_exp <- cumsum(c(0, inner * df_time))
#   these <- t_new %in% time
#   rev(vals_exp[these]) / rev(exp_eval[these]) - rev(vals_pop[these]) / rev(surv_eval[these])
# }

calcMean <- function(rel_surv, exp_function, time, tau, pars, expected, gaussxw){
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_rel_t <- rel_surv(time, pars)
  eval_gen_t <- exp_function(time, expected)
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * gaussxw$nodes + scale2[i]
    eval_gen <- exp_function(points, expected)
    eval_rel <- rel_surv(points, pars)
    eval[i] <- sum(gaussxw$weights * (eval_gen * eval_rel / eval_rel_t[i]))
  }
  scale * eval / eval_gen_t
}

calcExpMean <- function(exp_function, time, tau, expected, gaussxw){
  scale <- (tau - time) / 2
  scale2 <- (tau + time) / 2
  eval_gen_t <- exp_function(time, expected)
  eval <- rep(NA, length(time))
  for(i in 1:length(time)){
    points <- scale[i] * gaussxw$nodes + scale2[i]
    eval_gen <- exp_function(points, expected)
    eval[i] <- sum(gaussxw$weights * (eval_gen))
  }
  scale * eval / eval_gen_t
}

# .calcArea.MRL <- function(rel_surv, exp_function, time, tau, pars, expected, gaussxw){
#   t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
#   df_time <- -diff(t_new)
#   mid_points <- t_new[-length(t_new)] + diff(t_new) / 2
#   vals_pop <- c(0, cumsum(rel_surv(mid_points, pars) * exp_function(mid_points, expected) * df_time))
#   vals_pop <- rev(vals_pop[t_new %in% time])
#   vals_pop / (rel_surv(time, pars) * exp_function(time, expected))
# }

exp_function <- function(t, expected){
  s <- summary(expected, t)
  names(s$surv) <- s$time
  a <- s$surv[as.character(t)]
  names(a) <- NULL
  a
}
