#' Compute the time to statistical cure using the loss of lifetime function
#' @param fit Fitted model to do predictions from. Possible classes are
#' \code{fmc}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param q Threshold to estimate statistical cure according to.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param max.time Upper boundary of the interval [0, \code{max.time}] in which to search for solution.
#' @param tau The upper limit of the integral. Default is 100.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are computed.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}
#' @param expected Object of class \code{list} containing objects of class \code{survexp},
#' with the expected survival of each row in newdata. If not specified, the function computes the expected
#' survival.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @param type Type of life expectancy measure. Possible values are "ll" for the loss of lifetime
#' and "mrl" for the mean residual lifetime.
#' @return The estimated cure points.
#' @export


calc.LL.quantile <- function(fit, q = 2, newdata = NULL, max.time = 20, ci = TRUE,
                             expected = NULL, rmap = NULL, ratetable = survexp.dk,
                             tau = 100, type = "ll"){

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
      }else{
        data <- fit$data
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- lapply(1:nrow(newdata), function(x){
        do.call("survexp",
                list(formula = ~ 1, rmap = substitute(rmap),
                     data = newdata[x, ], ratetable = ratetable,
                     scale = ayear, times = times * ayear))
      })
    }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    f <- function(time, q) calc.LL(fit, time = time, ci = F, expected = expected[i],
                                   newdata = newdata[i,,drop = F], tau = tau)$Ests[[1]][, type] - q
    uni <- rootSolve::uniroot.all(f, lower = 0, upper = max.time, q = q)
    if(ci){
      gr <- grad(f, x = uni, q = 0)
      VAR <- calc.LL(fit, time = uni, expected = expected[i],
                     newdata = newdata[i,, drop = F], tau = tau)$Ests[[1]]$Var
      VAR2 <- gr ^ (-2) * VAR
      upper <- uni + sqrt(VAR2) * qnorm(0.975)
      lower <- uni - sqrt(VAR2) * qnorm(0.975)
      data.frame(Est = uni, var = VAR2, lower.ci = lower, upper.ci = upper)
    } else {
      data.frame(Est = uni)
    }
  })
  do.call(rbind, ests)
}




quantile.calc.LL2 <- function(fit, q = 2, newdata = NULL, max.time = 20, ci = TRUE,
                             expected = NULL, rmap, ratetable = survexp.dk,
                             last.point = 100, type = "ll"){

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, last.point + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
        coefs <- fit@coef
        covar <- fit@vcov
      }else{
        data <- fit$data
        coefs <- c(fit$coefs, fit$coefs.spline)
        covar <- fit$covariance
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- lapply(1:nrow(newdata), function(x){
        do.call("survexp",
                list(formula = ~ 1, rmap = substitute(rmap),
                     data = newdata[x, ], ratetable = ratetable,
                     scale = ayear, times = times * ayear))
      })
    }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){

    g <- function(pars, q){
      f <- function(time) calc.LL(fit = fit, time = time, ci = F, expected = expected[i],
                                  newdata = newdata[i,, drop = F], pars = pars)$Ests[[1]][, type] - q
      rootSolve::uniroot.all(f, lower = 0, upper = max.time)
    }

    uni <- g(pars = NULL, q = q)
    if(ci){
      gr <- grad(g, x = coefs, q = q)
      VAR <- gr %*% covar %*% gr
      upper <- uni + sqrt(VAR) * qnorm(0.975)
      lower <- uni - sqrt(VAR) * qnorm(0.975)
      data.frame(Est = uni, var = VAR, lower.ci = lower, upper.ci = upper)
    } else {
      data.frame(Est = uni)
    }
  })
  do.call(rbind, ests)
}
