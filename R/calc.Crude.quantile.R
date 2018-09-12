#' Compute the time to statistical cure using the conditional probability of disease-related death
#'
#' The following function estimates the time to statistical cure using the conditional
#' probability of disease-related death.
#' @param fit Fitted model to do predictions from. Possible classes are
#' \code{gfmc}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param q Threshold to estimate statistical cure according to.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param max.time Upper boundary of the interval [0, \code{max.time}] in which to search for solution (see details).
#' Default is 20.
#' @param tau Upper bound of integral (see ?calc.Crude). Default is 100.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}
#' @param exp.fun Object of class \code{list} containing functions for the expected survival
#' of each row in \code{newdata}. If not specified, the function computes the expected
#' survival based on \code{newdata} using the \code{survival::survexp} function. If \code{newdata} is not provided,
#' the expected survival is based on the data which the model was fitted on.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package if \code{exp.fun = NULL}.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @param reverse Logical passed on to \code{calc.Crude}. If \code{TRUE} (default), 1 - probability is provided.
#' Only applicable for \code{type = condother}.
#' @param scale Numeric. Passed to the \code{survival::survexp} function and defaults to 365.24.
#' That is, the time scale is assumed to be in years.
#' @details The cure point is calculated as the time point at which the conditional probability of disease-related
#' death reaches the threshold, \code{q}. If \code{q} is not reached within \code{max.time}, no solution is reported.
#' @return The estimated cure point.
#' @example inst/calc.Crude.quantile.ex.R
#' @export

calc.Crude.quantile <- function(fit, q = 0.05, newdata = NULL, max.time = 20, exp.fun = NULL,
                                var.type = c("ci", "se", "n"), rmap, ratetable = survexp.dk,
                                tau = 100, reverse = TRUE, scale = ayear){
  var.type <- match.arg(var.type)


  is_null_newdata <- is.null(newdata)
  if(is_null_newdata){
    if(any(class(fit) %in% c("stpm2", "pstpm2"))){
      data <- fit@data
      newdata <- data.frame(arbritary_var = 0)
    }else{
      data <- fit$data
    }
  }

  if(is.null(exp.fun)){
    #The time points for the expected survival
    times <- seq(0, tau + 1, by = 0.1)

    #Extract expected survival function
    if(is_null_newdata){
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = scale, times = times * scale)))
    }else{
      expected <- vector("list", nrow(newdata))
      for(i in 1:length(expected)){
        expected[[i]] <- do.call("survexp",
                                 list(formula = ~ 1, rmap = substitute(rmap),
                                      data = newdata[i, ], ratetable = ratetable,
                                      scale = scale, times = times * scale))
      }
    }
    exp.fun <- lapply(1:length(expected), function(i){
      smooth.obj <- smooth.spline(x = expected[[i]]$time, y = expected[[i]]$surv, all.knots = T)
      function(time) predict(smooth.obj, x = time)$y
    })
  }



  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    ci.new <- F
    f <- function(time, q) calc.Crude(fit, time = time, type = "condother",
                                      var.type = "n", newdata = newdata[i,,drop = F], tau = tau,
                                      exp.fun = exp.fun[i], reverse = reverse, link = "identity")[[1]]$Estimate - q
    lower <- 0
    if(f(lower, q = q) > 0 & f(max.time, q = q) < 0){
      uni <- rootSolve::uniroot.all(f, lower = lower, upper = max.time, q = q)
    }else{
      if(f(lower, q = q) <= 0){
        uni <- 0
        ci.new <- T
      } else if(f(max.time, q = q) >= 0){
        uni <- NA
        ci.new <- T
      }
    }
    if(var.type %in% c("ci", "se")){
      if(!ci.new){
        gr <- numDeriv::grad(f, x = uni, q = 0)
        VAR <- calc.Crude(fit, time = uni, exp.fun = exp.fun[i], newdata = newdata[i,,drop = F],
                          tau = tau, var.type = "se", type = "condother", link = "identity",
                          reverse = reverse)[[1]]$SE ^ 2
        SE <- sqrt(gr ^ (-2) * VAR / (uni ^ 2))
        upper <- log(uni) + SE * qnorm(0.975)
        lower <- log(uni) - SE * qnorm(0.975)
        df <- data.frame(Estimate = uni, SE = SE * uni, lower.ci = exp(lower), upper.ci = exp(upper))
      } else {
        df <- data.frame(Estimate = uni, SE = NA, lower.ci = NA, upper.ci = NA)
      }
      #if(var.type == "ci") subset(df, select = -SE) else subset(df, select = -c(lower.ci, upper.ci))
      df
    } else{
      data.frame(Estimate = uni)
    }
  })

  do.call(rbind, ests)
}
