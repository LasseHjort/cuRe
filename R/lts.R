#' Long term survival predictions
#'
#' Function for computing survival estimates using a relative survival model and the expected background survival.
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fcm}, \code{gfcm}, \code{stpm2},
#' \code{pstpm2}, and \code{cm}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and the maximum follow-up time is selected.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param ratetable Object of class \code{ratetable} used to compute the general population survival.
#' Default is \code{survexp.dk}.
#' @param exp.fun Object of class \code{list} containing functions for the expected survival
#' of each row in \code{newdata}. If not specified, the function computes the expected
#' survival using the \code{survival::survexp} function and smoothing by \code{smooth.spline}.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package if \code{exp.fun = NULL}.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @return A object of class \code{lts} containing the loss of lifetime estiamtes
#' of each individual in \code{newdata}.
#' @export
#' @example inst/predict.lts.ex.R


lts <- function(fit, type = c("surv", "hazard", "cumhaz", "loghaz"),
                newdata = NULL, time = NULL, var.type = c("ci", "se", "n"),
                exp.fun = NULL, ratetable = survexp.dk, rmap){

  var.type <- match.arg(var.type)
  type <- match.arg(type)

  if(is.null(time)){
    if(any(class(fit) %in% c("stpm2", "pstpm2"))){
      time <- seq(1e-05, max(fit@data[[fit@timeVar]]), length.out = 100)
    } else {
      time <- seq(1e-05, max(fit$time), length.out = 100)
    }
  }

  if(is.null(exp.fun)){
    #The time points for the expected survival
    times <- seq(0, max(time), by = 0.1)

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
    exp.fun <- lapply(1:length(expected), function(i){
      smooth.obj <- smooth.spline(x = expected[[i]]$time, y = expected[[i]]$surv, all.knots = T)
      function(time) predict(smooth.obj, x = time)$y
    })
  }

  #Extract relative survival function
  if(any(class(fit) %in% c("stpm2", "pstpm2"))){
    if(class(fit) == "stpm2"){
      response_name <- all.vars(fit@call.formula)[1]
    }else{
      response_name <- all.vars(fit@fullformula)[1]
    }

    fit_tmp <- fit
    rel_surv <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        #res <- rep(NA, length(t))
        fit_tmp@fullcoef <- pars
        #wh <- which(t != 0)
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        #res[wh] <-
        as.numeric(predict(fit_tmp, newdata = newdata_tmp, type = "surv"))
        #res[-wh] <- 1
        #res
      }
    })

    excess_haz <- lapply(1:length(exp.fun), function(i){
      function(t, pars){
        fit_tmp@fullcoef <- pars
        suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], t))
        names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
        as.numeric(predict(fit_tmp, newdata = newdata_tmp, type = "hazard"))
      }
    })
    model.params <- fit@fullcoef
    cov <- fit@vcov
  }else{
    if ("cuRe" %in% class(fit)) {
      rel_surv <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, var.type = "n")[[1]]$Estimate
      })

      excess_haz <- lapply(1:length(exp.fun), function(i){
        function(t, pars) predict(fit, newdata = newdata[i,, drop = F],
                                  time = t, pars = pars, type = "hazard",
                                  var.type = "n")[[1]]$Estimate
      })
    }
    model.params <- c(unlist(fit$coefs), fit$coefs.spline)
    cov <- fit$covariance
  }

  expected_haz <- lapply(1:length(exp.fun), function(i){
    cum_haz_smooth <- function(t) -log(exp.fun[[i]](t))
    function(t, pars) numDeriv::grad(func = cum_haz_smooth, t)
  })

  #pred <- predict(fit, newdata = newdata, type = type, var.type = var.type)

  res <- vector("list", length(exp.fun))
  for(i in 1:length(exp.fun)){
    Est <- if(type == "surv"){
      exp.fun[[i]](time) * rel_surv[[i]](time, model.params)
      # gr <- jacobian(rel_surv[[i]], x = model.params, t = time)
      # SE <- sqrt(apply(gr, 1, function(x) x %*% cov %*% x))
    } else if(type == "hazard"){
      expected_haz[[i]](time) + excess_haz[[i]](time, model.params)
    } else if(type == "surv"){
      -log(exp.fun[[i]](time) * rel_surv[[i]](time, model.params))
    } else if(type == "loghazard"){
      log(expected_haz[[i]](time) + excess_haz[[i]](time, model.params))
    }

    res[[i]] <- data.frame(Estimate = Est)
  }
  attributes(res) <- list(time = time, type = type)
  class(res) <- "lts"
  res
}
