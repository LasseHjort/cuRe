#' Long term survival predictions
#'
#' Function for computing loss of lifetime estimates from an estimated relative survival model
#'
#' @param fit Fitted model to do predictions from. Possible classes are \code{fcm}, \code{gfcm}, \code{stpm2},
#' \code{pstpm2}, and \code{cm}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions. If empty, a grid of 100 time points between 0
#' and the maximum follow-up time is selected.
#' @param ci Logical indicating whether confidence intervals should be computed.
#' @param ratetable Object of class \code{ratetable} to compute the general population survival from.
#' @param expected Object of class \code{list} containing objects of class \code{survexp}
#' denoting the expected survival of each row in newdata. If not specified, the function computes the expected
#' survival.
#' @param rmap List to be passed to \code{survexp} from the \code{survival} package.
#' Detailed documentation on this argument can be found by \code{?survexp}.
#' @return A object of class \code{lts} containing the loss of lifetime estiamtes
#' of each individual in \code{newdata}.
#' @export
#' @example inst/predict.lts.ex.R


lts <- function(fit, newdata = NULL, time = NULL, ci = T, expected = NULL, ratetable = survexp.dk,
                        rmap = rmap){

  if(!is.null(time)) time <- sort(time)
  if(is.null(expected)){
    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
        if(is.null(time)) time <- seq(0, max(eval(fit@timeExpr, fit@data)), length.out = 100)
      }else{
        data <- fit$data
        if(is.null(time)) time <- seq(0, max(fit$times), length.out = 100)
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = time * ayear)))

    }else{
      expected <- vector("list", nrow(newdata))
      for(i in 1:length(expected)){
        expected[[i]] <- do.call("survexp",
                                 list(formula = ~ 1, rmap = substitute(rmap),
                                      data = newdata[i, ], ratetable = ratetable,
                                      scale = ayear, times = time * ayear))
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

    rel_surv <- lapply(1:length(expected), function(i){
      wh <- which(time != 0)
      suppressWarnings(newdata_tmp <- cbind(newdata[i,,drop = F], time[wh]))
      names(newdata_tmp)[ncol(newdata_tmp)] <- response_name
      pred <- predict(fit, newdata = newdata_tmp, type = "surv", se.fit = ci)
      if(ci){
        D <- rbind(c(1,1,1), pred)
        names(D) <- c("Est", "ci.lower", "ci.upper")
      }else{
        D <- rbind(1, pred)
        names(D) <- "Est"
      }
      D
    })
  }else{
    rel_surv <- lapply(1:length(expected), function(i){
      predict(fit, newdata = newdata[i,, drop = F],
              time = time, ci = ci)$res[[1]]
    })
  }

  survs <- vector("list", length(expected))
  for(i in 1:length(expected)){
    Est <- rel_surv[[i]]$Est * expected[[i]]$surv
    if(ci){
      ci.lower <- rel_surv[[i]]$ci.lower * expected[[i]]$surv
      ci.upper <- rel_surv[[i]]$ci.upper * expected[[i]]$surv
      survs[[i]] <- data.frame(Est = Est, lower.ci = ci.lower, upper.ci = ci.upper)
    }else{
      survs[[i]] <- data.frame(Est = Est)
    }
  }

  all_res <- list(Ests = survs, time = time)
  class(all_res) <- "lts"
  all_res
}
