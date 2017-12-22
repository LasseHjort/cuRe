#' Plot function for Flexible mixture cure model
#'
#' Plot function associated with the class \code{fmcm}
#'
#' @param fit Object of class \code{cuRe}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, and \code{probcure} for the
#' conditional probability of being cured.
#' @param time Optional time points at which to compute predictions.
#' This argument is not used if type is \code{curerate}.
#' @param xlim Limits of the x-axis
#' @param ylim Limits of the y-axis.
#' @param xlabel Label of the x-axis. Default is "Time".
#' @param ylabel Label of the y-axis. Default is "Relative survival" if \code{type = relsurv},
#' "Excess hazard" if \code{type = ehaz},
#' "Conditional probability of cure" if \code{type = probcure},
#' "Disease-specific survival of the uncured" if \code{type = survuncured}.
#' @param col Colour of each line.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are added.
#' @param non.parametric Logical. If \code{TRUE}, a non-parametic
#' estimate (Ederer II method) is added to the plot (default is FALSE).
#' @param col.non.para Colour of the non-parametric estimate.
#' @param include.knots Logical. If \code{TRUE} vertical lines are added at
#' each knot of the baseline spline function (default is \code{FALSE}).
#' @return A plot containing the predictions of each observation in \code{newdata}.
#' @export
#' @import relsurv

plot.fcim <- function(fit, newdata = NULL, type = "cuminc",
                      time = NULL, cause = NULL, xlim = NULL, ylim = c(0, 1),
                      xlab = "Time", ylab = NULL, col = 1, ci = NULL,
                      non.parametric = F, col.non.para = 2,
                      include.knots = F, add = F, ...){

  if(is.null(ylab)){
    ylab <- switch(type,
                   cuminc = "Cumulative incidence",
                   subhaz = "Subdistribution hazard")
  }

  if(is.null(ci)) ci <- fit$ci

  if(length(col) == 1 & !is.null(newdata)){
    col <- rep(col, nrow(newdata))
  }

  if(is.null(time)){
    if(is.null(xlim)){
      xlim <- c(0, max(fit$times))
      time <- seq(xlim[1], xlim[2], length.out = 100)
    }
  }else{
    xlim <- range(time)
  }

  predict_rs <- predict(fit, newdata, time, type = type, ci = ci)

  if(is.null(cause)) cause <- 1:length(predict_rs$res[[1]])
  nr.samples <- length(predict_rs$res)

  if(type == "subhaz"){
    ylim <- c()
    for(j in 1:length(cause)){
      ylim <- c(ylim, range(unlist(lapply(predict_rs$res, function(x) x[[j]][,-2])), na.rm = T, finite = T))
    }
    ylim <- range(ylim)
  }

  for(i in 1:nr.samples){
    for(j in 1:length(cause)){
      if(i == 1 & !add){
        plot(Est ~ time, data = predict_rs$res[[i]][[j]], type = "l", ylim = ylim, xlim = xlim,
             xlab = xlab, ylab = ylab, col = col[i], ...)
      }else{
        lines(Est ~ time, data = predict_rs$res[[i]][[j]], type = "l", col = col[i], ...)
      }
      if(ci){
        lines(ci.upper ~ time, data = predict_rs$res[[i]][[j]], type = "l", col = col[i], lty = 2, ...)
        lines(ci.lower ~ time, data = predict_rs$res[[i]][[j]], type = "l", col = col[i], lty = 2, ...)
      }
      add <- TRUE
    }
  }
  if(include.knots){
    for(j in 1:length(cause)){
      abline(v = exp(fit$knots[[j]]), lty = 2, col = j)
    }
  }
}
