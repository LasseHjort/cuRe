#' Plot the long term survival
#'
#' Plot function for the class \code{lts}.
#'
#' @param obj Object of class \code{lts}.
#' @param ylim Limit of the y-axis. If empty, these are chosen occording to the extremes in \code{obj}.
#' @param xlim Limit of x-axis. If empty, these are chosen according to the time extremes of \code{obj}.
#' @param ci Logical indicating wether confidence intervals should be included in the plot (default is TRUE).
#' @param col Numeric or character indicating the colours of the curves.
#' @param ylab Label to be written on the y-axis.
#' @param xlab Label to be written on the x-axis.
#' @param add Logical indicating whether to add to current plot window (default is FALSE).
#' @param ... Further argument passed to \code{plot} and \code{lines}.
#' @export

plot.lts <- function(obj, ylim = NULL, xlim = NULL, ci = F, col = 1,
                    ylab = NULL, xlab = "Time", add = F, ...){
  time <- attr(obj, "time")
  if(is.null(ylab)){
    ylab <- switch(attr(obj, "type"),
                   surv = "Survival probability",
                   hazard = "Hazard",
                   cumhaz = "Cumulative hazard",
                   loghaz = "Log-hazard")
  }
  if(length(col) == 1){
    col <- rep(col, length(obj))
  }
  if(is.null(ylim)){
    ylim <- c(0, 1)
  }
  for(i in 1:length(obj)){
    if(i == 1 & !add){
      plot(obj[[i]]$Estimate ~ time, ylim = ylim, xlim = xlim,
           type = "l", col = col[i], xlab = xlab, ylab = ylab, ...)
    }else{
      lines(obj[[i]]$Estimate ~ time, col = col[i], ...)
    }
    if(ci){
      lines(lower.ci ~ time, data = obj[[i]], lty = 2, col = col[i], ...)
      lines(upper.ci ~ time, data = obj[[i]], lty = 2, col = col[i], ...)
    }
  }
}
