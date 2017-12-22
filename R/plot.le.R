#' Plot the loss of lifetime function
#'
#' Plot function for the class \code{le}
#'
#' @param obj Object of class \code{le}.
#' @param ylim Limit of the y-axis. If empty, these are chosen occording to the extremes in \code{obj}.
#' @param xlim Limit of x-axis. If empty, these are chosen according to the time extremes of \code{obj}.
#' @param ci Logical indicating wether confidence intervals should be included in the plot (default is TRUE).
#' @param col Numeric or character indicating the colours of the curves.
#' @param ylab Label to be written on the y-axis.
#' @param xlab Label to be written on the x-axis.
#' @param add Logical indicating whether to add to current plot window (default is FALSE).
#' @param ... Further argument passed to \code{plot}.
#' @export

plot.le <- function(obj, ylim = NULL, xlim = NULL, ci = T, col = 1,
                    ylab = NULL, xlab = "Time", add = F, ...){
  if(is.null(ylab)){
    ylab <- switch(obj$type,
                   ll = "Loss of lifetime",
                   mrl = "Mean residual lifetime")
  }

  ci <- ci & obj$ci

  if(length(col) == 1){
    col <- rep(col, length(obj$Ests))
  }
  if(is.null(ylim)){
    if(ci){
      ylim <- range(unlist(lapply(obj$Ests, function(x) x[, c("lower.ci", "upper.ci")])))
    }else{
      ylim <- range(unlist(lapply(obj$Ests, function(x) x[, obj$type])))
    }
  }
  for(i in 1:length(obj$Ests)){
    if(i == 1 & !add){
      plot(obj$Ests[[i]][, obj$type] ~ obj$time, ylim = ylim, xlim = xlim,
           type = "l", col = col[i], xlab = xlab, ylab = ylab)
    }else{
      lines(obj$Ests[[i]][, obj$type] ~ obj$time, data = obj$Est[[i]], col = col[i])
    }
    if(ci){
      lines(lower.ci ~ obj$time, data = obj$Ests[[i]], lty = 2, col = col[i])
      lines(upper.ci ~ obj$time, data = obj$Ests[[i]], lty = 2, col = col[i])
    }
  }
}
