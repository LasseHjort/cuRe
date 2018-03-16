#' Plot the loss of lifetime function
#'
#' Plot function for the class \code{le}
#'
#' @param obj Object of class \code{le}.
#' @param ylim Limit of the y-axis.
#' @param xlim Limit of x-axis.
#' @param ci Logical. If \code{TRUE}, confidence intervals are added to the plot.
#' @param col Numeric or character indicating the colours of the curves.
#' @param ylab Label to be written on the y-axis.
#' @param xlab Label to be written on the x-axis.
#' @param add Logical. If \code{TRUE}, the curve is added to the current plot window.
#' @param ... Further argument passed to \code{plot} and \code{lines}.
#' @export

plot.le <- function(obj, ylim = NULL, xlim = NULL, ci = T, col = 1,
                    ylab = NULL, xlab = "Time", add = F, ...){
  att <- attributes(obj)

  if(is.null(ylab)){
    ylab <- switch(att$type,
                   ll = "Loss of lifetime",
                   mrl = "Mean residual lifetime")
  }

  ci <- ci & att$var.type == "ci"

  if(length(col) == 1){
    col <- rep(col, length(obj))
  }
  if(is.null(ylim)){
    if(ci){
      ylim <- range(unlist(lapply(obj, function(x) x[, c("lower.ci", "upper.ci")])))
    }else{
      ylim <- range(unlist(lapply(obj, function(x) x[, obj$type])))
    }
  }
  for(i in 1:length(obj)){
    if(i == 1 & !add){
      plot(Estimate ~ att$time, data = obj[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], xlab = xlab, ylab = ylab, ...)
    }else{
      lines(Estimate ~ att$time, data = obj[[i]], col = col[i], ...)
    }
    if(ci){
      lines(lower.ci ~ att$time, data = obj[[i]], lty = 2, col = col[i], ...)
      lines(upper.ci ~ att$time, data = obj[[i]], lty = 2, col = col[i], ...)
    }
  }
}
