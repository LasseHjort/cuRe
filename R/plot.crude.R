#' Plot crude event probabilities
#'
#' Plot function for the computed crude event probabilties.
#'
#' @param obj Object of class \code{crude} in which crude probabilities are calculated.
#' @param ylim Limits of y-axis.
#' @param xlim Limits of x-axis.
#' @param col Colour of each curve.
#' @param ylab Label of the y-axis. If \code{NULL}, the function uses its default labels depending on \code{obj$type}.
#' @param xlab Label of the x-axis (default is "Time").
#' @param add Logical indicating wether the cruves should be added to the current plot window (default is \code{FALSE}).
#' @param ci Logical denoting whether confidence intervals should be plotted (default is \code{TRUE}).
#' @param ... Further arguments passed to \code{plot} and \code{lines}.
#' @export

plot.crude <- function(obj, ylim = c(0, 1), xlim = NULL, ci = T,
                       col = 1, ylab = NULL, xlab = "Time", add = F, ...){

  attr <- attributes(obj)
  if(is.null(ylab)){
    ylab <- switch(attr$type, cancer = "Cumulative incidence of disease-related death",
                   other = "Cumulative incidence of non-disease-related death",
                   condother = "Probability of non-disease-related death")
    if(attr$reverse) ylab <- "Probability of disease-related death"
  }
  if(is.null(xlim)) xlim <- range(attr$time)

  if(length(col) == 1){
    col <- rep(col, length(obj))
  }
  if(ci & !attr$ci){
    ci <- FALSE
  }

  for(i in 1:length(obj)){
    if(i == 1 & !add){
      plot(Estimate ~ attr$time, data = obj[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], ylab = ylab, xlab = xlab, ...)
    }else{
      lines(Estimate ~ attr$time, data = obj[[i]], col = col[i], ...)
    }
    if(ci){
      lines(lower.ci ~ attr$time, data = obj[[i]], lty = 2, col = col[i], ...)
      lines(upper.ci ~ attr$time, data = obj[[i]], lty = 2, col = col[i], ...)
    }
  }
}
