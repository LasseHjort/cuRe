#' Plot crude event probabilities
#'
#' Plot function for the computed crude event probabilties.
#'
#' @param obj Object of class \code{crude} in which crude probabilities are calculated.
#' @param ylim Limits of y-axis.
#' @param xlim Limits of x-axis.
#' @param col Colour of each curve.
#' @param ylab Label of the y-axis. If type is \code{cancer}, the label is "Cumulative incidence of
#' cancer related death", if \code{other}, the label is "Cumulative incidence of death from other causes than cancer",
#' and if type is \code{othertime}, the label is "probability of eventually dying from other causes than cancer".
#' @param xlab Label of the x-axis (default is "Time").
#' @param add Logical indicating wether the cruves should be added to the current plot window (default is FALSE).
#' @param ci Logical denoting whether confidence intervals should be plotted (default is TRUE).
#' @param ... Further arguments passed to \code{plot} and \code{lines}.
#' @export

plot.crude <- function(obj, ylim = c(0, 1), xlim = NULL, ci = T,
                       col = 1, ylab = NULL, xlab = "Time", add = F, ...){
  if(is.null(ylab)){
    ylab <- switch(obj$type, cancer = "Cumulative incidence of cancer related death",
                   other = "Cumulative incidence of death from other causes than cancer",
                   othertime = "Probability of eventually dying from other causes than cancer")
  }
  if(is.null(xlim)) xlim <- range(obj$time)

  if(length(col) == 1){
    col <- rep(col, length(obj$prob))
  }
  if(ci & !obj$ci){
    ci <- FALSE
  }

  for(i in 1:length(obj$prob)){
    if(i == 1 & !add){
      plot(prob ~ obj$time, data = obj$prob[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], ylab = ylab, xlab = xlab, ...)
    }else{
      lines(prob ~ obj$time, data = obj$prob[[i]], col = col[i], ...)
    }
    if(ci){
      lines(lower.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i], ...)
      lines(upper.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i], ...)
    }
  }
}
