#' Plot crude event probabilities
#'
#' Plot function for the computed crude event probabilties.
#'
#' @param obj Object of class \code{crude} in which crude probabilities are calculated
#' @export

plot.crude <- function(obj, ylim = c(0, 1), xlim = NULL, ci = T, col = 1, ylab = NULL, xlab = "Time"){
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
    if(i == 1){
      plot(prob ~ obj$time, data = obj$prob[[i]], ylim = ylim, xlim = xlim,
           type = "l", col = col[i], ylab = ylab, xlab = xlab)
    }else{
      lines(prob ~ obj$time, data = obj$prob[[i]], col = col[i])
    }
    if(ci){
      lines(lower.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i])
      lines(upper.ci ~ obj$time, data = obj$prob[[i]], lty = 2, col = col[i])
    }
  }
}
