#' Plot function for Flexible mixture cure model
#'
#' Plot function associated with the classes \code{gfcm} and \code{cm}
#'
#' @param object Object of class \code{cuRe}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Character. Defines the desired scale to plot. See ?predict.gfcm for possible values.
#' @param time Optional time points at which to compute predictions.
#' This argument is not used if type is \code{curerate}.
#' @param xlim Limits of the x-axis
#' @param ylim Limits of the y-axis.
#' @param xlab Label of the x-axis. Default is "Time".
#' @param ylab Label of the y-axis. If \code{NULL}, the function uses its default y-labels
#' depending on \code{object$type} and \code{object$excess}.
#' @param col Colour of each line.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are added to the plot.
#' @param add Loglca. If \code{TRUE} the curve is added to the existing plot.
#' @param ... Further arguments passed to \code{plot} and \code{lines}.
#' @return A plot containing the predictions of each observation in \code{newdata}.
#' @export

plot.cuRe <- function(object, newdata = NULL, type = c("surv", "probcure", "survuncured", "hazarduncured",
                                                       "cumhazuncured", "densityuncured", "failuncured",
                                                       "oddsuncured", "loghazarduncured", "hazard",
                                                       "density", "fail", "loghazard", "odds", "cumhaz"),
                      time = NULL, xlim = NULL, ylim = c(0, 1),
                      xlab = "Time", ylab = NULL, col = 1, ci = T,
                      add = F, ...){

  type <- match.arg(type)

  if(is.null(ylab)){
    if(!object$excess){
      ylab <- switch(type,
                     linkS = "Linear predictor", probcure = "Probability of cure",
                     survuncured = "Survival of the uncured", hazarduncured = "Hazard of the uncured",
                     cumhazuncured = "Cumulative hazard of the uncured",
                     densityuncured = "Density of the uncured", failuncured = "Distribution of the uncured",
                     oddsuncured = "Odds survival of the uncured", loghazarduncured = "Log-hazard of the uncured",
                     surv = "Survival probability", hazard = "Hazard", density = "Density", fail = "Distribution",
                     loghazard = "Log-hazard", odds = "Odds", cumhaz = "Cumulative incidence")
    } else {
      ylab <- switch(type,
                     linkS = "Linear predictor", probcure = "Probability of cure",
                     survuncured = "Relative survival of the uncured",
                     hazarduncured = "Excess hazard of the uncured",
                     cumhazuncured = "Cumulative excess hazard of the uncured",
                     densityuncured = "Excess density of the uncured",
                     failuncured = "Net distribution of the uncured",
                     oddsuncured = "Net odds survival of the uncured",
                     loghazarduncured = "Log-excess hazard of the uncured",
                     surv = "Relative survival", hazard = "Excess hazard", density = "Excess density",
                     fail = "Net distribution", loghazard = "Log-excess hazard",
                     odds = "Net odds", cumhaz = "Cumulative excess hazard")
    }
  }

  if(length(col) == 1 & !is.null(newdata)){
    col <- rep(col, nrow(newdata))
  }

  if(is.null(time)){
    if(is.null(xlim)){
      xlim <- c(1e-05, max(object$time))
      time <- seq(xlim[1], xlim[2], length.out = 100)
    }
  }else{
    xlim <- range(time)
  }

  if(!(ci & !add & object$ci)){
    ci <- F
  }

  pred <- predict(object, newdata = newdata, time = time, type = type,
                  var.type = if(ci) "ci" else "n")

  if(type == "hazard"){
    ylim <- range(unlist(lapply(pred, function(x) x)), na.rm = T, finite = T)
  }

  nr.samples <- length(pred)
  for(i in 1:nr.samples){
    if(i == 1 & !add){
      plot(Estimate ~ time, data = pred[[i]], type = "l", ylim = ylim, xlim = xlim,
           xlab = xlab, ylab = ylab, col = col[i], ...)
    }else{
      lines(Estimate ~ time, data = pred[[i]], type = "l", col = col[i], ...)
    }
    if(ci){
      lines(upper ~ time, data = pred[[i]], type = "l", col = col[i], lty = 2, ...)
      lines(lower ~ time, data = pred[[i]], type = "l", col = col[i], lty = 2, ...)
    }
  }
}

