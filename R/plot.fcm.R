#' Plot function for Flexible mixture cure model
#'
#' Plot function associated with the class \code{fmcm}
#'
#' @param fit Object of class \code{fmcm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, and \code{probcure} for the
#' conditional probability of being cured.
#' @param time Optional time points at which to compute predictions. This argument is not used if type is \code{curerate}.
#' @param ci Logical indicating whether confidence intervals should be computed
#' @param pars Numerical vector containing the parameters values of the model.
#' In general, this argument can be ignored by the user
#' @param non.parametric Logical indicating if a non-parametic estimate should be included in the plot (default is FALSE).
#' @param non.param.arg List of arguments to be passed to function \code{rs.surv} of the \code{relsurv} package.
#' @param include.knots Logical for including the baseline knots of the model.
#' @return A list containing the predictions of each individual in \code{newdata}.
#' @export

plot.cuRe <- function(fit, newdata = NULL, type = "relsurv",
                               time = NULL, ylim = c(0, 1), xlim = NULL,
                               xlab = "Time", ylab = NULL, non.parametric = F,
                               col = 1, col.non.para = 2, ci = T,
                               include.knots = F, add = F, ...){

  if(is.null(ylab)){
    ylab <- switch(type,
                   relsurv = "Relative survival",
                   ehaz = "Excess hazard",
                   probcure = "Conditional probability of cure",
                   crude_prob = "Probability of eventually dying from other causes than cancer",
                   lol = "Loss of lifetime",
                   survuncured = "Disease specific survival of the uncured")
  }

  if(length(col) == 1 & !is.null(newdata)){
    col <- rep(col, nrow(newdata))
  }
  if(is.null(time)){
    if(is.null(xlim)){
      xlim <- c(0, max(fit$data$FU_years))
      time <- seq(xlim[1], xlim[2], length.out = 100)
    }
  }else{
    xlim <- range(time)
  }

  predict_rs <- predict(fit, newdata, time, type = type, ci = ci)
  nr.samples <- length(predict_rs$res)
  if(type == "ehaz"){
    ylim <- range(unlist(lapply(predict_rs$res, function(x) x[,-2])), na.rm = T, finite = T)
  }

  for(i in 1:nr.samples){
    if(i == 1 & !add){
      plot(Est ~ time, data = predict_rs$res[[i]], type = "l", ylim = ylim, xlim = xlim,
           xlab = xlab, ylab = ylab, col = col[i], ...)
    }else{
      lines(Est ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], ...)
    }
    if(ci){
      lines(ci.upper ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], lty = 2, ...)
      lines(ci.lower ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], lty = 2, ...)
    }
  }
  if(non.parametric){
    rsfit <- relsurv::rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                     data = fit$data, ratetable = survexp.dk, method = "ederer2")
    rsfit$time <- rsfit$time / year
    lines(rsfit$surv ~ rsfit$time, type = "s", col = col.non.para, ...)
  }
  if(include.knots){
    abline(v = exp(fit$knots), lty = 2)
  }
}
