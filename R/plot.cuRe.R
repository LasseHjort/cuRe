#' Plot function for Flexible mixture cure model
#'
#' Plot function associated with the class \code{fmcm}
#'
#' @param fit Object of class \code{cuRe}.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Character. Defines the function to plot.
#' Possible values are \code{relsurv} (default) for the relative survival, \code{ehaz} for the excess hazard,
#' \code{probcure} for the conditional probability of being cured, and \code{survuncured}
#' for the survival of the uncured patients.
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

# plot.cuRe <- function(fit, newdata = NULL, type = c("relsurv", "ehaz", "probcure", "survuncured"),
#                       time = NULL, xlim = NULL, ylim = c(0, 1),
#                       xlab = "Time", ylab = NULL, col = 1, ci = NULL,
#                       non.parametric = F, col.non.para = 2,
#                       include.knots = F, add = F, ...){
#
#   type <- match.arg(type)
#   if(is.null(ylab)){
#     ylab <- switch(type,
#                    relsurv = "Relative survival",
#                    ehaz = "Excess hazard",
#                    probcure = "Conditional probability of cure",
#                    survuncured = "Disease specific survival of the uncured")
#   }
#
#   if(is.null(ci)){
#     if(add){
#       ci <- "n"
#     } else {
#       ci <- ifelse(fit$ci, "ci", "n")
#     }
#   }
#
#   if(length(col) == 1 & !is.null(newdata)){
#     col <- rep(col, nrow(newdata))
#   }
#   if(is.null(time)){
#     if(is.null(xlim)){
#       xlim <- c(0, max(fit$times))
#       time <- seq(xlim[1], xlim[2], length.out = 100)
#     }
#   }else{
#     xlim <- range(time)
#   }
#
#   predict_rs <- predict(fit, newdata, time, type = type, var.type = ci)
#   nr.samples <- length(predict_rs)
#   if(type == "ehaz"){
#     ylim <- range(unlist(lapply(predict_rs, function(x) x[,-2])), na.rm = T, finite = T)
#   }
#
#   for(i in 1:nr.samples){
#     if(i == 1 & !add){
#       plot(Estimate ~ time, data = predict_rs[[i]], type = "l", ylim = ylim, xlim = xlim,
#            xlab = xlab, ylab = ylab, col = col[i], ...)
#     }else{
#       lines(Estimate ~ time, data = predict_rs[[i]], type = "l", col = col[i], ...)
#     }
#     if(ci == "ci"){
#       lines(ci.upper ~ time, data = predict_rs[[i]], type = "l", col = col[i], lty = 2, ...)
#       lines(ci.lower ~ time, data = predict_rs[[i]], type = "l", col = col[i], lty = 2, ...)
#     }
#   }
#   if(non.parametric){
#     rsfit <- relsurv::rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
#                      data = fit$data, ratetable = survexp.dk, method = "ederer2")
#     rsfit$time <- rsfit$time / year
#     lines(rsfit$surv ~ rsfit$time, type = "s", col = col.non.para, ...)
#   }
#   if(include.knots){
#     abline(v = exp(fit$knots), lty = 2)
#   }
# }




plot.cuRe <- function(object, newdata = NULL, type = c("surv", "probcure", "survuncured", "hazarduncured",
                                                       "cumhazuncured", "densityuncured", "failuncured",
                                                       "oddsuncured", "loghazarduncured", "hazard",
                                                       "density", "fail", "loghazard", "odds", "cumhaz"),
                      time = NULL, xlim = NULL, ylim = c(0, 1),
                      xlab = "Time", ylab = NULL, col = 1, ci = NULL,
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

  if(is.null(ci)){
    if(add){
      ci <- "n"
    } else {
      ci <- ifelse(object$ci, "ci", "n")
    }
  }

  pred <- predict(object, newdata, time, type = type,
                  var.type = ci)

  if(type == "hazard"){
    ylim <- range(unlist(lapply(pred, function(x) x[,-2])), na.rm = T, finite = T)
  }

  nr.samples <- length(pred)
  for(i in 1:nr.samples){
    if(i == 1 & !add){
      plot(Estimate ~ time, data = pred[[i]], type = "l", ylim = ylim, xlim = xlim,
           xlab = xlab, ylab = ylab, col = col[i], ...)
    }else{
      lines(Estimate ~ time, data = pred[[i]], type = "l", col = col[i], ...)
    }
    if(ci == "ci"){
      lines(upper ~ time, data = pred[[i]], type = "l", col = col[i], lty = 2, ...)
      lines(lower ~ time, data = pred[[i]], type = "l", col = col[i], lty = 2, ...)
    }
  }
}

