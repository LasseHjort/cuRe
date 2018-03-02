#' Compute the time to statistical cure using the conditional probability of cure
#'
#' The following function estimates the time to statistical cure using the
#' conditional probability of cure
#' @param fit Fitted model to do predictions from. Possible classes are
#' \code{fmc}, \code{cm}, \code{stpm2}, and \code{pstpm2}.
#' @param q Threshold to estimate statistical cure according to.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param max.time Upper boundary of the interval [0, \code{max.time}] in which to search for solution.
#' @param ci Logical. If \code{TRUE} (default), confidence intervals are computed.
#' @param reverse Logical. Wether to use 1 - prob (default) or prob as measure.
#' @return The estimated cure points.
#' @export

calc.cure.quantile <- function(fit, q = 0.95, newdata = NULL, max.time = 20, ci = TRUE, reverse = TRUE){

  if("gfcm" %in% class(fit)){
    pred <- function(time, newdata, var.type = "n"){
      predict(fit, time = time, type = "probcure", newdata = newdata[i,,drop = F],
              var.type = var.type, link = "I")[[1]]
    }
  } else if ("cm" %in% class(fit)){
    pred <- function(time, newdata, var.type = "n"){
      predict(fit, time = time, type = "probcure", newdata = newdata[i,,drop = F],
              var.type = var.type, link = "identity")[[1]]
    }
  } else if("stpm2" %in% class(fit)){
    pred <- function(time, newdata, var.type = "n"){
      times <- c(time, max(fit@data[[fit@timeVar]]))
      p <- predict.stpm2(fit, newdata = data.frame(FU_years = times), type = "probcure",
                         exposed = function(x) x[which.max(x$FU_years),, drop = F],
                         se.fit = ifelse(var.type == "n", FALSE, TRUE), link = "I")
      col_names = if(var.type == "n") "Estimate" else c("Estimate", "SE")
      D <- as.data.frame(p)
      colnames(D) <- col_names
      D[-nrow(D),, drop = F]
    }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    if(reverse){
      f <- function(time, q) 1 - pred(time, newdata = newdata[i,,drop = F])$Estimate - q
    }else {
      f <- function(time, q) pred(time, newdata = newdata[i,,drop = F])$Estimate - q
    }
    uni <- rootSolve::uniroot.all(f, lower = 1e-05, upper = max.time, q = q)
    if(ci){
      gr <- grad(f, x = uni, q = 0)
      VAR <- pred(time = uni, newdata = newdata[i,,drop = F], var.type = "se")$SE ^ 2
      VAR2 <- gr ^ (-2) * VAR / (uni ^ 2)
      upper <- log(uni) + sqrt(VAR2) * qnorm(0.975)
      lower <- log(uni) - sqrt(VAR2) * qnorm(0.975)
      data.frame(Est = uni, var = VAR2 * uni ^ 2, lower.ci = exp(lower), upper.ci = exp(upper))
    } else{
      data.frame(Est = uni)
    }
  })

  do.call(rbind, ests)
}

