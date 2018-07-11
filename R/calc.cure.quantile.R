#' Compute the time to statistical cure using the conditional probability of cure
#'
#' The following function estimates the time to statistical cure using the
#' conditional probability of cure.
#' @param fit Fitted model to do predictions from. Possible classes are
#' \code{gfcm}, \code{cm}, and \code{stpm2}.
#' @param q Threshold to estimate statistical cure according to.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param max.time Upper boundary of the interval [0, \code{max.time}] in which to search for solution (see details).
#' Default is 20.
#' @param var.type Character. Possible values are "\code{ci}" (default) for confidence intervals,
#' "\code{se}" for standard errors, and "\code{n}" for neither.
#' @param reverse Logical. Whether to use the conditional probability of not being cured (default) or
#' the conditional probability of cure.
#' @param bdr.knot Time point from which cure is assumed. Only relevant for class \code{stpm2}.
#' @return The estimated cure point.
#' @details The cure point is calculated as the time point at which the conditional probability of disease-related
#' death reaches the threshold, \code{q}. If \code{q} is not reached within \code{max.time}, no solution is reported.
#' @example inst/calc.cure.quantile.ex.R
#' @export

calc.cure.quantile <- function(fit, q = 0.05, newdata = NULL, max.time = 20, var.type = c("ci", "n"),
                               reverse = TRUE, bdr.knot = NULL){
  var.type <- match.arg(var.type)

  if(any(c("gfcm", "cm") %in% class(fit))){
    pred <- function(time, newdata, var.type = "n"){
      predict(fit, time = time, type = "probcure", newdata = newdata,
              var.type = var.type, link = "I")[[1]]
    }
  } else if("stpm2" %in% class(fit)){
    if(is.null(bdr.knot)){
      bdr.knot <- exp(max(as.list(as.list(attr(fit@terms, "predvars"))[[3]])$Boundary.knots))
    }
    exposed <- function(newdata){
      newdata[[fit@timeVar]] <- bdr.knot
      newdata
    }

    pred <- function(time, newdata, var.type = "n"){
      data.x <- data.frame(time)
      names(data.x) <- fit@timeVar
      newdata[[fit@timeVar]] <- NULL
      newdata <- merge(newdata, data.x)
      colnames(newdata)[ncol(newdata)] <- fit@timeVar
      p <- predict(fit, newdata = data.frame(FU_years = time), type = "probcure",
                   exposed = exposed, se.fit = ifelse(var.type == "n", FALSE, TRUE), link = "I")
      if(var.type != "n"){
        p$SE <- (p$upper - p$Estimate) / qnorm(0.975)
        return(p)
      } else {
        D <- as.data.frame(p)
        colnames(D) <- "Estimate"
        return(D)
      }
    }

    # pred <- function(time, newdata, var.type = "n"){
    #   times <- c(time, max(fit@data[[fit@timeVar]]))###Change this
    #   p <- predict.stpm2(fit, newdata = data.frame(FU_years = times), type = "probcure",
    #                      exposed = function(x) x[which.max(x$FU_years),, drop = F],
    #                      se.fit = ifelse(var.type == "n", FALSE, TRUE), link = "I")
    #   col_names = if(var.type == "n") "Estimate" else c("Estimate", "SE")
    #   D <- as.data.frame(p)
    #   colnames(D) <- col_names
    #   D[-nrow(D),, drop = F]
    # }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- vector("list", n.obs)
  for(i in 1:n.obs){
    ci.new <- F
    if(reverse){
      f <- function(time, q) 1 - pred(time, newdata = newdata[i,,drop = F])$Estimate - q
    }else {
      f <- function(time, q) pred(time, newdata = newdata[i,,drop = F])$Estimate - q
    }
    lower <- 1e-05
    if(f(lower, q = q) > 0 & f(max.time, q = q) < 0){
      uni <- rootSolve::uniroot.all(f, lower = lower, upper = max.time, q = q)
    }else{
      if(f(lower, q = q) <= 0){
        uni <- 0
        ci.new <- T
      } else if(f(max.time, q = q) >= 0){
        uni <- NA
        ci.new <- T
      }
    }
    if(var.type %in% c("ci", "se")){
      if(!ci.new){
        gr <- grad(f, x = uni, q = 0)
        VAR <- pred(time = uni, newdata = newdata[i,,drop = F], var.type = "se")$SE ^ 2
        SE <- sqrt(gr ^ (-2) * VAR / (uni ^ 2))
        upper <- log(uni) + SE * qnorm(0.975)
        lower <- log(uni) - SE * qnorm(0.975)
        df <- data.frame(Estimate = uni, SE = SE * uni, lower.ci = exp(lower), upper.ci = exp(upper))
      }else{
        df <- data.frame(Estimate = uni, SE = NA, lower.ci = NA, upper.ci = NA)
      }
      #if(var.type == "ci") subset(df, select = -SE) else subset(df, select = -c(lower.ci, upper.ci))
      ests[[i]] <- df
    } else{
      ests[[i]] <- data.frame(Estimate = uni)
    }
  }

  do.call(rbind, ests)
}

