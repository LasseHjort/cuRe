predict.fcm2 <- function (object, newdata = NULL,
                          type = c("surv", "curerate", "probcure", "survuncured", "hazarduncured",
                                   "cumhazuncured", "densityuncured", "failuncured", "oddsuncured",
                                   "loghazarduncured", "hazard", "density", "fail",
                                   "loghazard", "odds", "cumhaz"),
                          indi = TRUE, time = NULL, var.type = c("ci", "se", "n"), pars = NULL,
                          link = NULL, var = NULL,
                          keep.attributes = FALSE, use.gr = FALSE, ...)
{
  type <- match.arg(type)
  args <- object$args

  if(!is.null(pars)){
    object$coefs <- pars[1:length(object$coefs)]
    object$coefs.spline <- pars[(length(object$coefs) + 1):length(pars)]
  }

  calcX <- !is.null(newdata)
  if (is.null(newdata)) {
    if(indi){
      vars <- c(all.vars(formula), all.vars(object$logH.formula))
      vars <- vars[!vars %in% c(as.character(object$timeExpr), as.character(object$eventExpr))]
      if(length(vars) != 0){
        stop("'newdata' needs to be specified with option 'indi = TRUE' when covariates are present")
      }
      newdata <- data.frame(x = 1)
      colnames(newdata) <- "(Intercept)"
    } else {
      X <- object$args$X
      XD <- object$args$XD
      X.cr <- object$args$X.cr
      y <- object$args$time
      time <- y
      newdata <- as.data.frame(object$data)
    }
  }

  lpfunc <- function(delta, fit, data, var) {
    data[[var]] <- data[[var]] + delta
    rstpm2:::lpmatrix.lm(fit, data)
  }

  if (is.null(time)) {
    if(indi){
      time <- seq(1e-05, max(object$time), length.out = 100)
    }else{
      time <- eval(object$timeExpr, newdata, parent.frame())
    }
    if(type == "curerate"){
      time <- 1
    }
  }

  if(indi){
    newdata.list <- split(newdata, f = 1:nrow(newdata))
    for(i in 1:length(newdata.list)){
      newdata.list[[i]] <- newdata.list[[i]][rep(1, length(time)),, drop = F]
      newdata.list[[i]][, object$timeVar] <- time
    }

    X <- lapply(1:length(newdata.list), function(i){
      object$transX(rstpm2:::lpmatrix.lm(object$lm.obj, newdata.list[[i]]), newdata.list[[i]])
    })

    XD <- lapply(1:length(newdata.list), function(i){
      XD.tmp <- rstpm2:::grad(lpfunc, 0, object$lm.obj, newdata.list[[i]], object$timeVar)
      object$transXD(matrix(XD.tmp, nrow = nrow(X[[i]])))
    })

    X.cr <- lapply(1:length(newdata.list), function(i){
      model.matrix(object$cr.formula, data = newdata.list[[i]])
    })

  } else {
    if(calcX){
      X <- object$transX(rstpm2:::lpmatrix.lm(object$lm.obj, newdata),
                         newdata)
      XD <- rstpm2:::grad(lpfunc, 0, object$lm.obj, newdata, object$timeVar)
      XD <- object$transXD(matrix(XD, nrow = nrow(X)))
      X.cr <- model.matrix(object$cr.formula, data = newdata)
    }
  }


  var.type <- match.arg(var.type)
  pred <- if (!var.type %in% c("ci", "se")) {
    if(is.list(X)){
      lapply(1:length(X), function(i){
        data.frame(Estimate = local(object, newdata, type, var.link = function(x) x,
                                    X = X[[i]], XD = XD[[i]], X.cr = X.cr[[i]]))
      })
    } else {
      data.frame(Estimate = local(object, newdata, type, var.link = function(x) x,
                                  X = X, XD = XD, X.cr = X.cr))
    }
  } else {
    gd <- NULL
    #beta <- object$coefs.spline
    if (is.null(link)){
      if(!object$excess){
        link <- switch(type, linkS = "I", linkpi = "I", curerate = "cloglog",
                       probcure = "cloglog", survuncured = "cloglog",
                       hazarduncured = "log", cumhazuncured = "log",
                       densityuncured = "log", failuncured = "cloglog",
                       oddsuncured = "cloglog", loghazarduncured = "I",
                       surv = "cloglog", hazard = "log", density = "log", fail = "cloglog",
                       loghazard = "I", odds = "cloglog", cumhaz = "log")
      } else {
        link <- switch(type, linkS = "I", linkpi = "I", curerate = "cloglog",
                       probcure = "cloglog", survuncured = "log",
                       hazarduncured = "I", cumhazuncured = "I",
                       densityuncured = "I", failuncured = "log2",
                       oddsuncured = "cloglog", loghazarduncured = "I",
                       surv = "log", hazard = "I", density = "I", fail = "log2",
                       loghazard = "I", odds = "cloglog", cumhaz = "I")
      }
    }

    var.link <- switch(link, I = function(x) x, log = function(x) log(x),
                       cloglog = function(x) log(-log(x)), log2 = function(x) -log(x))
    var.link.inv <- switch(link, I = function(x) x, log = function(x) exp(x),
                           cloglog = function(x) exp(-exp(x)), log2 = function(x) exp(-x))

    if (use.gr) {
      if (type == "hazard" && link %in% c("I", "log")) {
        betastar <- beta
        gd <- switch(link, I = t(object$link.surv$gradh(X %*%
                                                          betastar, XD %*% betastar, list(X = X, XD = XD))),
                     log = t(object$link.surv$gradh(X %*% betastar, XD %*%
                                                      betastar, list(X = X, XD = XD))/object$link.surv$h(X %*%
                                                                                                           betastar, XD %*% betastar)))
      }
    }

    lapply(1:length(X), function(i){
      res <- predictnl.default(object, local, var.link = var.link, newdata = newdata[i,, drop = F],
                               type = type, gd = if (use.gr) gd else NULL,
                               X = X[[i]], XD = XD[[i]], X.cr = X.cr[[i]])
      if(var.type == "ci"){
        lower <- var.link.inv(res$Estimate - res$SE * qnorm(0.975))
        upper <- var.link.inv(res$Estimate + res$SE * qnorm(0.975))
        res$lower <- pmin(lower, upper)
        res$upper <- pmax(lower, upper)
        res <- subset(res, select = -SE)
      }

      res$Estimate <- var.link.inv(res$Estimate)
      res
    })
  }
  if (keep.attributes)
    attr(pred, "newdata") <- newdata
  return(pred)
}



predictnl.default <- function (object, fun, newdata = NULL, gd = NULL, ...)
{
  if (is.null(newdata) && !is.null(object$data))
    newdata <- object$data
  localf <- function(coef, ...) {
    object$coefs <- coef[1:length(object$coefs)]
    object$coefs.spline <- coef[(length(object$coefs) + 1):length(coef)]
    fun(object, ...)
  }
  numDeltaMethod(object, localf, newdata = newdata, gd = gd,
                 ...)
}


numDeltaMethod <- function (object, fun, gd = NULL, ...)
{
  coef <- c(object$coefs, object$coefs.spline)
  est <- fun(coef, ...)
  Sigma <- object$covariance
  if (is.null(gd))
    gd <- rstpm2:::grad(fun, coef, ...)
  se.est <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}

local <- function(object, newdata, type = "surv", var.link = function(x) x,
                  X = X, XD = XD, X.cr = X.cr) {
  gamma <- object$coefs
  beta <- object$coefs.spline
  #tt <- object$lm.obj$terms
  link.surv <- object$link.surv
  link.type.cr <- object$link.type.cr
  eta_pi <- as.vector(X.cr %*% gamma)
  eta <- as.vector(X %*% beta)
  etaD <- as.vector(XD %*% beta)

  pi <- get.link(link.type.cr)(eta_pi)
  Su <- link.surv$ilink(eta)
  hazu <- link.surv$h(eta, etaD)
  Hu = link.surv$H(eta)
  S <- if(object$type == "mixture") pi + (1 - pi) * Su else pi ^ (1 - Su)
  dSu <- link.surv$gradS(eta, etaD)
  haz <- - (1 - pi) * dSu / Su
  if (!object$excess && any(h < 0))
    warning(sprintf("Predicted hazards less than zero (n=%i).",
                    sum(h < 0)))
  H <- -log(S)
  Sigma = object$covariance
  est <- switch(type, linkS = eta, linkpi = eta_pi, curerate = pi,
                probcure = pi / S, survuncured = Su,
                hazarduncured = hazu, cumhazuncured = Hu,
                densityuncured = Su * hazu, failuncured = 1 - Su,
                oddsuncured = (1 - Su)/Su, loghazarduncured = log(hazu),
                surv = S, hazard = haz, density = S * haz, fail = 1 - S,
                loghazard = log(haz), odds = (1 - S) / S, cumhaz = H)

  est <- var.link(est)
  return(est)
}





plot.fcm2 <- function(object, newdata = NULL, type = c("surv", "probcure", "survuncured", "hazarduncured",
                                                       "cumhazuncured", "densityuncured", "failuncured",
                                                       "oddsuncured", "loghazarduncured", "hazard",
                                                       "density", "fail", "loghazard", "odds", "cumhaz"),
                      time = NULL, xlim = NULL, ylim = c(0, 1),
                      xlab = "Time", ylab = NULL, col = 1, ci = NULL,
                      add = F, ...){

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

  pred <- predict2.fcm(object, newdata, time, type = type,
                       var.type = ci, indi = TRUE)

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
