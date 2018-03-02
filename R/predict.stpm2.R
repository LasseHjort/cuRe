predict.stpm2 <- function (object, newdata = NULL, type = c("surv", "cumhaz",
                                                            "hazard", "density", "hr", "sdiff", "hdiff", "loghazard",
                                                            "link", "meansurv", "meansurvdiff", "odds", "or", "margsurv",
                                                            "marghaz", "marghr", "meanhaz", "af", "fail", "margfail",
                                                            "meanmargsurv", "uncured", "probcure"), grid = FALSE, seqLength = 300,
                           se.fit = FALSE, link = NULL, exposed = NULL, var = NULL,
                           keep.attributes = TRUE, use.gr = TRUE, ...)
{
  type <- match.arg(type)
  args <- object@args
  if (type %in% c("fail", "margfail")) {
    out <- 1 - predict.stpm2.base(object, newdata = newdata,
                                  type = switch(type, fail = "surv", margfail = "margsurv"),
                                  grid, seqLength, se.fit, link, exposed, var, keep.attributes,
                                  use.gr, ...)
    if (se.fit) {
      temp <- out$lower
      out$lower <- out$upper
      out$upper <- temp
    }
    return(out)
  }
  if (is.null(exposed) && is.null(var) & type %in% c("hr",
                                                     "sdiff", "hdiff", "meansurvdiff", "or", "marghr", "af",
                                                     "uncured", "probcure"))
    stop("Either exposed or var required for type in (\"hr\",\"sdiff\", \"hdiff\",\"meansurvdiff\",\"or\",\"marghr\",\"af\",\"uncured\", \"probcure\")")
  if (type %in% c("margsurv", "marghaz", "marghr", "margfail",
                  "meanmargsurv") && !object@args$frailty)
    stop("Marginal prediction only for frailty models")
  if (is.null(newdata) && type %in% c("hr", "sdiff", "hdiff",
                                      "meansurvdiff", "or", "marghr", "uncured", "probcure"))
    stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','or','marghr','uncured', 'probcure') requires newdata to be specified.")
  calcX <- !is.null(newdata)
  time <- NULL
  if (is.null(newdata)) {
    X <- object@x
    XD <- object@xd
    y <- object@y
    time <- as.vector(y[, ncol(y) - 1])
    newdata <- as.data.frame(object@data)
  }
  lpfunc <- if (inherits(object, "pstpm2"))
    function(x, ...) {
      newdata2 <- newdata
      newdata2[[object@timeVar]] <- x
      predict(object@gam, newdata2, type = "lpmatrix")
    }
  else function(delta, fit, data, var) {
    data[[var]] <- data[[var]] + delta
    rstpm2:::lpmatrix.lm(fit, data)
  }
  if (grid) {
    Y <- object@y
    event <- Y[, ncol(Y)] == 1 | object@args$interval
    time <- object@data[[object@timeVar]]
    eventTimes <- time[event]
    tt <- seq(min(eventTimes), max(eventTimes), length = seqLength)[-1]
    data.x <- data.frame(tt)
    names(data.x) <- object@timeVar
    newdata[[object@timeVar]] <- NULL
    newdata <- merge(newdata, data.x)
    calcX <- TRUE
  }
  if (calcX) {
    if (inherits(object, "stpm2")) {
      X <- object@args$transX(rstpm2:::lpmatrix.lm(object@lm, newdata),
                              newdata)
      XD <- rstpm2:::grad(lpfunc, 0, object@lm, newdata, object@timeVar)
      XD <- object@args$transXD(matrix(XD, nrow = nrow(X)))
    }
    if (inherits(object, "pstpm2")) {
      X <- object@args$transX(predict(object@gam, newdata,
                                      type = "lpmatrix"), newdata)
      XD <- object@args$transXD(rstpm2:::grad1(lpfunc, newdata[[object@timeVar]]))
    }
  }
  if (is.null(time)) {
    time <- eval(object@timeExpr, newdata, parent.frame())
  }
  if (type %in% c("hr", "sdiff", "hdiff", "meansurvdiff", "or",
                  "marghr", "af", "uncured", "probcure")) {
    newdata2 <- exposed(newdata)
    if (inherits(object, "stpm2")) {
      X2 <- object@args$transX(rstpm2:::lpmatrix.lm(object@lm, newdata2),
                               newdata2)
      XD2 <- rstpm2:::grad(lpfunc, 0, object@lm, newdata2, object@timeVar)
      XD2 <- object@args$transXD(matrix(XD2, nrow = nrow(X)))
    }
    if (inherits(object, "pstpm2")) {
      X2 <- object@args$transX(predict(object@gam, newdata2,
                                       type = "lpmatrix"), newdata2)
      XD2 <- object@args$transXD(rstpm2:::grad1(lpfunc, newdata2[[object@timeVar]]))
    }
  }
  colMeans <- function(x) colSums(x)/apply(x, 2, length)
  if (object@frailty && type %in% c("af", "meansurvdiff") &&
      args$RandDist == "Gamma" && !object@args$interval &&
      !object@args$delayed) {
    times <- newdata[[object@timeVar]]
    utimes <- sort(unique(times))
    n <- nrow(X)/length(utimes)
    n.cluster <- length(unique(args$cluster))
    link <- object@link
    beta <- coef(object)
    npar <- length(beta)
    logtheta <- beta[npar]
    theta <- exp(beta[npar])
    beta <- beta[-npar]
    Hessian <- solve(vcov(object))
    eta <- as.vector(X %*% beta)
    eta2 <- as.vector(X2 %*% beta)
    S <- link$ilink(eta)
    S2 <- link$ilink(eta2)
    H <- -log(S)
    H2 <- -log(S2)
    marg <- function(logtheta, H) (1 + exp(logtheta) * H)^(-1/exp(logtheta))
    margS <- marg(logtheta, H)
    margS2 <- marg(logtheta, H2)
    dmarg.dlogtheta <- function(logtheta, H) {
      theta <- exp(logtheta)
      marg(logtheta, H) * (exp(-logtheta) * log(1 + theta *
                                                  H) - H/(1 + theta * H))
    }
    meanS <- tapply(margS, times, mean)
    meanS2 <- tapply(margS2, times, mean)
    fit <- switch(type, af = 1 - (1 - meanS2)/(1 - meanS),
                  meansurvdiff = meanS - meanS2)
    se.fit <- vector("numeric", length(utimes))
    for (i in 1:length(utimes)) {
      index <- which(times == utimes[i])
      newobj <- object
      newobj@args$X <- X[index, , drop = FALSE]
      newobj@args$XD <- XD[index, , drop = FALSE]
      gradli <- residuals(newobj, type = "gradli")
      res <- cbind(margS[index] - mean(margS[index]), margS2[index] -
                     mean(margS2[index]))
      res <- apply(res, 2, function(col) tapply(col, args$cluster,
                                                sum))
      res <- cbind(res, gradli)
      meat <- stats::var(res, na.rm = TRUE)
      colnames(meat) <- rownames(meat) <- c("S", "S0",
                                            names(beta), "logtheta")
      S.hessian <- cbind(-diag(2) * n/n.cluster, rbind(colSums(margS[index] *
                                                                 (-link$gradH(eta[index], list(X = X[index, ,
                                                                                                     drop = FALSE]))/(1 + theta * H[index])))/n.cluster,
                                                       colSums(margS2[index] * (-link$gradH(eta2[index],
                                                                                            list(X = X2[index, , drop = FALSE]))/(1 + theta *
                                                                                                                                    H2[index])))/n.cluster), c(sum(dmarg.dlogtheta(logtheta,
                                                                                                                                                                                   H[index]))/n.cluster, sum(dmarg.dlogtheta(logtheta,
                                                                                                                                                                                                                             H2[index]))/n.cluster))
      par.hessian <- cbind(matrix(0, nrow = npar, ncol = 2),
                           -Hessian/n.cluster)
      bread <- rbind(S.hessian, par.hessian)
      ibread <- solve(bread)
      sandwich <- (ibread %*% meat %*% t(ibread)/n.cluster)[1:2,
                                                            1:2]
      gradient <- switch(type, af = as.matrix(c(-(1 - meanS2[i])/(1 -
                                                                    meanS[i])^2, 1/(1 - meanS[i])), nrow = 2, ncol = 1),
                         meansurvdiff = matrix(c(1, -1), nrow = 2))
      AF.var <- t(gradient) %*% sandwich %*% gradient
      se.fit[i] <- sqrt(AF.var)
    }
    pred <- data.frame(Estimate = fit, lower = fit - 1.96 *
                         se.fit, upper = fit + 1.96 * se.fit)
    if (keep.attributes)
      attr(pred, "newdata") <- newdata
    return(pred)
  }
  if (object@frailty && type %in% c("meanmargsurv") && args$RandDist ==
      "Gamma" && !object@args$interval && !object@args$delayed) {
    times <- newdata[[object@timeVar]]
    utimes <- sort(unique(times))
    n <- nrow(X)/length(utimes)
    n.cluster <- length(unique(args$cluster))
    link <- object@link
    beta <- coef(object)
    npar <- length(beta)
    logtheta <- beta[npar]
    theta <- exp(beta[npar])
    beta <- beta[-npar]
    Hessian <- solve(vcov(object))
    eta <- as.vector(X %*% beta)
    S <- link$ilink(eta)
    H <- -log(S)
    marg <- function(logtheta, H) (1 + exp(logtheta) * H)^(-1/exp(logtheta))
    margS <- marg(logtheta, H)
    dmarg.dlogtheta <- function(logtheta, H) {
      theta <- exp(logtheta)
      marg(logtheta, H) * (exp(-logtheta) * log(1 + theta *
                                                  H) - H/(1 + theta * H))
    }
    meanS <- tapply(margS, times, mean)
    fit <- meanS
    se.fit <- vector("numeric", length(utimes))
    for (i in 1:length(utimes)) {
      index <- which(times == utimes[i])
      newobj <- object
      newobj@args$X <- X[index, , drop = FALSE]
      newobj@args$XD <- XD[index, , drop = FALSE]
      gradli <- residuals(newobj, type = "gradli")
      res <- tapply(margS[index] - mean(margS[index]),
                    args$cluster, sum)
      res <- cbind(res, gradli)
      meat <- stats::var(res, na.rm = TRUE)
      colnames(meat) <- rownames(meat) <- c("S", names(beta),
                                            "logtheta")
      S.hessian <- c(-n/n.cluster, colSums(margS[index] *
                                             (-link$gradH(eta[index], list(X = X[index, ,
                                                                                 drop = FALSE]))/(1 + theta * H[index])))/n.cluster,
                     sum(dmarg.dlogtheta(logtheta, H[index]))/n.cluster)
      par.hessian <- cbind(matrix(0, nrow = npar, ncol = 1),
                           -Hessian/n.cluster)
      bread <- rbind(S.hessian, par.hessian)
      ibread <- solve(bread)
      sandwich <- (ibread %*% meat %*% t(ibread)/n.cluster)[1,
                                                            1]
      se.fit[i] <- sqrt(sandwich)
    }
    pred <- data.frame(Estimate = fit, lower = fit - 1.96 *
                         se.fit, upper = fit + 1.96 * se.fit)
    if (keep.attributes)
      attr(pred, "newdata") <- newdata
    return(pred)
  }
  local <- function(object, newdata = NULL, type = "surv",
                    exposed) {
    beta <- coef(object)
    tt <- object@terms
    link <- object@link
    if (object@frailty) {
      theta <- exp(beta[length(beta)])
      beta <- beta[-length(beta)]
      if (object@args$RandDist == "LogN") {
        gauss_x <- object@args$gauss_x
        gauss_w <- object@args$gauss_w
        Z <- model.matrix(args$Z.formula, newdata)
        if (ncol(Z) > 1)
          stop("Current implementation only allows for a single random effect")
        Z <- as.vector(Z)
      }
    }
    eta <- as.vector(X %*% beta)
    etaD <- as.vector(XD %*% beta)
    S <- link$ilink(eta)
    h <- link$h(eta, etaD)
    if (!object@args$excess && any(h < 0))
      warning(sprintf("Predicted hazards less than zero (n=%i).",
                      sum(h < 0)))
    H = link$H(eta)
    Sigma = vcov(object)
    if (type == "link") {
      return(eta)
    }
    if (type == "cumhaz") {
      return(H)
    }
    if (type == "density")
      return(S * h)
    if (type == "surv") {
      return(S)
    }
    if (type == "fail") {
      return(1 - S)
    }
    if (type == "odds") {
      return((1 - S)/S)
    }
    if (type == "sdiff")
      return(link$ilink(as.vector(X2 %*% beta)) - S)
    if (type == "hazard") {
      return(h)
    }
    if (type == "loghazard") {
      return(log(h))
    }
    if (type == "hdiff") {
      eta2 <- as.vector(X2 %*% beta)
      etaD2 <- as.vector(XD2 %*% beta)
      h2 <- link$h(eta2, etaD2)
      return(h2 - h)
    }
    if (type == "uncured") {
      S2 <- link$ilink(as.vector(X2 %*% beta))
      return((S - S2)/(1 - S2))
    }
    if(type == "probcure"){
      S2 <- link$ilink(as.vector(X2 %*% beta))
      return(S2 / S)
    }
    if (type == "hr") {
      eta2 <- as.vector(X2 %*% beta)
      etaD2 <- as.vector(XD2 %*% beta)
      h2 <- link$h(eta2, etaD2)
      return(h2/h)
    }
    if (type == "or") {
      S2 <- link$ilink(as.vector(X2 %*% beta))
      return((1 - S2)/S2/((1 - S)/S))
    }
    if (type == "meansurv") {
      return(tapply(S, newdata[[object@timeVar]], mean))
    }
    if (type == "meanhaz") {
      return(tapply(S * h, newdata[[object@timeVar]], sum)/tapply(S,
                                                                  newdata[[object@timeVar]], sum))
    }
    if (type == "meansurvdiff") {
      eta2 <- as.vector(X2 %*% beta)
      S2 <- link$ilink(eta2)
      return(tapply(S2, newdata[[object@timeVar]], mean) -
               tapply(S, newdata[[object@timeVar]], mean))
    }
    if (type == "af") {
      eta2 <- as.vector(X2 %*% beta)
      S2 <- link$ilink(eta2)
      meanS <- tapply(S, newdata[[object@timeVar]], mean)
      meanS2 <- tapply(S2, newdata[[object@timeVar]], mean)
      if (object@frailty) {
        if (object@args$RandDist == "Gamma") {
          meanS <- tapply((1 + theta * (-log(S)))^(-1/theta),
                          newdata[[object@timeVar]], mean)
          meanS2 <- tapply((1 + theta * (-log(S2)))^(-1/theta),
                           newdata[[object@timeVar]], mean)
        }
        else {
          meanS <- tapply(sapply(1:length(gauss_x), function(i) link$ilink(eta +
                                                                             Z * sqrt(2) * sqrt(theta) * gauss_x[i])) %*%
                            gauss_w/sqrt(pi), newdata[[object@timeVar]],
                          mean)
          meanS2 <- tapply(sapply(1:length(gauss_x),
                                  function(i) link$ilink(eta2 + Z * sqrt(2) *
                                                           sqrt(theta) * gauss_x[i])) %*% gauss_w/sqrt(pi),
                           newdata[[object@timeVar]], mean)
        }
      }
      return((meanS2 - meanS)/(1 - meanS))
    }
    if (type == "meanmargsurv") {
      stopifnot(object@frailty && object@args$RandDist %in%
                  c("Gamma", "LogN"))
      if (object@args$RandDist == "Gamma")
        return(tapply((1 + theta * H)^(-1/theta), newdata[[object@timeVar]],
                      mean))
      if (object@args$RandDist == "LogN") {
        return(tapply(sapply(1:length(gauss_x), function(i) link$ilink(eta +
                                                                         Z * sqrt(2) * sqrt(theta) * gauss_x[i])) %*%
                        gauss_w/sqrt(pi), newdata[[object@timeVar]],
                      mean))
      }
    }
    if (type == "margsurv") {
      stopifnot(object@args$frailty && object@args$RandDist %in%
                  c("Gamma", "LogN"))
      if (object@args$RandDist == "Gamma")
        return((1 + theta * H)^(-1/theta))
      if (object@args$RandDist == "LogN") {
        return(sapply(1:length(gauss_x), function(i) link$ilink(eta +
                                                                  Z * sqrt(2) * sqrt(theta) * gauss_x[i])) %*%
                 gauss_w/sqrt(pi))
      }
    }
    if (type == "marghaz") {
      stopifnot(object@frailty && object@args$RandDist %in%
                  c("Gamma", "LogN"))
      if (object@args$RandDist == "Gamma") {
        return(h/(1 + H * theta))
      }
      if (object@args$RandDist == "LogN") {
        return(sapply(1:length(gauss_x), function(i) link$h(eta +
                                                              Z * sqrt(2) * sqrt(theta) * gauss_x[i], etaD)) %*%
                 gauss_w/sqrt(pi))
      }
    }
    if (type == "marghr") {
      stopifnot(object@frailty && object@args$RandDist %in%
                  c("Gamma", "LogN"))
      eta2 <- as.vector(X2 %*% beta)
      etaD2 <- as.vector(XD2 %*% beta)
      if (object@args$RandDist == "Gamma") {
        H2 <- link$H(eta2)
        h2 <- link$h(eta2, etaD2)
        margsurv <- (1 + theta * H)^(-1/theta)
        marghaz <- h * margsurv^theta
        margsurv2 <- (1 + theta * H2)^(-1/theta)
        marghaz2 <- h2 * margsurv2^theta
      }
      if (object@args$RandDist == "LogN") {
        marghaz <- sapply(1:length(gauss_x), function(i) as.vector(link$h(eta +
                                                                            Z * sqrt(2) * sqrt(theta) * gauss_x[i], etaD))) %*%
          gauss_w/sqrt(pi)
        marghaz2 <- sapply(1:length(gauss_x), function(i) as.vector(link$h(eta2 +
                                                                             Z * sqrt(2) * sqrt(theta) * gauss_x[i], etaD2))) %*%
          gauss_w/sqrt(pi)
      }
      return(marghaz2/marghaz)
    }
  }
  pred <- if (!se.fit) {
    local(object, newdata, type = type, exposed = exposed,
          ...)
  }
  else {
    gd <- NULL
    beta <- coef(object)
    if (is.null(link))
      link <- switch(type, surv = "cloglog", cumhaz = "log",
                     hazard = "log", hr = "log", sdiff = "I", hdiff = "I",
                     loghazard = "I", link = "I", odds = "log", or = "log",
                     margsurv = "cloglog", marghaz = "log", marghr = "log",
                     meansurv = "I", meanhaz = "I", af = "I",
                     uncured = "cloglog", probcure = "cloglog")
    if (use.gr) {
      colMeans <- function(x) apply(x, 2, mean)
      collapse <- function(gd) do.call("cbind", tapply(1:nrow(gd),
                                                       newdata[[object@timeVar]], function(index) colMeans(gd[index,
                                                                                                              , drop = FALSE])))
      collapse1 <- function(S) as.vector(tapply(S, newdata[[object@timeVar]],
                                                mean))
      fd <- function(f, x, eps = 1e-05) t(sapply(1:length(x),
                                                 function(i) {
                                                   upper <- lower <- x
                                                   upper[i] = x[i] + eps
                                                   lower[i] = x[i] - eps
                                                   (f(upper) - f(lower))/2/eps
                                                 }))
      if (type == "hazard" && link %in% c("I", "log")) {
        betastar <- if (args$frailty)
          beta[-length(beta)]
        else beta
        gd <- switch(link, I = t(object@link$gradh(X %*%
                                                     betastar, XD %*% betastar, list(X = X, XD = XD))),
                     log = t(object@link$gradh(X %*% betastar, XD %*%
                                                 betastar, list(X = X, XD = XD))/object@link$h(X %*%
                                                                                                 betastar, XD %*% betastar)))
      }
      if (type == "meansurv" && !object@frailty) {
        gd <- collapse(object@link$gradS(X %*% beta,
                                         X))
      }
      if (type == "meansurvdiff" && !object@frailty) {
        gd <- collapse(object@link$gradS(X2 %*% beta,
                                         X2) - object@link$gradS(X %*% beta, X))
      }
      if (type == "margsurv" && link %in% c("I", "cloglog") &&
          args$RandDist == "Gamma") {
        theta <- exp(beta[length(beta)])
        betastar <- beta[-length(beta)]
        eta <- as.vector(X %*% betastar)
        H <- as.vector(object@link$H(eta))
        gradH <- object@link$gradH(eta, list(X = X))
        S0 <- 1 + theta * H
        margS <- S0^(-1/theta)
        if (link == "I")
          gd <- t(cbind(-margS * gradH/(1 + theta * H),
                        margS * (1/theta * log(1 + theta * H) - H/(1 +
                                                                     theta * H))))
        if (link == "cloglog")
          gd <- t(cbind(-(theta^2 * S0^(-theta - 1) *
                            gradH/(S0^(-theta) * (-theta) * log(S0))),
                        (theta * log(S0) * S0^(-theta) - theta^2 *
                           H * S0^(-1 - theta))/(S0^(-theta) * (-theta *
                                                                  log(S0)))))
      }
      if (type == "marghaz" && link %in% c("I", "log") &&
          args$RandDist == "Gamma") {
        theta <- exp(beta[length(beta)])
        betastar <- beta[-length(beta)]
        eta <- as.vector(X %*% betastar)
        etaD <- as.vector(XD %*% betastar)
        H <- as.vector(object@link$H(eta))
        h <- as.vector(object@link$h(eta, etaD))
        gradH <- object@link$gradH(eta, list(X = X))
        gradh <- object@link$gradh(eta, etaD, list(X = X,
                                                   XD = XD))
        S0 <- 1 + theta * H
        margS <- S0^(-1/theta)
        if (link == "I")
          gd <- t(cbind((S0 * gradh - theta * h * gradH)/(theta^2 *
                                                            H^2 + S0), -(theta * H * h/theta^2 * H^2 +
                                                                           S0)))
        if (link == "log")
          gd <- t(cbind((S0 * gradh - theta * h * gradH)/(S0 *
                                                            h), -theta * H/S0))
      }
      if (type == "af" && !object@frailty) {
        meanS <- collapse1(as.vector(object@link$ilink(X %*%
                                                         beta)))
        meanS2 <- collapse1(as.vector(object@link$ilink(X2 %*%
                                                          beta)))
        gradS <- collapse(object@link$gradS(X %*% beta,
                                            X))
        gradS2 <- collapse(object@link$gradS(X2 %*% beta,
                                             X2))
        gd <- t((t(gradS2 - gradS) * (1 - meanS) + (meanS2 -
                                                      meanS) * t(gradS))/(1 - meanS)^2)
      }
    }
    predictnl.rstpm2(object, local, newdata = newdata,
                     type = type, gd = if (use.gr)
                       gd
                     else NULL, exposed = exposed, ...)
  }
  if (keep.attributes)
    attr(pred, "newdata") <- newdata
  return(pred)
}


predictnl.rstpm2 <- function (object, fun, newdata = NULL, gd = NULL, ...)
{
  if (is.null(newdata) && !is.null(object$data))
    newdata <- object$data
  localf <- function(coefs, ...) {
    if ("coefficients" %in% names(object)) {
      object$coefficients <- coefs
    }
    else if ("coef" %in% names(object)) {
      object$coef <- coefs
    }
    else object@fullcoef <- coefs
    fun(object, ...)
  }
  rstpm2:::numDeltaMethod(object, localf, newdata = newdata, gd = gd,
                          ...)
}
