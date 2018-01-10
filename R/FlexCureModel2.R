#' Fit spline-based mixture cure model
#'
#' The following function fits a generalized mixture or non-mixture cure model
#' using a link function for the cure rate and for the survival of the uncured, i.e.,
#' \deqn{R(t|z) = \pi(z) + (1 - \pi(z)) S_u(t|z)},
#' where
#' \deqn{g_1[S_u(t|z)] = \eta_1(z) and g_2(\pi(z)) = \eta_2(z)}.
#'
#' @param formula Formula for modelling the the cure rate. Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variables names in \code{formula}, \code{smooth.formula}.
#' @param bhazard Background hazard.
#' @param smooth.formula Formula for modelling the disease-specific survival of the uncured.
#' @param knots Knots used for the baseline hazard in the disease-specific survival function.
#' @param n.knots Number of knots for the disease-specific survival function.
#' The knots are calculated as the equidistant quantiles of the uncensored event-times.
#' If \code{knots} is supplied, this argument will be ignored.
#' @param knots.time A named list containing the knots of each of time-varying covariate effect.
#' @param n.knots.time A named list, containing the number of knots for the time-varying covariate effects.
#' The knots are calculated as the equidistant quantiles of the uncensored event-times.
#' If \code{knots.time} is supplied, this argument will be ignored.
#' @param covariance Logical. If \code{TRUE} (default), the covariance matrix is computed.
#' @param verbose Logical. If \code{TRUE} status messages of the function is outputted.
#' @param type A character indicating the type of cure model.
#' Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param linkpi Character giving the link function selected for the cure rate.
#' Possible values are \code{logit} (default), \code{identity}, \code{loglog}, and \code{probit}.
#' @param linksu Character giving the link function selected for the survival of the uncured.
#' Possible values are \code{loglog} (default), \code{logit}, and \code{probit}.
#' @param constr.optim Logical. If \code{TRUE} the model is fitted using constraints optimization yielding
#' a non-negative hazard of the uncured (default is \code{FALSE}).
#' This option is only implemented for \code{linksu = loglog}.
#' @param ortho Logical. If \code{TRUE} (default), all splines are orthogonalized using a QR-decomposition.
#' @param optim.args List with additional arguments passed to \code{optim}.
#' @param ini.types Character vector denoting the executed schemes for computing initial values.
#' @return An object of class \code{fcm}.
#' @export
#' @import survival
#' @import rstpm2
#' @import numDeriv
#' @example inst/FlexCureModel.ex.R


FlexCureModel2 <- function(formula, data, smooth.formula = NULL, smooth.args = NULL,
                           df = 3, logH.args = NULL, logH.formula = NULL, tvc = NULL,
                           tvc.formula = NULL, bhazard = NULL, cr.formula = ~ 1,
                           type = "mixture",
                           link.type.cr = c("logit", "loglog", "identity"),
                           link.type = c("PH", "PO", "probit", "AH", "AO"),
                           init = NULL, timeVar = "",
                           covariance = T, verbose = T,
                           control = list(maxit = 10000, method = c("Nelder-Mead")),
                           ini.types = c("cure", "flexpara"), cure = FALSE){

  if(!type %in% c("mixture", "nmixture"))
    stop("Wrong specication of argument 'type', must be either 'mixture' or 'nmixture'")

  link.type <- match.arg(link.type)
  link.surv <- switch(link.type, PH = rstpm2:::link.PH, PO = rstpm2:::link.PO, probit = rstpm2:::link.probit,
                      AH = rstpm2:::link.AH)

  link.type.cr <- match.arg(link.type.cr)

  if (!is.null(smooth.formula) && is.null(logH.formula))
    logH.formula <- smooth.formula
  if (!is.null(smooth.args) && is.null(logH.args))
    logH.args <- smooth.args

  eventInstance <- eval(rstpm2:::lhs(formula),envir = data)
  stopifnot(length(rstpm2:::lhs(formula)) >= 2)
  eventExpr <- rstpm2:::lhs(formula)[[length(rstpm2:::lhs(formula))]]
  delayed <- length(rstpm2:::lhs(formula)) >= 4
  surv.type <- attr(eventInstance, "type")
  if (surv.type %in% c("interval2", "left", "mstate"))
    stop("stpm2 not implemented for Surv type ", surv.type,
         ".")

  counting <- attr(eventInstance, "type") == "counting"
  interval <- attr(eventInstance, "type") == "interval"
  timeExpr <- rstpm2:::lhs(formula)[[ifelse(delayed, 3, 2)]]
  if (timeVar == "")
    timeVar <- all.vars(timeExpr)

  if (is.null(logH.formula) && is.null(logH.args)) {
    logH.args$df <- df
    if (cure)
      logH.args$cure <- cure
  }
  if (is.null(logH.formula))
    logH.formula <- as.formula(call("~", as.call(c(quote(nsx),
                                                   call("log", timeExpr), rstpm2:::vector2call(logH.args)))))
  if (is.null(tvc.formula) && !is.null(tvc)) {
    tvc.formulas <- lapply(names(tvc), function(name) call(":",
                                                           as.name(name),
                                                           as.call(c(quote(nsx),
                                                                     call("log",
                                                                          timeExpr),
                                                                     rstpm2:::vector2call(if (cure) list(cure = cure,
                                                                                                         df = tvc[[name]]) else list(df = tvc[[name]]))))))

    if (length(tvc.formulas) > 1)
      tvc.formulas <- list(Reduce(rstpm2:::`%call+%`, tvc.formulas))
    tvc.formula <- as.formula(call("~", tvc.formulas[[1]]))
  }
  if (!is.null(tvc.formula)) {
    rstpm2:::rhs(logH.formula) <- rstpm2:::rhs(logH.formula) %call+% (rstpm2:::rhs(tvc.formula))
  }

  full.formula <- formula
  rstpm2:::rhs(full.formula) <- rstpm2:::rhs(formula) %call+% rstpm2:::rhs(logH.formula)

  .include <- apply(model.matrix(formula, data, na.action = na.pass),
                    1, function(row) !any(is.na(row))) & !is.na(eval(eventExpr,
                                                                     data)) & !is.na(eval(timeExpr, data))
  data <- data[.include, , drop = FALSE]
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "contrasts", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  time <- eval(timeExpr, data, parent.frame())
  if (any(time > 0 & time < 1e-04))
    warning("Some event times < 1e-4: consider transforming time to avoid problems with finite differences")
  time0Expr <- NULL
  if (delayed) {
    time0Expr <- rstpm2:::lhs(formula)[[2]]
    if (time0Var == "")
      time0Var <- all.vars(time0Expr)
    time0 <- eval(time0Expr, data, parent.frame())
    if (any(time0 > 0 & time0 < 1e-04))
      warning("Some entry times < 1e-4: consider transforming time to avoid problems with finite differences")
  }
  event <- eval(eventExpr, data)
  event <- if (length(unique(event)) == 1){
    rep(TRUE, length(event))
  } else {
    event <- event > min(event)
  }

  # if (!interval) {
  #   coxph.call <- mf
  #   coxph.call[[1L]] <- as.name("coxph")
  #   coxph.call$model <- TRUE
  #   coxph.obj <- eval(coxph.call, envir = parent.frame())
  #   y <- model.extract(model.frame(coxph.obj), "response")
  #   data$logHhat <- pmax(-18, link.surv$link(rstpm2:::Shat(coxph.obj)))
  # }
  # if (interval) {
  #   survreg.call <- mf
  #   survreg.call[[1L]] <- as.name("survreg")
  #   survreg.obj <- eval(survreg.call, envir = parent.frame())
  #   weibullShape <- 1/survreg.obj$scale
  #   weibullScale <- predict(survreg.obj)
  #   y <- model.extract(model.frame(survreg.obj), "response")
  #   data$logHhat <- pmax(-18, link$link(pweibull(time, weibullShape,
  #                                                weibullScale, lower.tail = FALSE)))
  # }
  # lm.call <- mf
  # lm.call[[1L]] <- as.name("lm")
  # lm.formula <- full.formula
  # rstpm2:::lhs(lm.formula) <- quote(logHhat)
  # lm.call$formula <- lm.formula
  # dataEvents <- data[event, ]
  # if (interval)
  #   dataEvents <- data
  # lm.call$data <- quote(dataEvents)
  # lm.obj <- eval(lm.call)
  # mt <- terms(lm.obj)
  # mf <- model.frame(lm.obj)


  lm.call <- mf
  lm.call[[1L]] <- as.name("lm")
  lm.formula <- full.formula
  rstpm2:::lhs(lm.formula) <- quote(arbri)
  lm.call$formula <- lm.formula
  dataEvents <- data[event, ]
  dataEvents$arbri <- rnorm(nrow(dataEvents))
  if (interval)
    dataEvents <- data
  lm.call$data <- quote(dataEvents)
  lm.obj <- eval(lm.call)


  #Create background hazard
  if(is.null(bhazard)){
    bhazard <- rep(0, nrow(data))
  }else {
    if(!is.numeric(bhazard)){
      bhazard <- data[, bhazard]
    }
  }
  excess <- !all(bhazard == 0)


  if(length(bhazard) != nrow(data))
    stop("Length of bhazard is not the same as nrow(data)")

  lpfunc <- function(delta, fit, dataset, var) {
    dataset[[var]] <- dataset[[var]] + delta
    rstpm2:::lpmatrix.lm(fit, dataset)
  }

  transX <- function(X, data) X
  transXD <- function(XD) XD


  if (!interval) {
    X <- rstpm2:::lpmatrix.lm(lm.obj, data)
    if (link.type == "AH") {
      datat0 <- data
      datat0[[timeVar]] <- 0
      index0 <- which.dim(X - lpmatrix.lm(lm.obj, datat0))
      transX <- function(X, data) {
        datat0 <- data
        datat0[[timeVar]] <- 0
        Xt0 <- lpmatrix.lm(lm.obj, datat0)
        (X - Xt0)[, index0, drop = FALSE]
      }
      transXD <- function(XD) XD[, index0, drop = FALSE]
      init <- init[index0]
    }
    X <- transX(X, data)
    XD <- rstpm2:::grad(lpfunc, 0, lm.obj, data, timeVar)
    XD <- transXD(matrix(XD, nrow = nrow(X)))
    X1 <- matrix(0, nrow(X), ncol(X))
    X0 <- matrix(0, 1, ncol(X))
    if (delayed && all(time0 == 0))
      delayed <- FALSE
    if (delayed) {
      ind0 <- time0 > 0
      map0 <- vector("integer", nrow(X))
      map0[ind0] <- as.integer(1:sum(ind0))
      map0[!ind0] <- NaN
      which0 <- 1:nrow(X)
      which0[!ind0] <- NaN
      data0 <- data[ind0, , drop = FALSE]
      data0[[timeVar]] <- data0[[time0Var]]
      X0 <- transX(lpmatrix.lm(lm.obj, data0), data0)
      wt0 <- wt[ind0]
      rm(data0)
    }
  }
  else {
    ttype <- eventInstance[, 3]
    X1 <- transX(lpmatrix.lm(lm.obj, data), data)
    data0 <- data
    data0[[timeVar]] <- data0[[time0Var]]
    X <- transX(lpmatrix.lm(lm.obj, data0), data0)
    XD <- grad(lpfunc, 0, lm.obj, data0, timeVar)
    XD <- transXD(matrix(XD, nrow = nrow(X)))
    X0 <- matrix(0, nrow(X), ncol(X))
    rm(data0)
  }

  X.cr <- model.matrix(cr.formula, data = data)

  if(is.null(init)){
    if(verbose) cat("Finding initial values... ")
    init <- vector("list", length(ini.types))
    for(i in 1:length(init)){
      args <- list(formula = formula, data = data, smooth.formula = smooth.formula,
                   logH.formula = logH.formula, tvc.formula = tvc.formula, cr.formula = cr.formula,
                   full.formula = full.formula, X = X, X.cr = X.cr,
                   bhazard = bhazard, type = type, link.type.cr = link.type.cr,
                   link.surv = link.surv, time = time, timeExpr = as.character(timeExpr),
                   lm.obj = lm.obj, method = ini.types[i])
      init[[i]] <- do.call(get.init, args)
    }
  } else {
    if(verbose) cat("Initial values provided by the user... ")
  }


  #Extract minus log likelihood function
  minusloglik <- switch(type,
                        mixture = flexible_mixture_minuslog_likelihood2,
                        nmixture = flexible_nmixture_minuslog_likelihood2)

  #Prepare optimization arguments
  args <- list(event = event, X = X, XD = XD, X.cr = X.cr,
               bhazard = bhazard, link.type.cr = link.type.cr,
               link.surv = link.surv)

  if(is.null(control$maxit)){
    control$maxit <- 10000
  }

  optim.args <- c(control = list(control), args)

  #Test if initial values are within the feasible region
  ini.eval <- sapply(init, function(inival) do.call(minusloglik, c(args, list(inival))))
  run.these <- !is.na(ini.eval)

  if(all(!run.these))
    stop("Initial values are outside feasible region")

  if(verbose) cat("Completed!\nFitting the model... ")

  optim.args$fn <- minusloglik
  res_list <- lapply(init[run.these], function(inival){
    optim.args$par <- inival
    suppressWarnings(do.call(optim, optim.args))
  })

  #Choose the best model according to the maximum likelihood estimate
  MLs <- sapply(res_list, function(x) tail(x$value, 1))
  wh <- which.min(MLs)
  res <- res_list[[wh]]

  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(verbose) cat("Completed!\n")

  #Compute the covariance matrix matrix
  if(covariance){
    args$x0 <- res$par
    args$f <- minusloglik
    hes <- do.call(pracma::hessian, args)
    cov <- solve(hes)
    if(any(is.na(cov))){
      warning("Hessian is not invertible!")
    }
  }else{
    cov <- NULL
  }

  #Output the results
  L <- list(formula = formula, smooth.formula = smooth.formula, tvc.formula = tvc.formula,
            logH.formula = logH.formula, cr.formula = cr.formula,
            coefs = res$par[1:ncol(X.cr)],
            coefs.spline = res$par[(ncol(X.cr) + 1):length(res$par)],
            data = data, NegMaxLik = min(MLs), covariance = cov, ci = covariance,
            type = type, NegMaxLiks = MLs, optim.pars = optim.args[c("control", "fn")],
            args = args, timeExpr = timeExpr, lm.obj = lm.obj, link.type.cr = link.type.cr,
            link.surv = link.surv, excess = excess, timeVar = timeVar, transX = transX, transXD = transXD,
            time = time, event = event, eventExpr = eventExpr)

  class(L) <- c("fcm2", "cuRe")
  L
}


#Function for computing initial values
get.init <- function(formula, data, smooth.formula, logH.formula, tvc.formula, cr.formula, full.formula,
                     bhazard, type, link.type.cr, link.surv, timeExpr, time, lm.obj, X, X.cr, method){

  if(!method %in% c("cure", "flexpara")){
    stop("Argument method should be either 'cure' or 'flexpara'")
  }

  if(method == "cure"){

    formula.pi <- cr.formula
    rstpm2:::lhs(formula.pi) <- rstpm2:::lhs(formula)
    formula.k1 <- formula
    rstpm2:::lhs(formula.k1) <- NULL

    if(length(attr(terms(formula.k1), "term.labels"))){
      a <- Reduce(paste, deparse(formula.k1))
      a <- gsub("-1", "1", a)
      formula.k1 <- as.formula(a)
    } else {
      formula.k1 <- ~ 1
    }

    #Fit mixture or non-mixture cure model
    fit <- fit.cure.model(formula = formula.pi, data = data, bhazard = bhazard, covariance = F,
                          formula.k1 = formula.k1, formula.k2 = ~ 1, type = type)

    #Scale by link function
    pi_hat <- predict(fit, type = "curerate", newdata = data)[,1]

    #Predict survival of the uncured
    lp <- exp(model.matrix(formula.k1, data = data) %*% fit$coefs[[2]])
    suhat <- exp(-lp * time ^ exp(fit$coefs[[3]]))

  }else if(method == "flexpara"){

    formula.2 <- formula
    vars1 <- attr(terms(cr.formula), "term.labels")
    vars2 <- attr(terms(formula.2), "term.labels")

    wh <- which(!vars1 %in% vars2)
    if(length(wh)){
      formula.pi <- as.formula(paste0("~ ", paste(vars1, collapse = "+")))
      rstpm2:::rhs(formula.2) <- rstpm2:::rhs(formula) %call+% rstpm2:::rhs(formula.pi)
    }

    if(length(attr(terms(formula.2), "term.labels"))){
      a <- Reduce(paste, deparse(formula.2))
      a <- gsub("-1", "1", a)
      formula.2 <- as.formula(a)
    } else {
      rstpm2:::rhs(formula.2) <- 1
    }

    #Fit relative survival model
    suppressWarnings(fit <- rstpm2::stpm2(formula.2, data = data,
                                          bhazard = bhazard))


    #Predict survival function
    shat <- predict(fit, newdata = data, se.fit = F, keep.attributes = F)

    #If predictions are all 1, we manually change these
    shat[shat == 1] <- shat[shat == 1] - 0.01

    #Change follow-up times and predict cure rate
    data2 <- data
    data2[, timeExpr] <- max(data2[, timeExpr]) + 0.1
    pi_hat <- predict(fit, newdata = data2, se.fit = F, keep.attributes = F)

    #Change cases with increasing relative survival
    wh <- which(pi_hat >= shat)
    pi_hat[wh] <- shat[wh] - 0.01

    #Run linear model for S_u(t) to obtain initial values for either mixture or non-mixture models
    if(type == "mixture"){
      suhat <- (shat - pi_hat) / (1 - pi_hat)
    } else {
      suhat <- 1 - log(shat) / log(pi_hat)
    }
  }

  #Run linear model for pi to obtain initial values
  gpi_hat <- get.inv.link(link.type.cr)(pi_hat)
  ini_pi <- lm(gpi_hat ~ -1 + X.cr)$coefficients
  names(ini_pi) <- colnames(X.cr)

  #Run linear model for survival to obtain initial values
  gsuhat <- link.surv$link(suhat)
  finites <- is.finite(gsuhat)
  suppressWarnings(ini_surv <- lm(gsuhat[finites] ~ -1 + X[finites,])$coefficients)
  names(ini_surv) <- colnames(X)

  c(ini_pi, ini_surv)
}


#' @export
#Print function for class fcm
print.fcm <- function(fit){
  cat("Call pi:\n")
  print(fit$formula)
  cat("Call S_u(t):\n")
  print(fit$formula_main)
  cat("\nCoefficients:\n")
  print(list(pi = fit$coefs,
             surv = fit$coefs.spline))
}

#' @export
#Summary function for class fcm
summary.fcm <- function(fit){
  se <- sqrt(diag(fit$covariance))
  tval <- c(fit$coefs, fit$coefs.spline) / se
  coefs <- c(fit$coefs, fit$coefs.spline)
  TAB1 <- cbind(Estimate = fit$coefs,
                StdErr = se[1:length(fit$coefs)],
                t.value = tval[1:length(fit$coefs)],
                p.value = ifelse(is.na(tval[1:length(fit$coefs)]), rep(NA, length(fit$coefs)),
                                 2*pt(-abs(tval[1:length(fit$coefs)]), df = fit$df)))

  TAB2 <- cbind(Estimate = fit$coefs.spline,
                StdErr = se[1:length(fit$coefs.spline)],
                t.value = tval[1:length(fit$coefs.spline)],
                p.value = ifelse(is.na(tval[1:length(fit$coefs.spline)]), rep(NA, length(fit$coefs.spline)),
                                 2*pt(-abs(tval[1:length(fit$coefs.spline)]), df = fit$df)))


  results <- list(pi = TAB1, surv = TAB2)
  results$type <- fit$type
  results$linkpi <- fit$linkpi
  results$linksu <- fit$linksu
  results$ML <- fit$NegMaxLik
  results$formula <- fit$formula
  results$formula.fix <- fit$formula_main
  results$formula.tvc <- fit$tvc.formula
  class(results) <- "summary.fcm"
  results
}

#' @export
#Print for class summary.fcm
print.summary.fcm <- function(x)
{
  cat("Call - pi:\n")
  print(x$formula)
  #    cat("\n")
  printCoefmat(x$pi, P.value = TRUE, has.Pvalue = T)
  cat("\nCall - surv - baseline: ")
  print(as.formula(deparse(x$formula.fix)))
  if(length(all.vars(x$formula.tvc))){
    cat("Call - surv - tvc: ")
    print(deparse(x$formula.tvc))
  }
  printCoefmat(x$surv, P.value = TRUE, has.Pvalue = T)
  cat("\n")
  cat("Type =", x$type, "\n")
  cat("Link - pi =", x$linkpi, "\n")
  cat("Link - surv = ", x$linksu, "\n")
  cat("LogLik(model) =", x$ML, "\n")
}

