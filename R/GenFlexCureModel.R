#' Fit generalized mixture cure model
#'
#' The following function fits a generalized mixture or non-mixture cure model
#' using a link function for the cure rate and for the survival of the uncured.
#' For a mixture cure model, the model is specified by
#' \deqn{S(t|z) = \pi(z) + [1 - \pi(z)] S_u(t|z),}
#' where
#' \deqn{g_1[S_u(t|z)] = \eta_1(t, z)\qquad and \qquad g_2[\pi(z)] = \eta_2(z).}
#' The function implements multiple link functions for both \eqn{g_1} and \eqn{g_2}. The default time-effect
#' is natural cubic splines, but the function allows for the use of other smoothers.
#'
#' @param formula Formula for modelling the survival of the uncured. Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variables names in \code{formula},
#' \code{smooth.formula}, \code{tvc.formula}, and \code{cr.formula}.
#' @param smooth.formula Formula for describing the time-effect of the survival of the uncured.
#' If \code{NULL} (default), the function uses the natural cubic splines modelled on the log-time scale.
#' @param smooth.args List. Optional arguments to the time-effect of the survival
#' of the uncured (default is \code{NULL}).
#' @param df Integer. Degrees of freedom (default is 3) for the time-effect of the survival of the uncured.
#' Is not used if \code{smooth.formula} is provided.
# @param logH.args
# @param logH.formula blabal
#' @param tvc Named list of integers. Specifies the degrees of freedom for time-varying covariate effects.
#' For instance, \code{tvc = list(x = 3)} creates a time-varying spline-effect of the covariate "x" with
#' 3 degrees of freedom using the \code{rstpm2::nsx} function.
#' @param tvc.formula Formula for the time-varying covariate effects.
#' For time-varying effects, a linear term of the covariate has to be included in \code{formula}.
#' @param cr.formula Formula for the cure proportion.
#' The left hand side of the formula is not used and should therefore not be specified.
#' @param bhazard Background hazard.
#' @param type A character indicating the type of cure model.
#' Possible values are \code{mixture} for mixture cure models (default) and \code{nmixture}
#' for non-mixture cure models.
#' @param covariance Logical. If \code{TRUE} (default), the covariance matrix is computed.
#' @param verbose Logical. If \code{TRUE} status messages of the function is outputted.
#' @param link.type.cr Character providing the link function for the cure proportion.
#' Possible values are \code{logit} (default), \code{loglog}, \code{identity}, and \code{probit}.
#' @param link.type Character providing the link function for the survival of the uncured.
#' Possible values are \code{PH} for a proportional hazards model (default), \code{PO} for a proportion odds model,
#' and \code{probit} for a probit model.
#' @param init Initial values for the optimization procedure.
#' If not specified, the function will create initial values internally.
#' @param timeVar Optional character giving the name of the variable specifying the time component of the \code{Surv} object.
#' Should currently not be used.
#' @param time0Var Optional character giving the name of the variable specifying the time start time component used for delayed entry.
#' Should currently not be used.
#' @param baseoff Logical. If \code{TRUE}, the time-effect is modelled only using \code{tvc.formula} rather
#' than merging with \code{smooth.formula}.
#' @param control Named list with control arguments passed to \code{optim}.
#' @param method Character passed to \code{optim} indicating the method for optimization.
#' See \code{?optim} for details.
#' @param constraint Logical. Indicates whether non-negativity constraints should be forced upon
#' the hazard of the uncured patients (see details).
#' @param ini.types Character vector denoting the executed schemes for computing initial values (see details).
# @param cure Logical. Indicates whether a cure model specification is needed for the survival of the uncured.
# This is usually \code{FALSE} (default).
#' @return An object of class \code{gfcm}.
#' @details The default smoother is natural cubic splines established by the \code{rstpm2::nsx} function.
#' Functions such as \code{ns}, \code{bs} are readily available for usage. Also the \code{cb} function in this package
#' can be used. Initial values are calculated by two procedures and the model is fitted under each set of initial values.
#' The model producing the highest likelihood is selected.\cr
#'
#' Using \code{link.type = 'PH'}, the link function \eqn{g_1(x) = \log(-\log(x))} is used.
#' Using \code{link.type = 'PO'}, the link function \eqn{g_1(x) = \log(\frac{x}{1 - x})} is used.\cr
#'
#' If \code{constraint = TRUE}, a non-negative hazard of the uncured is ensured by a general penalization scheme.
#' If \code{constraint = FALSE}, penalization is still employed, but on the all-cause hazard instead.
#' @export
#' @import survival
#' @import rstpm2
#' @import numDeriv
#' @import date
#' @example inst/GenFlexCureModel.ex.R

#Note that the order of the terms in the smooth.formula matters.
#Meaning: smooth.formula = nsx(df = 4, x = something), will produce an error, while nsx(x = smothing, df= 4) will not
GenFlexCureModel <- function(formula, data, smooth.formula = NULL, smooth.args = NULL,
                             df = 3, tvc = NULL,
                             tvc.formula = NULL, bhazard = NULL, cr.formula = ~ 1,
                             type = "mixture",
                             link.type.cr = c("logit", "loglog", "identity", "probit"),
                             link.type = c("PH", "PO", "probit"),
                             init = NULL, baseoff = FALSE, timeVar = "", time0Var = "",
                             covariance = T, verbose = T,
                             control = list(maxit = 10000), method = "Nelder-Mead",
                             constraint = TRUE,
                             ini.types = c("cure", "flexpara")){

  if(!type %in% c("mixture", "nmixture"))
    stop("Wrong specication of argument 'type', must be either 'mixture' or 'nmixture'")

  logH.args <- NULL
  logH.formula <- NULL
  cure <- FALSE

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

  ## set up the formulae
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

  if (baseoff)
    rstpm2:::rhs(logH.formula) <- rstpm2:::rhs(tvc.formula)

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

  ## initialise values specific to either delayed entry or interval-censored
  ind0 <- FALSE
  map0 <- 0L
  which0 <- 0
  #wt0 <- 0
  ttype <- 0
  transX <- function(X, data) X
  transXD <- function(XD) XD

  if (!interval) {
    X <- rstpm2:::lpmatrix.lm(lm.obj, data)
    if (link.type == "AH") {
      datat0 <- data
      datat0[[timeVar]] <- 0
      index0 <- which.dim(X - rstpm2:::lpmatrix.lm(lm.obj, datat0))
      transX <- function(X, data) {
        datat0 <- data
        datat0[[timeVar]] <- 0
        Xt0 <- rstpm2:::lpmatrix.lm(lm.obj, datat0)
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
      X0 <- transX(rstpm2:::lpmatrix.lm(lm.obj, data0), data0)
      #wt0 <- wt[ind0]
      rm(data0)
    }
  } else {
    ttype <- eventInstance[, 3]
    X1 <- transX(rstpm2:::lpmatrix.lm(lm.obj, data), data)
    data0 <- data
    data0[[timeVar]] <- data0[[time0Var]]
    X <- transX(rstpm2:::lpmatrix.lm(lm.obj, data0), data0)
    XD <- rstpm2:::grad(lpfunc, 0, lm.obj, data0, timeVar)
    XD <- transXD(matrix(XD, nrow = nrow(X)))
    X0 <- matrix(0, nrow(X), ncol(X))
    rm(data0)
  }

  #Create linear object for the cure rate part
  lm.call <- mf
  lm.call[[1L]] <- as.name("lm")
  lm.formula <- cr.formula
  rstpm2:::lhs(lm.formula) <- quote(arbri)
  lm.call$formula <- lm.formula
  dataEvents <- data
  dataEvents$arbri <- rnorm(nrow(dataEvents))
  if (interval)
    dataEvents <- data
  lm.call$data <- quote(dataEvents)
  lm.obj.cr <- eval(lm.call)

  X.cr <- rstpm2:::lpmatrix.lm(lm.obj.cr, data)
  #X.cr <- model.matrix(cr.formula, data = data)

  if(is.null(init)){
    if(verbose) cat("Finding initial values... ")
    #ini.types <- if(delayed) ini.types[2] else ini.types
    init <- vector("list", length(ini.types))
    for(i in 1:length(init)){
      args <- list(formula = formula, data = data, smooth.formula = smooth.formula,
                   logH.formula = logH.formula, tvc.formula = tvc.formula, cr.formula = cr.formula,
                   full.formula = full.formula, X = X, X0 = X0, X.cr = X.cr, delayed = delayed,
                   bhazard = bhazard, type = type, link.type.cr = link.type.cr,
                   link.surv = link.surv, time = time, timeExpr = as.character(timeExpr),
                   lm.obj = lm.obj, method = ini.types[i])
      init[[i]] <- do.call(get.init, args)
    }
  } else {
    if(verbose) cat("Initial values provided by the user... ")
  }


  # #Extract minus log likelihood function
  # if(delayed) {
  #   minusloglik <- switch(type,
  #                         mixture = GenFlexMixMinLogLikDelayed,
  #                         nmixture = GenFlexNmixMinLogLikDelayed)
  # } else {
  #   minusloglik <- switch(type,
  #                         mixture = GenFlexMixMinLogLik,
  #                         nmixture = GenFlexNmixMinLogLik)
  # }


  cure.type <- switch(type,
                      mixture = mix,
                      nmixture = nmix)

  #minusloglik <- switch(type,
  #                      mixture = GenFlexMixMinLogLikDelayed,
  #                      nmixture = GenFlexNmixMinLogLikDelayed)
  minusloglik <- GenFlexMinLogLikDelayed

  #Prepare optimization arguments
  args <- list(event = event, X = X, XD = XD, X.cr = X.cr, X0 = X0, ind0 = ind0,
               bhazard = bhazard, link.type.cr = link.type.cr,
               link.surv = link.surv, kappa = 0, constraint = FALSE, cure.type = cure.type)

  if(is.null(control$maxit)){
    control$maxit <- 10000
  }


  #Test if initial values are within the feasible region
  ini.eval <- sapply(init, function(inival) do.call(minusloglik, c(args, list(inival))))
  run.these <- !is.na(ini.eval)

  if(all(!run.these))
    stop("Initial values are outside feasible region")

  if(verbose) cat("Completed!\nFitting the model... ")

  optim.args <- c(control = list(control), args)
  optim.args$kappa <- 1
  optim.args$fn <- minusloglik
  optim.args$method <- method
  optim.args$constraint <- constraint
  optim.args$cure.type <- cure.type
  res_list <- vector("list", length(init[run.these]))
  for(i in 1:length(res_list)){
    neghaz <- T
    while(neghaz){
      #cat(optim.args$kappa, "\n")
      optim.args$par <- init[run.these][[i]]
      res.optim <- do.call(optim, optim.args)
      beta <- res.optim$par[(ncol(X.cr) + 1):length(res.optim$par)]
      eta <- X %*% beta
      etaD <- XD %*% beta
      if(constraint){
        haz.const <- link.surv$h(eta, etaD)
      } else {
        gamma <- res.optim$par[1:ncol(X.cr)]
        eta.pi <- X.cr %*% gamma
        pi <- get.link(link.type.cr)(eta.pi)
        surv <- link.surv$ilink(eta)
        rsurv <- cure.type$surv(pi, surv)
        ehaz <- cure.type$haz(pi, link.surv$gradS(eta, etaD), rsurv)
        haz.const <- bhazard + ehaz
      }
      neghaz <- any(haz.const < 0)
      optim.args$kappa <- optim.args$kappa * 10
    }
    optim.args$kappa <- 1
    res_list[[i]] <- res.optim
  }

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
    args$kappa <- 0
    args$x <- res$par
    args$func <- minusloglik
    args$constraint <- FALSE
    hes <- do.call(numDeriv::hessian, args)
    cov <- if (!inherits(vcov <- try(solve(hes)), "try-error"))  vcov
    # cov <- solve(hes)
    if(!is.null(cov) && any(is.na(cov))){
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
            args = args, timeExpr = timeExpr, lm.obj = lm.obj, lm.obj.cr = lm.obj.cr,
            link.type.cr = link.type.cr, link.surv = link.surv, excess = excess,
            timeVar = timeVar, transX = transX, transXD = transXD,
            time = time, event = event, eventExpr = eventExpr, cure.type = cure.type)

  class(L) <- c("gfcm", "cuRe")
  L
}


#Function for computing initial values
get.init <- function(formula, data, smooth.formula, logH.formula, tvc.formula, cr.formula, full.formula, delayed,
                     bhazard, type, link.type.cr, link.surv, timeExpr, time, lm.obj, X, X0, X.cr, method){

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
                          formula.surv = list(formula.k1, ~ 1), type = type)

    #Scale by link function
    pi_hat <- do.call(rbind, predict(fit, type = "curerate", newdata = data))$Estimate

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
    suppressWarnings(fit <- do.call(rstpm2::stpm2, list(formula = formula.2, data = data, bhazard = bhazard)))


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
print.gfcm <- function(fit){
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
summary.gfcm <- function(fit){
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
print.summary.gfcm <- function(x)
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
