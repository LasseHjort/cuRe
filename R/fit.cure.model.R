#' Parametric cure model
#'
#' This function fits parametric cure models using simple parametric distributions.
#'
#' @param formula Formula for modelling the cure proportion. The left hand side
#' has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variable names in \code{formula} and \code{formula.surv}.
#' @param bhazard Background hazard.
#' @param formula.surv List of formulas for each parameter in the parametric distribution (see details).
#' @param type A character indicating the type of cure model.
#' Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param dist The parametric distribution of the survival of the uncured.
#' @param link Character. Specifies the link function of the cure proportion.
#' @param covariance Logical. If \code{TRUE} (default), the covariance matrix is computed.
#' @param link.mix Character. Specifies the link function for the mixture parameter in a
#' weibull-weibull mixture model and weibull-exponential model.\cr Only used when \code{dist = "weiwei"} and \code{dist = "weiexp"}.
#' @param control List of control parameters passed to \code{optim}.
#' @param method Optimization method passed to \code{optim}.
#' @param init Initial values for the maximum likelihood optimization.
#' If not provided, the optimization will start in 0.
#' @return An object of class \code{cm} containing the estimated parameters of the cure model.
#' The appropriate link functions taken on \eqn{\pi} and the \eqn{\theta_i}'s are linear in the covariates corresponding to their respective parameter estimates.
#' @details If \code{type = "mixture"}, the function fits the model,
#' \deqn{S(t|z) = \pi(z) + [1 - \pi(z)] S_u(t|z),}
#' and if \code{type = "nmixture"}, the function fits the model,
#' \deqn{S(t|z) = \pi(z)^{\widetilde F(t)},}
#' where z is a vector of covariates. The \code{formula.surv} argument is used to model
#' \eqn{S_u(t)} (1 - \eqn{\widetilde F(t)}). It is a \code{list} of formulas with as many entries as there are
#' parameters in the chosen parametric distribution. If not specified, all formulas are assumed to be \code{~1}.
#' The ith formula, i.e., \code{formula.surv[[i]]} refers to \eqn{\theta_i} in the below survival functions.\cr\cr
#' Exponential model:
#' \deqn{S_u(t) = \exp\left(-t \theta_1\right).}
#' Weibull model:
#' \deqn{S_u(t) = \exp\left(-\theta_1 t^{\theta_2}\right).}
#' Log-normal model:
#' \deqn{S_u(t) = 1 - \Phi\left(\frac{\log(t) - \theta_1}{\theta_2}\right).}
#' Weibull-exponential mixture model:
#' \deqn{S_u(t) = \theta_1\exp\left(-\theta_2 t^{\theta_3}\right) + (1 - \theta_1)\exp\left(-\theta_4 t\right).}
#' Weibull-Weibull mixture model:
#' \deqn{S_u(t) = \theta_1\exp\left(-\theta_2 t^{\theta_3}\right) + (1 - \theta_1)\exp\left(-\theta_4 t^{\theta_5}\right).}
#' Generalized modified Weibull distribution:
#' \deqn{S_u(t) = 1-\left(1 - \exp\left(-\theta_1 t ^ \theta_2 \exp(\theta_3 t)\right)\right) ^ \theta_4.}
#' In the Weibull-exponential and Weibull-Weibull mixture models, the link function for the mixture component is controlled by \code{link.mix}.
#' The remaining parameters are modelled using an exponential link function except \eqn{\theta_1} in the log-normal model,
#' which is modelled using the identity.
#' @export
#' @example inst/fit.cure.model.ex.R


fit.cure.model <- function(formula, data, formula.surv = NULL, type = c("mixture", "nmixture"),
                           dist = c("weibull", "exponential", "lognormal", "weiwei", "weiexp", "gmw"),
                           link = c("logit", "loglog", "identity", "probit"), bhazard = NULL,
                           covariance = TRUE, link.mix = c("logit", "loglog", "identity", "probit"),
                           control = list(maxit = 10000), method = "Nelder-Mead", init = NULL){

  type <- match.arg(type)
  dist <- match.arg(dist)
  link <- match.arg(link)
  link.mix <- match.arg(link.mix)

  max.formulas <- switch(dist,
                         weibull = 2,
                         exponential = 1,
                         lognormal = 2,
                         weiwei = 5,
                         weiexp = 4,
                         gmw = 4)

  if(is.null(formula.surv)){
    formula.surv <- rep(list(~1), max.formulas)
  } else {
    if(length(formula.surv) != max.formulas){
      add.formulas <- rep(list(~1), max.formulas - length(formula.surv))
      formula.surv <- c(formula.surv, add.formulas)
    }
  }

  #Delete missing observations and extract response data
  all.formulas <- c(formula, formula.surv)
  all_vars <- unique(unlist(lapply(all.formulas, all.vars)))
  # data.c <- data[complete.cases(data[, all_vars]),]
  data.c <- data

  #Extract survival time and event variable
  eventExpr <- lhs(formula)[[length(lhs(formula))]]
  delayed <- length(lhs(formula)) >= 4
  timeExpr <- lhs(formula)[[ifelse(delayed, 3, 2)]]
  timeVar <- all.vars(timeExpr)
  time <- eval(timeExpr, data, parent.frame())

  time0Expr <- NULL # initialise
  if (delayed) {
    time0Expr <- lhs(formula)[[2]]
    time0Var <- all.vars(time0Expr)
    time0 <- eval(time0Expr, data, parent.frame())
  } else {
    time0 <- rep(0, nrow(data))
  }
  ind0 <- time0 > 0

  event <- eval(eventExpr, data)
  event <- if (length(unique(event)) == 1){
    rep(TRUE, length(event))
  } else {
    event <- event > min(event)
  }

  #Create background hazard
  if(is.null(bhazard)){
    bhazard <- rep(0, nrow(data))
  }else {
    if(!is.numeric(bhazard)){
      bhazard <- data[, bhazard]
    }
  }

  excess <- ifelse(any(bhazard != 0), T, F)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "contrasts", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  lm.objects <- lapply(all.formulas, function(formula){
    lm.call <- mf
    lm.call[[1L]] <- as.name("lm")
    lm.formula <- formula
    lhs(lm.formula) <- quote(arbri)
    lm.call$formula <- lm.formula
    data$arbri <- rnorm(nrow(data))
    lm.call$data <- quote(data)
    eval(lm.call)
  })

  #Compute design matrices
  X.all <- lapply(lm.objects, lpmatrix.lm, newdata = data)
  #X.all <- lapply(all.formulas, get_design, data = data)

  #Extract link functions
  link.fun <- get.link(link)
  surv.fun <- get.surv(dist, link.mix = get.link(link.mix))
  dens.fun <- get.dens(dist, link.mix = get.link(link.mix))

  #Extract likelihood function
  # if(delayed){
  #   minusloglik <- switch(type,
  #                         mixture = mixture_minuslog_likelihoodDelayed,
  #                         nmixture = nmixture_minuslog_likelihoodDelayed)
  # } else {
  #   minusloglik <- switch(type,
  #                         mixture = mixture_minuslog_likelihood,
  #                         nmixture = nmixture_minuslog_likelihood)
  # }

  cure.type <- switch(type,
                      mixture = mix,
                      nmixture = nmix)

  minusloglik <- minuslog_likelihoodDelayed
  #minusloglik <- switch(type,
  #                      mixture = mixture_minuslog_likelihoodDelayed,
  #                      nmixture = nmixture_minuslog_likelihoodDelayed)

  #Prepare optimization arguments
  args <- list(time = time,
               time0 = time0,
               event = event,
               Xs = X.all,
               ind0 = ind0,
               link = link.fun,
               surv.fun = surv.fun,
               dens.fun = dens.fun,
               bhazard = bhazard,
               cure.type = cure.type)

  if(is.null(control$maxit)){
    control$maxit <- list(maxit = 10000)
  }

  n.param.formula <- sapply(X.all, ncol)
  #Get initial values
  if(is.null(init)){
    n.param     <- sum(n.param.formula)
    init        <- rep(0, n.param)
    names(init) <- unlist(lapply(X.all, colnames))
  }

  optim.args        <- c(control = list(control), args)
  optim.args$method <- method
  optim.args$fn     <- minusloglik
  optim.args$par    <- init

  #Run optimization
  optim.out <- do.call(optim, optim.args)

  #Check for convergence
  if(optim.out$convergence != 0){
    warning("Convergence not reached")
  }

  #Compute covariance
  if(covariance){
    args$x <- optim.out$par
    args$func <- minusloglik
    hes <- do.call(numDeriv::hessian, args)
    cov <- if (!inherits(vcov <- try(solve(hes)), "try-error"))  vcov
    if(!is.null(cov) && any(is.na(cov))){
      warning("Hessian is not invertible!")
    }
  }else{
    cov <- NULL
  }

  groups <- factor(rep(1:(max.formulas + 1), n.param.formula), 1:(max.formulas + 1))
  coefs <- split(optim.out$par, f = groups)

  #Output list
  L <- list(data = data, all.formulas = all.formulas,
            coefs = coefs, dist = dist, link = link,
            type = type, ci = covariance,
            ML = optim.out$value, covariance = cov,
            df = nrow(data) - length(optim.out$par),
            optim = optim.out, n.param.formula = n.param.formula,
            surv.fun = surv.fun, dens.fun = dens.fun, optim.args = optim.args,
            time = time, event = event, timeVar = timeVar, link.mix = link.mix,
            excess = excess, cure.type = cure.type, args = args, lm.objects = lm.objects)
  class(L) <- c("cm", "cuRe")
  L
}


get_design <- function(formula, data){
  if(!is.null(formula))
    model.matrix(formula, data = data)
  else
    data.frame()
}

#' @export
#' @method print cm
print.cm <- function(fit){
  type <- switch(fit$type,
                 mixture = "mixture",
                 nmixture = "non-mixture")
  cat("Model:\n")
  print(paste0("Parametric ", type, " cure model"))
  cat("Family survival / curerate: \n")
  print(paste0(fit$dist, " / ", fit$link))
  is.not.null <- sapply(fit$formulas, function(f) !is.null(f))
  coef.names <- unlist(lapply(fit$formulas[is.not.null], function(x) Reduce(paste, deparse(x))))
  cat("\nCoefficients:\n")
  coefs <- fit$coefs[sapply(fit$coefs, function(coef) length(coef) != 0)]
  names(coefs) <- coef.names
  print(coefs)
}

replace_names <- function(coef_list){
  label_names <- c("1" = "pi", "2" = "theta1",
                   "3" = "theta2", "4" = "theta3",
                   "5" = "theta4", "6" = "theta5")

  names(coef_list) <- label_names[names(coef_list)]
  coef_list
}

#' @export
#' @method summary cm
summary.cm <- function(fit){
  se <- sqrt(diag(fit$cov))
  coef_list <- fit$coefs
  coef_list <- replace_names(coef_list)
  coefs <- unlist(coef_list)
  z <- coefs / se
  TAB1 <- cbind(Estimate = coefs,
                StdErr = se,
                z.value = z,
                p.value = ifelse(is.na(z), rep(NA, length(coefs)),
                                 2 * pnorm(-abs(z))))


  results <- list(coefs = TAB1)
  results$type <- fit$type
  results$link <- fit$link
  results$ML <- fit$ML
  formulas <- fit$all.formulas
  names(formulas) <- c("gamma", paste0("k", 1:(length(formulas) - 1)))
  results$formulas <- formulas[sapply(formulas, function(x) !is.null(x))]
  class(results) <- "summary.cm"
  results
}

#' @export
#' @method print summary.cm
print.summary.cm <- function(x)
{
  cat("Calls:\n")
  print(x$formulas)
  #    cat("\n")
  printCoefmat(x$coefs, P.values = TRUE, has.Pvalue = T)
  cat("\n")
  cat("Type =", x$type, "\n")
  cat("Link =", x$link, "\n")
  cat("LogLik(model) =", x$ML, "\n")
}
