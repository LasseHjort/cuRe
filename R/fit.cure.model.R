#' Parametric cure model
#'
#' This function is used to fit parametric cure models on the relative survival.
#'
#' @param formula Formula for modelling the cure proportion. Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variable names in \code{formula} and \code{formula.surv}.
#' @param bhazard Background hazard.
#' @param formula.surv List of formulas for each parameter in the parametric distribution (see details).
#' @param type A character indicating the type of cure model.
#' Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param dist The parametric distribution of the survival of the uncured.
#' @param link Character. Specifies the link function of the cure proportion.
#' @param covariance Logical. If \code{TRUE} (default), the covariance matrix is computed.
#' @param link.mix Character. Specifies the link function for the mixture parameter in a
#' weibull-weibull mixture model and weibull-exponential model. Only used when \code{dist = "weiwei"} and \code{dist = "weiexp"}.
#' @param control List of control parameters passed to \code{optim}.
#' @param method Optimization method passed to \code{optim}.
#' @param init Initial values for the maximum likelihood optimization.
#' If not provided, the optimization will start in 0.
#' @return An object of class \code{cm} containing the parameters of the cure model.
#' @details The function fits the model,
#' \deqn{S(t) = \pi + (1 - \pi) S_u(t).}
#' The \code{formula.surv} argument is used to model \eqn{S_u(t)}. It is a \code{list} of formulas with as many entries as there are
#' parameters in the chosen parametric distribution. If not specified, all formulas are assumed to be \code{~1}.
#' The ith formula, i.e., \code{formula.surv[[i]]} refers \eqn{\theta_i} in the below survival functions.\cr\cr
#' Exponential model:
#' \deqn{S_u(t) = \exp\left(-t \theta_1\right).}
#' Weibull model:
#' \deqn{S_u(t) = \exp\left(-\theta_1 t^{\theta_2}\right).}
#' Log-normal model:
#' \deqn{S_u(t) = 1 - \Phi\left(\frac{\log(t) - \theta_1}{\theta_2}\right)}
#' Weibull-exponential mixture model:
#' \deqn{S_u(t) = \theta_1\exp\left(-\theta_2 t^{\theta_3}\right) + (1 - \theta_1)\exp\left(-\theta_4 t\right).}
#' Weibull-Weibull mixture model:
#' \deqn{S_u(t) = \theta_1\exp\left(-\theta_2 t^{\theta_3}\right) + (1 - \theta_1)\exp\left(-\theta_4 t^{\theta_5}\right).}
#' In the the mixture models, the link function for the mixture component is controlled by \code{link.mix}.
#' The remaining parameters are modelled using an exponential link function except \eqn{\theta_1} in the log-normal model,
#' which is modeled using the identity.
#' @export
#' @example inst/fit.cure.model.ex.R


fit.cure.model <- function(formula, data, bhazard = NULL, formula.surv = NULL, type = c("mixture", "nmixture"),
                           dist = c("weibull", "exponential", "lognormal", "weiwei", "weiexp"),
                           link = c("logit", "loglog", "identity", "probit"),
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
                         weiexp = 4)

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
  all.vars <- unique(unlist(lapply(all.formulas, all.vars)))
  data.c <- data[complete.cases(data[, all.vars]),]

  #Extract survival time and event variable
  eventExpr <- rstpm2:::lhs(formula)[[length(rstpm2:::lhs(formula))]]
  delayed <- length(rstpm2:::lhs(formula)) >= 4
  timeExpr <- rstpm2:::lhs(formula)[[ifelse(delayed, 3, 2)]]
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

  #Compute design matrices
  X.all <- lapply(all.formulas, get_design, data = data)

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

  minusloglik <- switch(type,
                        mixture = mixture_minuslog_likelihoodDelayed,
                        nmixture = nmixture_minuslog_likelihoodDelayed)

  #Prepare optimization arguments
  args <- list(time = time,
               time0 = time0,
               event = event,
               Xs = X.all,
               link = link.fun,
               surv.fun = surv.fun,
               dens.fun = dens.fun,
               bhazard = bhazard)

  if(is.null(control$maxit)){
    control$maxit <- list(maxit = 10000)
  }

  n.param.formula <- sapply(X.all, ncol)
  #Get initial values
  if(is.null(init)){
    n.param <- sum(n.param.formula)
    init <- rep(0, n.param)
    names(init) <- unlist(lapply(X.all, colnames))
  }

  optim.args <- c(control = list(control), args)
  optim.args$method <- method
  optim.args$fn <- minusloglik
  optim.args$par <- init

  #Run optimization
  optim.out <- do.call(optim, optim.args)

  #Check for convergence
  if(optim.out$convergence != 0){
    warning("Convergence not reached")
  }

  #Compute covariance
  if(covariance){
    args$x0 <- optim.out$par
    args$f <- minusloglik
    hes <- do.call(pracma::hessian, args)
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
            time = time, event = event, timeVar = timeVar, link.mix = link.mix, excess = excess)
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

#' @export
summary.cm <- function(fit){
  se <- sqrt(diag(fit$cov))
  coefs <- unlist(fit$coefs)
  tval <- coefs / se
  TAB1 <- cbind(Estimate = coefs,
                StdErr = se,
                t.value = tval,
                p.value = ifelse(is.na(tval), rep(NA, length(coefs)),
                                 2 * pt(-abs(tval), df = fit$df)))

  results <- list(coefs = TAB1)
  results$type <- fit$type
  results$link <- fit$link
  results$ML <- fit$ML
  formulas <- fit$formulas
  names(formulas) <- c("gamma", "k1", "k2", "k3")
  results$formulas <- formulas[sapply(formulas, function(x) !is.null(x))]
  class(results) <- "summary.cm"
  results
}

#' @export
print.summary.cm <- function(x)
{
  cat("Calls:\n")
  print(x$formulas)
  #    cat("\n")
  printCoefmat(x$coefs, P.value = TRUE, has.Pvalue = T)
  cat("\n")
  cat("Type =", x$type, "\n")
  cat("Link =", x$link, "\n")
  cat("LogLik(model) =", x$ML, "\n")
}

