#' Parametric cure model
#'
#' This function is used to fit parametric cure models on the relative survival.
#'
#' @param formula Formula for modelling the the cure fraction Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variables names in \code{formula},
#' \code{formula.k1}, \code{formula.k2}, and \code{formula.k3}.
#' @param bhazard Background hazard.
#' @param formula.k1 Formula for the first linear predictor of the parametric distribution (see details).
#' @param formula.k2 Formula for the second linear predictor of the parametric distribution (see details).
#' @param formula.k3 Formula for the third linear predictor of the parametric distribution (see details).
#' @param type A character indicating the type of cure model.
#' Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param dist The parametric distribution of the disease-specific survival function. Possible values are
#' \code{weibull} (default), \code{exponential}, and \code{lognormal}.
#' @param link Specifies the link function of the cure fraction.
#' Possible values are \code{logit} (default), \code{identity}, and \code{loglog}.
#' @param init Initial values for the maximum likelihood optimization.
#' If not provided, the optimization will start in 0.
#' @param covariance Logical. If \code{TRUE} (default), the covariance matrix is computed.
#' @param optim.args List with additional arguments passed to \code{optim}.
#' @return An object of class \code{CureModel}.
#' @details The arguments for modelling the parameters of the cure model have different meanings dependent on the chosen distribution. \cr
#' For the exponential distribution, k1 denotes the rate, for the weibull model, k1 denotes the scale parameter and k2 denotes the shape parameter, sfor the log normal distribution k1 denotes the mu and sigma.
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
  time <- eval(timeExpr, data, parent.frame())
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
  minusloglik <- switch(type,
                        mixture = mixture_minuslog_likelihood,
                        nmixture = nmixture_minuslog_likelihood)

  #Prepare optimization arguments
  args <- list(time = time,
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
            times = time, link.mix = link.mix, excess = excess)
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

