#' Parametric cure model
#'
#' This function is used to fit parametric cure models on the relative survival.
#'
#' @param formula An object of type formula for the cure rate.
#' @param data The provided data. Only objects of class \code{data.frame} is allowed.
#' @param bhazard The expected survival rate at the last follow-up of the individual.
#' @param formula.k1 Formula for the first parameter of the parametric distribution (see details).
#' @param formula.k2 Formula for the second parameter of the parametric distribution (see details).
#' @param formula.k3 Formula for the third parameter of the parametric distribution (see details).
#' @param type Determines the type of cure model. Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param dist The parametric distribution of the disease-specific survival function. Possible values are\cr
#' \code{weibull} (default), \code{exponential}, and \code{lognormal}.
#' @param link Specifies the link function of the cure rate. Possible values are \code{logistic} (default), \code{identity}, and \code{loglog}.
#' @param init Initial values for the maximum likelihood optimization
#' @param covariance Logical variable determining whether to compute the hessian or not. Default is \code{TRUE}
#' @param optim.args List with additional arguments to pass to \code{optim}.
#' @return An object of class \code{CureModel}.
#' @details The arguments for modelling the parameters of the cure model have different meanings dependent on the chosen distribution. \cr
#' For the exponential distribution, k1 denotes the rate, for the weibull model, k1 denotes the scale parameter and k2 denotes the shape parameter, sfor the log normal distribution k1 denotes the mu and sigma.
#' @export
#' @example inst/fit.cure.model.ex.R


fit.cure.model <- function(formula, data, bhazard, formula.k1 = ~ 1, formula.k2 = NULL,
                           formula.k3 = NULL, type = "mixture",
                           dist = "weibull", link = "logit",
                           covariance = TRUE,
                           optim.args = NULL){

  #Delete missing observations and extract response data
  formulas <- list(formula, formula.k1, formula.k2, formula.k3)
  all.vars <- unique(unlist(lapply(formulas, all.vars)))
  ccs <- complete.cases(data[, all.vars])
  data <- data[ccs,]
  Surv_object <- eval(formula[[2]], envir = data)
  time <- Surv_object[,1]
  event <- Surv_object[, 2]
  bhazard.eval <- data[, bhazard]

  #Check if formula is NULL
  if(dist == "weibull" & is.null(formula.k2)){
    stop("formula.k2 has to be specified for dist = weibull")
  }

  #Compute design matrices
  X.all <- lapply(formulas, get_design, data = data)

  #Extract link functions
  link_fun <- get.link(link)
  surv_fun <- get.surv(dist)
  dens_fun <- get.dens(dist)

  #Extract likelihood function
  minuslog_likelihood <- switch(type,
                                mixture = mixture_minuslog_likelihood,
                                nmixture = nmixture_minuslog_likelihood)

  #Prepare optimization arguments
  likelihood.pars <- list(fn = minuslog_likelihood,
                           time = time,
                           event = event,
                           Xs = X.all,
                           link = link_fun,
                           surv_fun = surv_fun,
                           dens_fun = dens_fun,
                           bhazard = bhazard.eval)

  if(is.null(optim.args$control$maxit)){
    optim.args$control <- list(maxit = 10000)
  }

  #Get initial values
  if(is.null(optim.args$par)){
    types <- c("cure", "flexpara")
    n.param.formula <- sapply(X.all, ncol)
    n.param <- sum(n.param.formula)
    optim.args$par <- rep(0, n.param)
    names(optim.args$par) <- unlist(lapply(X.all, colnames))
  }

  optim.pars <- c(optim.args, likelihood.pars)


  #Run optimization
  optim.out <- do.call(optim, optim.pars)

  #Check for convergence
  if(optim.out$convergence != 0){
    warning("Convergence not reached")
  }

  #Compute covariance
  if(covariance){
    cov <- solve(pracma::hessian(minuslog_likelihood,
                                 optim.out$par,
                                 time = time,
                                 event = event,
                                 X = X.all,
                                 link = link_fun,
                                 surv_fun = surv_fun,
                                 dens_fun = dens_fun,
                                 bhazard = bhazard.eval))
  }else{
    cov <- NULL
  }

  groups <- factor(rep(1:4, n.param.formula), 1:4, labels = c("gamma", "k1", "k2", "k3"))
  coefs <- split(optim.out$par, f = groups)

  #Output list
  L <- list(data = data, formulas = formulas,
            coefs = coefs, dist = dist, link = link,
            type = type,
            ML = optim.out$value, covariance = cov,
            df = nrow(data) - length(optim.out$par),
            optim = optim.out, n.param.formula = n.param.formula,
            surv_fun = surv_fun, dens_fun = dens_fun, optim.pars = optim.pars,
            times = time)
  class(L) <- c("cm", "cuRe")
  L
}

get_design <- function(formula, data){
  if(!is.null(formula))
    model.matrix(formula, data = data)
  else
    data.frame()
}


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

