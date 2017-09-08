get_design <- function(formula, data){
  if(!is.null(formula))
    model.matrix(formula, data = data)
  else
    data.frame()
}

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
#' @return An object of class \code{CureModel}.
#' @details The arguments for modelling the parameters of the cure model have different meanings dependent on the chosen distribution. \cr
#' For the exponential distribution, k1 denotes the rate, for the weibull model, k1 denotes the scale parameter and k2 denotes the shape parameter, sfor the log normal distribution k1 denotes the mu and sigma.
#' @export


fit.cure.model <- function(formula, data, bhazard, formula.k1 = ~ 1, formula.k2 = NULL,
                           formula.k3 = NULL, type = "mixture",
                           dist = "weibull", link = "logit",
                           init = NULL, covariance = TRUE){

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
  link_fun <- get_link(link)
  surv_fun <- get_surv(dist)
  dens_fun <- get_dens(dist)

  #Set initial values for the optimization
  n.param.formula <- sapply(X.all, ncol)
  n.param <- sum(n.param.formula)
  init <- rep(0, n.param)
  names(init) <- unlist(lapply(X.all, colnames))

  #Extract likelihood function
  minuslog_likelihood <- switch(type,
                                mixture = mixture_minuslog_likelihood,
                                nmixture = nmixture_minuslog_likelihood)

  #Run optimization
  optim.out <- optim(init,
                     fn = minuslog_likelihood,
                     time = time,
                     status = event,
                     Xs = X.all,
                     link = link_fun,
                     surv_fun = surv_fun,
                     dens_fun = dens_fun,
                     bhazard = bhazard.eval,
                     control = list(maxit = 10000), hessian = T)

  #Check for convergence
  if(optim.out$convergence != 0){
    warning("Convergence not reached")
  }

  #Compute covariance
  if(covariance){
    cov <- solve(pracma::hessian(minuslog_likelihood,
                                 optim.out$par,
                                 time = time,
                                 status = event,
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
            ML = res$value, cov = cov,
            df = nrow(data) - length(res$par),
            optim = optim.out, n.param.formula = n.param.formula,
            surv_fun = surv_fun, dens_fun = dens_fun)
  class(L) <- c("curemodel", "cuRe")
  L
}


print.CureModel <- function(fit){
  cat("Call:\n")
  print(fit$formula)
  cat("\nCoefficients:\n")
  print(fit$coefs)
}
