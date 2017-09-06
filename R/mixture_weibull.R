
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
#' @param hes Logical variable determining whether to compute the hessian or not. Default is \code{TRUE}
#' @return An object of class \code{CureModel}.
#' @details The arguments for modelling the parameters of the cure model have different meanings dependent on the chosen distribution. \cr
#' For the exponential distribution, k1 denotes the rate, for the weibull model, k1 denotes the scale parameter and k2 denotes the shape parameter, sfor the log normal distribution k1 denotes the mu and sigma.
#' @export


MixtureCureModel <- function(formula, data, bhazard, formula.k1 = ~ 1, formula.k2 = NULL,
                             formula.k3 = NULL, type = "mixture",
                             dist = "weibull", link = "logistic",
                             init = NULL, hes = TRUE){
  all_vars <- c(all.vars(formula), all.vars(formula.k1), all.vars(formula.k2), all.vars(formula.k3))
  ccs <- complete.cases(data[, all_vars])
  data <- data[ccs,]
  fu <- eval(formula[[2]][[2]], envir = data)
  status <- eval(formula[[2]][[3]], envir = data)
  bhazard <- data[, bhazard]


  link_fun <- get_link(link)
  if(dist == "weibull" & is.null(formula.k2)){
    formula.k2 <- ~ 1
  }


  X <- model.matrix(formula, data = data)
  X.k1 <- model.matrix(formula.k1, data = data)

  if(!is.null(formula.k2)){
    X.k2 <- model.matrix(formula.k2, data = data)
  }else{
    X.k2 <- data.frame()
  }

  if(!is.null(formula.k3)){
    X.k3 <- model.matrix(formula.k3, data = data)
  }else{
    X.k3 <- data.frame()
  }

  surv_fun <- get_surv(dist)
  dens_fun <- get_dens(dist)


  init <- rep(0, ncol(X) + ncol(X.k1) + ncol(X.k2) + ncol(X.k3))
  names(init) <- rep(c("gamma", "k1", "k2", "k3"), c(ncol(X), ncol(X.k1), ncol(X.k2), ncol(X.k3)))

  X_all <- list(X, X.k1, X.k2, X.k3)
  if(type == "mixture"){
    minuslog_likelihood <- mixture_minuslog_likelihood
  }else{
    minuslog_likelihood <- nmixture_minuslog_likelihood
  }

  res <- optim(init, fn = mixture_minuslog_likelihood,
               time = fu,
               status = status,
               X = X_all,
               link = link_fun,
               surv_fun = surv_fun,
               dens_fun = dens_fun,
               bhazard = bhazard,
               control = list(maxit = 10000))

  if(res$convergence != 0){
    warning("Convergence not reached")
  }

  if(hes){
    cov <- solve(pracma::hessian(mixture_minuslog_likelihood, res$par,
                                 time = fu,
                                 status = status,
                                 X = X_all,
                                 link = link_fun,
                                 surv_fun = surv_fun,
                                 dens_fun = dens_fun,
                                 bhazard = bhazard))
  }else{
    cov <- NULL
  }

  coefs <- vector("list", 4)
  tmp_coef <- c("gamma", paste0("k", 1:3))
  for(i in 1:length(coefs)){
    coefs[[i]] <- res$par[grepl(tmp_coef[i], names(res$par))]
    names(coefs[[i]]) <- colnames(X_all[[i]])
  }
  L <- list(data = data, formulas = list(formula.gamma = formula, formula.k1 = formula.k1, formula.k2 = formula.k2,
                                         formula.k3 = formula.k3),
            coefs = coefs, dist = dist, link = link, type = type,
            ML = res$value, cov = cov, df = nrow(data) - length(res$par))
  class(L) <- "CureModel"
  L
}


print.CureModel <- function(fit){
  cat("Call:\n")
  print(fit$formula)
  cat("\nCoefficients:\n")
  print(fit$coefs)
}

# MixtureCureModel2 <- function(formula, data, bhazard, formula.k1 = ~ 1, formula.k2 = NULL,
#                              formula.k3 = NULL, type = "mixture",
#                              dist = "weibull", link = "logistic",
#                              init = NULL, hes = TRUE){
#   fu <- eval(formula[[2]][[2]], envir = data)
#   status <- eval(formula[[2]][[3]], envir = data)
#   #od <- order(data[, fu])
#   #data <- data[od,]
#   #bhazard <- bhazard[od]
#
#   if(dist == "weibull" & is.null(formula.k2)){
#     formula.k2 <- ~ 1
#   }
#
#   link_fun <- get_link(link)
#
#   X <- model.matrix(formula, data = data)
#   X.k1 <- model.matrix(formula.k1, data = data)
#
#   if(!is.null(formula.k2)){
#     X.k2 <- model.matrix(formula.k2, data = data)
#   }else{
#     X.k2 <- data.frame()
#   }
#
#   if(!is.null(formula.k3)){
#     X.k3 <- model.matrix(formula.k3, data = data)
#   }else{
#     X.k3 <- data.frame()
#   }
#
#   surv_fun <- get_surv(dist)
#   dens_fun <- get_dens(dist)
#
#
#   fit_glm <- glm(status ~ -1 + X, data = data, family = binomial(link = "logit"))
#
#   formula.2 <- Surv(FU_years, status) ~ age_years + sex
#   fit_cox <- coxph(formula.2, data = data[status == 1, ])
#   cum_haz <- get_basehaz(fit_cox)
#   logH <- log(cum_haz$hazard)
#   fit_lm <- lm(logH ~ -1 + X.k1[status == 1,])
#
#   init <- c(fit_glm$coefficients)
#
#
#   init <- c(fit_glm$coefficients, rep(0, ncol(X.k1) + ncol(X.k2) + ncol(X.k3)))
#   names(init) <- rep(c("gamma", "k1", "k2", "k3"), c(ncol(X), ncol(X.k1), ncol(X.k2), ncol(X.k3)))
#   X_all <- list(X, X.k1, X.k2, X.k3)
#
#   if(type == "mixture"){
#     minuslog_likelihood <- mixture_minuslog_likelihood
#   }else{
#     minuslog_likelihood <- nmixture_minuslog_likelihood
#   }
#
#   res <- optim(init, fn = mixture_minuslog_likelihood,
#                time = fu,
#                status = status,
#                X = X_all,
#                link = link_fun,
#                surv_fun = surv_fun,
#                dens_fun = dens_fun,
#                bhazard = bhazard,
#                control = list(maxit = 10000, reltol = 1e-16))
#
#   if(res$convergence != 0){
#     warning("Convergence not reached")
#   }
#
#   if(hes){
#     cov <- solve(pracma::hessian(mixture_minuslog_likelihood, res$par,
#                                  time = fu,
#                                  status = status,
#                                  X = X_all,
#                                  link = link_fun,
#                                  surv_fun = surv_fun,
#                                  dens_fun = dens_fun,
#                                  bhazard = bhazard))
#   }else{
#     cov <- NULL
#   }
#
#   coefs <- vector("list", 4)
#   tmp_coef <- c("gamma", paste0("k", 1:3))
#   for(i in 1:length(coefs)){
#     coefs[[i]] <- res$par[grepl(tmp_coef[i], names(res$par))]
#     names(coefs[[i]]) <- colnames(X_all[[i]])
#   }
#   L <- list(data = data, formulas = list(formula.gamma = formula, formula.k1 = formula.k1, formula.k2 = formula.k2,
#                                          formula.k3 = formula.k3),
#             coefs = coefs, dist = dist, link = link, type = type,
#             ML = res$value, cov = cov, df = nrow(data) - length(res$par))
#   class(L) <- "CureModel"
#   L
# }
