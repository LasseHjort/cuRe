#' Fit spline-based mixture cure model
#'
#' The following function fits a generalized mixture or non-mixture cure model using a link function for the cure rate and for the survival of the uncured
#'
#' @param formula Formula for the covariates in the modelling of the cure rate, \eqn{\pi}. Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variables names in \code{formula}, \code{smooth.formula}, and \code{tvc.formula}.
#' @param bhazard Background hazard
#' @param smooth.formula Formula to control the modelling of the disease-specific survival of the uncured.
#' @param knots Knots used for the baseline hazard in the disease-specific survival function
#' @param n.knots Number of knots for the disease-specific survival function. The knots are calculated as the equidistant quantiles of the uncensored event-times.
#' If \code{knots} are supplied, this argument will be ignored.
#' @param tvc.formula Formula for the time-varying effects.
#' @param knots.time An object of class list containing the knots for each of the covariates in the time-varying covariate effects
#' @param n.knots.time An object of class list, containing the number of knots for the time-varying covariate effect.
#' If \code{knots.time} is supplied, this argument will be ignored.
#' @param hes Logical for computing the inverse hessian matrix (default is \code{TRUE}).
#' @param verbose Logical indicating whether to output messages from the function
#' @param type Character indicating which type of model is fitting.
#' Possible values are \code{mixture} (default) and \code{nmixture}.
#' @param linkpi Character for the type of link function selected for the cure rate.
#' Possible values are \code{logistic} (default), \code{identity}, \code{loglog}, and \code{probit}.
#' @param linksu Character for the type of link function selected for the survival of the uncured.
#' Possible values are \code{loglog}, \code{logistic}, and \code{probit}.
#' @param optim.args List with additional arguments to pass to \code{optim}.
#' @param ini.types Character vector denoting which procedures for calculating initial values has to be executed.
#' @return An object of class \code{fmcm}.
#' @export
#' @import survival
#' @import rstpm2
#' @example inst/FlexMixtureCureModel.ex.R



FlexMixtureCureModel <- function(formula, data, bhazard, smooth.formula = ~ 1,
                                 knots = NULL, n.knots = NULL,
                                 tvc.formula = NULL, knots.time = NULL, n.knots.time = NULL,
                                 covariance = T, type = "mixture", linkpi = "logit",
                                 linksu = "loglog", message = T, optim.args = NULL,
                                 ini.types = c("cure", "flexpara")){

  #Extract relevant variables
  fu <- eval(formula[[2]][[2]], envir = data)
  status <- eval(formula[[2]][[3]], envir = data)
  death_times <- fu[status == 1]

  #Caculate placement of knots and establish basis matrices
  if(is.null(knots)){
    bd_knots <- log(range(death_times))
    inner_knots <- log(quantile(death_times, 1 / (n.knots - 1)*1:(n.knots - 2)))
    knots <- sort(c(bd_knots, inner_knots))
  }else{
    knots <- sort(knots)
    bd_knots <- range(knots)
    inner_knots <- knots[-c(1, length(knots))]
  }

  b <- flexsurv::basis(knots = knots, log(fu))
  db <- flexsurv::dbasis(knots = knots, log(fu))

  if(!is.null(tvc.formula)){
    if(is.null(knots.time)){
      vars <- all.vars(tvc.formula)
      knots.time <- lapply(vars, function(x){
        bd_knots <- range(death_times)
        if(n.knots.time[[x]] > 2){
          inner_knots <- quantile(death_times, 1 / (n.knots.time[[x]] - 1)*1:(n.knots.time[[x]] - 2))
          log(sort(c(bd_knots, inner_knots)))
        }else{
          log(bd_knots)
        }
      })
      #names(knots.time) <- names(n.knots.time)
      b_list <- lapply(knots.time, flexsurv::basis, x = log(fu))
      db_list <- lapply(knots.time, flexsurv::dbasis, x = log(fu))
    }
  }else{
    tvc.formula <- ~ 1
  }

  X_time <- model.matrix(tvc.formula, data = data)[,-1, drop = FALSE]
  if(ncol(X_time) > 0){
    for(i in 1:ncol(X_time)){
      b_list[[i]] <- b_list[[i]] * X_time[,i]
      db_list[[i]] <- db_list[[i]] * X_time[,i]
    }
    b_time <- do.call(cbind, b_list)
    db_time <- do.call(cbind, db_list)
  }else{
    b_time <- NULL
    db_time <- NULL
    tvc.formula <- NULL
  }


  #Construct design matrix
  X <- model.matrix(smooth.formula, data = data)
  b <- cbind(b, X[, -1], b_time)
  db <- cbind(db, matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)), db_time)
  X <- model.matrix(formula, data = data)

  #Extract link function
  link_fun_pi <- get.link(linkpi)
  link_fun_su <- get.link(linksu)
  dlink_fun_su <- get.dlink(linksu)

  #Extract minus log likelihood function
  minusloglik <- switch(type,
                        mixture = flexible_mixture_minuslog_likelihood,
                        nmixture = flexible_nmixture_minuslog_likelihood)

  #Prepare optimization arguments
  likelihood.pars <- list(fn = minusloglik, time = fu,
                          status = status, X = X,
                          b = b, db = db, bhazard = data[, bhazard],
                          link_fun_pi = link_fun_pi,
                          link_fun_su = link_fun_su,
                          dlink_fun_su = dlink_fun_su)

  if(is.null(optim.args$control$maxit)){
    optim.args$control <- list(maxit = 10000)
  }

  #Generate initial values if these are not provided by the user
  if(is.null(optim.args$par)){
    if(message) cat("Finding initial values... ")
    inivalues <- lapply(ini.types, function(ini.type) get_ini_values(smooth.formula = smooth.formula,
                                                          tvc.formula =  tvc.formula,
                                                          data = data,
                                                          bhazard = bhazard,
                                                          linkpi = linkpi,
                                                          linksu = linksu,
                                                          formula = formula,
                                                          type = type,
                                                          fu = fu,
                                                          status = status,
                                                          n.knots.time = n.knots.time,
                                                          knots = knots,
                                                          X = X, b = b,
                                                          method = ini.type))
  }else{
    if(message) cat("Initial values provided by the user... ")
    inivalues <- optim.args$par
    optim.args <- optim.args[-which(names(optim.args) == "par")]
  }

  optim.pars <- c(optim.args, likelihood.pars)

  #Test if initial values are within the feasible region
  ini.eval <- sapply(inivalues, function(inival) do.call(minusloglik, c(likelihood.pars[-1], list(inival))))
  run.these <- !is.na(ini.eval)

  if(message) cat("Completed!\nFitting the model... ")

  #Fit each model
  res_list <- lapply(inivalues[run.these], function(inival){
    optim.pars$par <- inival
    do.call(optim, optim.pars)
  })

  #Choose the best model according to the maximum likelihood estimate
  MLs <- sapply(res_list, function(x) x$value)
  wh <- which.min(MLs)
  res <- res_list[[wh]]

  #names(res$par) <- gsub("spline_", "", names(res$par))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(message) cat("Completed!\n")

  #Compute the covariance matrix matrix
  if(covariance){
    cov <- solve(pracma::hessian(minusloglik, res$par,
                                 time = fu, status = status,
                                 X = X, b = b, db = db,
                                 bhazard = data[,bhazard],
                                 link_fun_pi = link_fun_pi,
                                 link_fun_su = link_fun_su,
                                 dlink_fun_su = dlink_fun_su))
  }else{
    cov <- NULL
  }

  #Output the results
  L <- list(formula = formula,
            data = data,
            coefs = res$par[1:ncol(X)],
            coefs.spline = res$par[(ncol(X) + 1):length(res$par)],
            knots = knots, knots.time = knots.time,
            NegMaxLik = res$value, covariance = cov, tvc.formula = tvc.formula, formula = formula,
            formula_main = smooth.formula, type = type,
            link_fun_pi = link_fun_pi,
            link_fun_su = link_fun_su,
            dlink_fun_su = dlink_fun_su,
            linkpi = linkpi, linksu = linksu,
            df = length(res$par) - 1, NegMaxLiks = MLs, optim.pars = optim.pars,
            times = fu)

  class(L) <- c("fcm", "cuRe")
  L
}



#Function for computing initial values
get_ini_values <- function(smooth.formula, tvc.formula, data, bhazard, linkpi, linksu, formula, type,
                           fu, status, n.knots.time, knots, X, b, method = "cure"){
  vars <- c(all.vars(smooth.formula), all.vars(tvc.formula))
  if(method == "cure"){
    formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                             response = NULL)
    fit <- fit.cure.model(formula, data = data, bhazard = bhazard, covariance = F,
                          formula.k1 = formula.2, formula.k2 = ~ 1, type = type)
    pi_hat <- get.link("logit")(X %*% fit$coefs[[1]])
    gpi_hat <- get.inv.link(linkpi)(pi_hat)
    ini_pi <- lm(gpi_hat ~ -1 + X)$coefficients
    lp <- exp(model.matrix(formula.2, data = data) %*% fit$coefs[[2]])
    shat <- exp(-lp * fu ^ exp(fit$coefs[[3]]))
    gshat <- get.inv.link(linksu)(shat)
    fit_lm <- lm(gshat ~ -1 + b)
  }else if(method == "deaths"){
    formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                             response = formula[[2]])
    status2 <- 1 - status
    fit_glm <- glm(status2 ~ -1 + X, family = binomial(link = "logit"))
    pi_hat <- get.link("logit")(predict(fit_glm))
    gpi_hat <- get.inv.link(linkpi)(pi_hat)
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    fit <- coxph(formula.2, data = data[status == 1,])
    cum_base_haz <- get_basehaz(fit)
    shat <- exp(-cum_base_haz$hazard) ^ exp(fit$linear.predictors)
    suppressWarnings(gshat <- get.inv.link(linksu)(shat))
    fit_lm <- lm(gshat ~ -1 + b[status == 1,])
  }else if(method == "flexpara"){
    tt <- terms(formula)
    formula.2 <- formula(delete.response(tt))
    vars <- unique(c("-1", all.vars(formula.2), all.vars(smooth.formula), all.vars(tvc.formula)))
    formula.3 <- reformulate(termlabels = vars, response = formula[[2]])
    fu_time <- all.vars(formula.3)[1]
    smooth.formula.paste <- as.formula(paste0("~basis(knots = knots, x = log(", fu_time, "))"))
    #fit <- stpm2(formula.3, data = data, smooth.formula = smooth.formula.paste,
    #             bhazard = data[, bhazard])

    fit <- stpm2(Surv(FU_years, status) ~ 1, data = data, df = length(knots) - 1, bhazard = data[, bhazard])
    shat <- predict(fit, newdata = data, se.fit = F)
    gshat <- get.inv.link(linksu)(shat)
    data2 <- data
    data2[, fu_time] <- max(data2[, fu_time]) + 0.1
    pi_hat <- predict(fit, newdata = data2, se.fit = F)
    wh <- which(pi_hat >= shat)
    pi_hat[wh] <- shat[wh] - 0.01
    gpi_hat <- get.link(linkpi)(pi_hat)
    formula.logistic <- reformulate(termlabels = if(length(all.vars(formula.2)) == 0) "1" else all.vars(formula.2),
                                    response = "pred_pi")
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    suhat <- (shat - pi_hat) / (1 - pi_hat)
    gsuhat <- get.inv.link(linksu)(suhat)
    fit_lm <- lm(gsuhat ~ -1 + b, data = data)
  }
  #Make proper naming of the initial values
  suppressWarnings(coefs <- summary(fit_lm)$coefficients[, "Estimate"])
  names(coefs)[1:length(knots)] <- paste0("gamma", 1:length(knots))
  fix_var <- all.vars(smooth.formula)
  if(length(fix_var) != 0){
    if(length(coefs) > length(knots)){
      names(coefs)[(length(knots) + 1):(length(knots) + length(fix_var))] <- paste0("spline_", fix_var)
    }
    if(length(coefs) > (length(knots) + length(fix_var))){
      names(coefs)[(length(knots) + length(fix_var) + 1):length(coefs)] <- paste0("gamma", 1:n.knots.time[[1]], ":",
                                                                                  all.vars(tvc.formula))
    }
  }else{
    if(length(coefs) > length(knots)){
      names(coefs)[(length(knots) + 1):length(coefs)] <- unlist(lapply(all.vars(tvc.formula),
                                                                       function(x) paste0("gamma",
                                                                                          1:n.knots.time[[x]],
                                                                                          ":", x)))
    }
  }
  ini_values <- c(ini_pi, coefs)
  names(ini_values)[1:ncol(X)] <- colnames(X)
  ini_values
}


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
  class(results) <- "summary.fmcm"
  results
}

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

