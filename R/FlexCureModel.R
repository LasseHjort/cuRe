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
#' @example inst/FlexCureModel.ex.R



FlexCureModel <- function(formula, data, bhazard, smooth.formula = ~ 1,
                          knots = NULL, n.knots = NULL,
                          knots.time = NULL, n.knots.time = NULL,
                          covariance = T, type = "mixture", linkpi = "logit",
                          linksu = "loglog", verbose = T, constr.optim = F,
                          optim.args = NULL, ortho = TRUE,
                          ini.types = c("cure", "flexpara")){

  if(!type %in% c("mixture", "nmixture"))
    stop("Wrong specication of argument type, must be either 'mixture' or 'nmixture'")

  #Extract relevant variables
  times <- eval(formula[[2]][[2]], envir = data)
  event <- eval(formula[[2]][[3]], envir = data)
  d.times <- times[event == 1]

  #Caculate placement of knots and establish basis matrices
  if(is.null(knots)){
    bd_knots <- log(range(d.times))
    inner_knots <- log(quantile(d.times, 1 / (n.knots - 1)*1:(n.knots - 2)))
    knots <- sort(c(bd_knots, inner_knots))
  }else{
    knots <- sort(knots)
    bd_knots <- range(knots)
    inner_knots <- knots[-c(1, length(knots))]
  }

  #Evaluate baseline model matrices
  b <- basis(knots = knots, x = log(times), ortho = ortho)
  colnames(b) <- paste0("basis(knots = knots, x = log(times), ortho = ortho)", 1:ncol(b))
  R.inv <- attr(b, "R.inv")
  db <- dbasis(knots = knots, log(times), ortho = ortho, R.inv = R.inv)

  #Calculate time-varying knots
  if( !is.null(n.knots.time) | !is.null(knots.time) ){
    if( is.null(knots.time) ){
      vars <- names(n.knots.time)
      knots.time <- lapply(vars, function(x){
        bd_knots <- range(d.times)
        if( n.knots.time[[x]] > 2 ){
          inner_knots <- quantile(d.times, 1 / (n.knots.time[[x]] - 1)*1:(n.knots.time[[x]] - 2))
          log(sort(c(bd_knots, inner_knots)))
        } else {
          log(bd_knots)
        }
      })
    }
    b_list <- lapply(knots.time, basis, x = log(times), ortho = ortho)
    R.inv_list <- lapply(b_list, attr, "R.inv")
    db_list <- lapply(1:length(knots.time), function(i){
      dbasis(x = log(times), knots = knots.time[[i]], ortho = ortho, R.inv = R.inv_list[[i]])
    })
    vars2 <- c(vars)
    tvc.formula <- as.formula(paste0("~ ", paste0(vars2, collapse = " + ")))
  } else {
    tvc.formula <- ~ 1
    R.inv_list <- NULL
  }

  #Get time-varying design matrices
  X_time <- model.matrix(tvc.formula, data = data)[,-1, drop = FALSE]
  if(ncol(X_time) > 0){
    for(i in 1:ncol(X_time)){
      colnames(b_list[[i]]) <- paste0("basis(knots = knots, x = log(times), ortho = ortho)",
                                      1:ncol(b_list[[i]]), ":",
                                      colnames(X_time)[i])
      b_list[[i]] <- b_list[[i]] * X_time[,i]
      db_list[[i]] <- db_list[[i]] * X_time[,i]
    }
    b_time <- do.call(cbind, b_list)
    db_time <- do.call(cbind, db_list)
  }else{
    b_time <- NULL
    db_time <- NULL
    #tvc.formula <- NULL
  }

  #Construct design matrix
  X <- model.matrix(smooth.formula, data = data)
  b <- cbind(b, X[,-1, drop = FALSE], b_time)
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
  likelihood.pars <- list(time = times,
                          event = event, X = X,
                          b = b, db = db, bhazard = data[, bhazard],
                          link_fun_pi = link_fun_pi,
                          link_fun_su = link_fun_su,
                          dlink_fun_su = dlink_fun_su)

  if(is.null(optim.args$control$maxit)){
    optim.args$control <- list(maxit = 10000)
  }


  #Generate initial values if these are not provided by the user
  if(is.null(optim.args$par)){
    if(verbose) cat("Finding initial values... ")
    inivalues <- lapply(ini.types,
                        function(ini.type) get.ini.values(smooth.formula = smooth.formula,
                                                          tvc.formula =  tvc.formula,
                                                          data = data,
                                                          bhazard = bhazard,
                                                          linkpi = linkpi,
                                                          linksu = linksu,
                                                          formula = formula,
                                                          type = type,
                                                          times = times,
                                                          event = event,
                                                          n.knots.time = n.knots.time,
                                                          knots.parse = knots,
                                                          X = X, b = b,
                                                          method = ini.type))
  }else{
    if(verbose) cat("Initial values provided by the user... ")
    inivalues <- optim.args$par
    optim.args <- optim.args[-which(names(optim.args) == "par")]
  }

  optim.pars <- c(optim.args, likelihood.pars)

  #Test if initial values are within the feasible region
  ini.eval <- sapply(inivalues, function(inival) do.call(minusloglik, c(likelihood.pars, list(inival))))
  run.these <- !is.na(ini.eval)

  if(verbose) cat("Completed!\nFitting the model... ")

  #Fit each model
  if(constr.optim){
    if(constr.optim & linksu != "loglog"){
      stop("Constrained optimization only works for linksu = 'loglog'")
    }
    optim.pars$ui <- cbind(matrix(0, nrow = nrow(X), ncol(X)), db)
    optim.pars$ci <- .Machine$double.eps
    optim.pars$f <- minusloglik
    optim.pars["grad"] <- list(NULL)

    res_list <- lapply(inivalues[run.these], function(inival){
      optim.pars$theta <- inival
      suppressWarnings(do.call(constrOptim, optim.pars))
    })

  }else{
    optim.pars$fn <- minusloglik
    res_list <- lapply(inivalues[run.these], function(inival){
      optim.pars$par <- inival
      suppressWarnings(do.call(optim, optim.pars))
    })
  }

  #Choose the best model according to the maximum likelihood estimate
  MLs <- sapply(res_list, function(x) tail(x$value, 1))
  wh <- which.min(MLs)
  res <- res_list[[wh]]

  #names(res$par) <- gsub("spline_", "", names(res$par))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(verbose) cat("Completed!\n")

  #Compute the covariance matrix matrix
  if(covariance){
    cov <- solve(numDeriv::hessian(minusloglik, res$par,
                                   time = times, event = event,
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
            NegMaxLik = min(MLs), covariance = cov,
            tvc.formula = tvc.formula, formula = formula,
            formula_main = smooth.formula, type = type,
            link_fun_pi = link_fun_pi,
            link_fun_su = link_fun_su,
            dlink_fun_su = dlink_fun_su,
            linkpi = linkpi, linksu = linksu,
            R.inv = R.inv, R.inv_list = R.inv_list, ortho = ortho,
            df = length(res$par) - 1, NegMaxLiks = MLs, optim.pars = optim.pars,
            times = times, constr.optim = constr.optim)

  class(L) <- c("fcm", "cuRe")
  L
}


#Function for computing initial values
get.ini.values <- function(smooth.formula, tvc.formula, data, bhazard, linkpi, linksu, formula, type,
                           times, event, n.knots.time, knots.parse, X, b, method = "cure"){
  if(method == "cure"){

    #Merge smooth.formula with variables in time-varying effects
    vars <- all.vars(tvc.formula)

    formula.2 <- as.formula(paste0(Reduce(paste, deparse(smooth.formula)),
                                   ifelse(length(vars), " + ", ""),
                                   vars,
                                   collapse = " + "))

    formula.new <- as.formula(paste0(deparse(lhs(formula)), " ~ -1 + X"))
    #Fit mixture or non-mixture cure model
    fit <- fit.cure.model(formula.new, data = data, bhazard = bhazard, covariance = F,
                          formula.k1 = formula.2, formula.k2 = ~ 1, type = type)

    #Scale by link function
    pi_hat <- get.link("logit")(X %*% fit$coefs[["gamma"]])
    gpi_hat <- get.inv.link(linkpi)(pi_hat)
    #Fit linear model to obtain initial values
    ini_pi <- lm(gpi_hat ~ -1 + X)$coefficients
    names(ini_pi) <- colnames(X)

    #Predict survival of the uncured
    lp <- exp(model.matrix(formula.2, data = data) %*% fit$coefs[[2]])
    shat <- exp(-lp * times ^ exp(fit$coefs[[3]]))
    gshat <- get.inv.link(linksu)(shat)

    finites <- is.finite(gshat)
    gshat <- gshat[finites]
    b <- b[finites,]
    fit_lm <- lm(gshat ~ -1 + b)

  }else if(method == "deaths"){
    formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                             response = formula[[2]])
    event2 <- 1 - event
    fit_glm <- glm(event2 ~ -1 + X, family = binomial(link = "logit"))
    pi_hat <- get.link("logit")(predict(fit_glm))
    gpi_hat <- get.inv.link(linkpi)(pi_hat)
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    fit <- coxph(formula.2, data = data[event == 1,])
    cum_base_haz <- get_basehaz(fit)
    shat <- exp(-cum_base_haz$hazard) ^ exp(fit$linear.predictors)
    suppressWarnings(gshat <- get.inv.link(linksu)(shat))
    fit_lm <- lm(gshat ~ -1 + b[event == 1,])
  }else if(method == "flexpara"){

    vars1 <- attr(terms(formula), "term.labels")
    vars2 <- c(attr(terms(smooth.formula), "term.labels"),
               attr(terms(tvc.formula), "term.labels"))

    wh.vars <- vars1[which(vars1 %in% vars2)]

    if(length(wh.vars)){
      rm.formula <- paste0(".~. -", paste0(wh.vars, collapse = " + "))
      formula.2 <- update(formula, rm.formula)
    }else{
      formula.2 <- formula
    }

    X_new <- lapply(list(formula.2, smooth.formula, tvc.formula), function(form){
      M <- model.matrix(form, data = data)[, -1, drop = FALSE]
      if(ncol(M)){
        wh <- grepl("basis\\(", colnames(M))
        if(any(wh)) colnames(M)[wh] <- paste0("b", 1:sum(wh))
      }
      M
    })

    X_new <- do.call(cbind, X_new)

    #Merge data into a single data frame
    data2 <- cbind(data, X_new)

    #Formula for survival of the uncred
    fuvar <- as.character(formula.2[[2]][[2]])
    base_formula <- as.formula(paste0("~basis(knots = knots.parse, x = log(",
                                      fuvar, "), ortho = F, intercept = F)"))

    rhs(formula.2) <- rhs(as.formula(paste0("~ ", paste0(c(1, colnames(X_new)), collapse = " + "))))
    environment(formula.2) <- environment()

    #Fit relative survival model
    suppressWarnings(fit <- rstpm2::stpm2(formula.2, data = data2,
                                          smooth.formula = base_formula,
                                          bhazard = data2[, bhazard]))

    #plot(fit, newdata = data.frame(age_years = 50, sex = c("male")), ylim = c(0,1))

    #Predict survival function
    shat <- predict(fit, newdata = data2, se.fit = F)

    #If predictions are all 1, we manually change these
    shat[shat == 1] <- shat[shat == 1] - 0.01
    #Scale by link function
    gshat <- get.inv.link(linksu)(shat)

    #Change follow-up times and predict cure rate
    fu_time <- all.vars(formula)[1]
    data2[, fu_time] <- max(data2[, fu_time]) + 0.1
    pi_hat <- predict(fit, newdata = data2, se.fit = F)

    #Change cases with increasing relative survival
    wh <- which(pi_hat >= shat)
    pi_hat[wh] <- shat[wh] - 0.01

    #Scale by link function
    data$gpi_hat <- get.inv.link(linkpi)(pi_hat)

    #Run linear model for pi to obtain initial values
    pi_fit <- lm(gpi_hat ~ -1 + X, data = data)
    ini_pi <- pi_fit$coefficients
    names(ini_pi) <- colnames(X)

    #Run linear model for S_u(t) to obtain initial values for either mixture or non-mixture models
    if(type == "mixture"){
      suhat <- (shat - pi_hat) / (1 - pi_hat)
    } else {
      suhat <- 1 - log(shat) / log(pi_hat)
    }
    gsuhat <- get.inv.link(linksu)(suhat)
    finites <- is.finite(gsuhat)
    gsuhat <- gsuhat[finites]
    b <- b[finites,]
    fit_lm <- lm(gsuhat ~ -1 + b)
  }

  #Make proper naming of the initial values
  suppressWarnings(coefs <- summary(fit_lm)$coefficients[, "Estimate"])
  names(coefs) <- colnames(b)
  ini_values <- c(ini_pi, coefs)
  #names(ini_values)[1:ncol(X)] <- colnames(X)
  ini_values
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

