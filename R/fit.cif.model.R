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
#' @import survival
#' @import rstpm2
#' @import numDeriv
#' @example inst/FlexCureModel.ex.R


#data <- read.csv("../../Kurser/AdvancedSurvivalAnalysis/Assignments/Data/bmt1715.txt", sep = " ")
#formula <- Surv(time, event) ~ 1
#n.knots <- list("1" = 3, "2" = 3)

fit.cif.model <- function(formula, data, bhazard,
                          knots = NULL, n.knots = NULL,
                          knots.time = NULL, n.knots.time = NULL,
                          covariance = T, link = "loglog",
                          verbose = T, constr.optim = F,
                          optim.args = NULL, ortho = TRUE,
                          ini.types = c("cure", "flexpara")){

  #Extract relevant variables
  times <- eval(formula[[2]][[2]], envir = data)
  event <- eval(formula[[2]][[3]], envir = data)
  d.times <- times[event == 1]
  cens <- 0
  events <- unique(sort(event[event != cens]))
  n.events <- length(events)

  if(is.null(n.knots) & is.null(knots)){
    n.knots <- rep(list(3), n.events)
    names(n.knots) <- events
  }

  #Caculate placement of knots and establish basis matrices
  if(is.null(knots)){
      knots <- vector("list", n.events)
      for(i in 1:n.events){
        d.times <- times[event == events[i]]
        bd_knots <- range(d.times)
        if(n.knots[[i]] > 2 ){
          inner_knots <- quantile(d.times, 1 / (n.knots[[i]] - 1)*1:(n.knots[[i]] - 2))
          knots[[i]] <- log(sort(c(bd_knots, inner_knots)))
        } else {
          knots[[i]] <- log(bd_knots)
        }
      }
  }

  #Evaluate baseline model matrices
  b <- lapply(knots, basis, x = log(times), ortho = ortho)
  #colnames(b) <- paste0("basis(knots = knots, x = log(times), ortho = ortho)", 1:ncol(b))
  R.inv <- attr(b, "R.inv")
  db <- lapply(knots, dbasis, x = log(times), ortho = ortho, R.inv = R.inv)


  #Construct design matrix
  X <- model.matrix(~ 1, data = data)
  for(i in 1:length(b)){
    b[[i]] <- cbind(b[[i]], X[,-1, drop = FALSE])
    db[[i]] <- cbind(db[[i]], matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)))
  }


  #Extract link function
  link_fun <- get.link(link)
  dlink_fun <- get.dlink(link)

  #Extract minus log likelihood function
  minusloglik <- minuslog_likelihood

  #Prepare optimization arguments
  likelihood.pars <- list(time = times,
                          event = event,
                          b = b, db = db,
                          link_fun = link_fun,
                          dlink_fun = dlink_fun)

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

  f <- function(t) 1 - exp(-exp(basis(x = log(t), knots = knots[[1]], ortho = F) %*% c(-1, 0.3, -0.0002)))
  curve(f, from = 0, to = 10, ylim = c(0, 1))

  inivalues <- rep(list(c(-1, 0.3, -0.0002)), 2)

  inivalues <- unlist(inivalues)
  optim.pars <- c(optim.args, likelihood.pars, list(inivalues))

  optim.pars$fn <- minuslog_likelihood
  debug(minuslog_likelihood)
  do.call(optim, optim.pars)

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
            NegMaxLik = min(MLs), covariance = cov, ci = covariance,
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
