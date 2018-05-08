data("colonDC")
colonDC <- colonDC[1:5000,]
colonDC$FU_years <- colonDC$FUyear
system.time(colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                                        data = colonDC, ratetable = survexp.dk))
bhazard <- "bhaz"
formula <- Surv(FU_years, status) ~ 1
data <- colonDC
smooth.formula <- ~ 1
knots <- NULL
n.knots <- 7
knots.time <- NULL
type <- "mixture"
linkpi <- c("logit")
linksu <- "loglog"
verbose <- T
optim.args <- NULL
ini.types <- c("cure", "flexpara")

FlexCureModel.bs <- function(formula, data, bhazard = NULL, smooth.formula = ~ 1,
                             knots = NULL, n.knots = 3,
                             knots.time = NULL, n.knots.time = NULL,
                             covariance = T, type = c("mixture", "nmixture"),
                             linkpi = c("logit", "identity", "loglog", "probit"),
                             linksu = c("loglog", "logit", "probit"),
                             verbose = T, constr.optim = F,
                             optim.args = NULL, ortho = TRUE,
                             ini.types = c("cure", "flexpara")){

  type <- match.arg(type)
  linkpi <- match.arg(linkpi)
  linksu <- match.arg(linksu)

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
  b <- bs(knots = knots[-c(1, length(knots))], x = log(times), Boundary.knots = range(knots))

  #Construct design matrix
  X <- model.matrix(smooth.formula, data = data)
  X.pi <- model.matrix(formula, data = data)


  #Extract link function
  link_fun_pi <- get.link(linkpi)
  link_fun_su <- get.link(linksu)
  dlink_fun_su <- get.dlink(linksu)

  #Extract minus log likelihood function
  minusloglik <- newmixturenegloglik

  #Create background hazard
  if(is.null(bhazard)){
    bhazard <- rep(0, nrow(data))
  }else {
    if(!is.numeric(bhazard)){
      bhazard <- data[, bhazard]
    }
  }

  if(length(bhazard) != nrow(data))
    stop("Length of bhazard is not the same as nrow(data)")

  #Prepare optimization arguments
  likelihood.pars <- list(time = times,
                          event = event, X = X[,-1, drop = FALSE],
                          b = b, X.pi = X.pi, bhazard = bhazard,
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

  fit <- fit.cure.model(formula, data = data, bhazard = bhazard)
  haz.eval <- exp(fit$coefs[[3]]) * exp(fit$coefs[[2]]) * times ^ (exp(fit$coefs[[3]]) - 1)
  log.haz <- log(haz.eval)
  basis <- bs(knots = knots[-c(1, length(knots))], x = log(times), Boundary.knots = range(knots))
  basis[, ncol(basis)] <- basis[, ncol(basis)] * 10
  new.fit <- lm(log.haz ~ -1 + basis[,-ncol(basis)] + offset(basis[,ncol(basis)]))

  pi_hat <- get.link("logit")(X.pi %*% fit$coefs[[1]])
  gpi_hat <- get.inv.link(linkpi)(pi_hat)
  #Fit linear model to obtain initial values
  ini_pi <- lm(gpi_hat ~ -1 + X.pi)$coefficients
  names(ini_pi) <- colnames(X)

  inivalues <- list(c(ini_pi, new.fit$coefficients))

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
                                   bhazard = bhazard,
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
                                          bhazard = bhazard))

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


#Prepare optimization arguments
likelihood.pars <- list(time = times,
                        event = event, X = X[,-1, drop = FALSE],
                        b = b, X.pi = X.pi, bhazard = bhazard,
                        link_fun_pi = link_fun_pi,
                        link_fun_su = link_fun_su,
                        dlink_fun_su = dlink_fun_su)

newmixturenegloglik <- function(param, event, times, X, b, X.pi, bhazard,
                                link_fun_pi, link_fun_su, dlink_fun_su){
  #Get parameters
  gamma <- param[1:ncol(X.pi)]
  beta <- param[(ncol(X.pi) + 1):length(param)]

  #Calculate linear predictors
  eta.pi <- X.pi %*% gamma
  pi <- link_fun_pi(eta.pi)
  hazu <- exp(b %*% c(beta, 10))

  #Calculate survival by Gaussian quadrature
  cumhaz <- rep(NA, length(times))
  for(i in 1:length(times)){

  }
  eta <- b %*% c(beta, 10)
  surv <- link_fun_su(eta)
  rsurv <- pi + (1 - pi) * surv
  likterms <- log(rsurv)

  #Add the hazard term only for events
  etaD <- exp(b %*% c(beta, 10))
  ehaz <- cure.type$haz(pi[event], link.surv$gradS(eta[event], etaD[event]), rsurv[event])
  haz <- haz.const <- bhazard[event] + ehaz
  haz[haz <= 0] <- .Machine$double.eps
  likterms[event] <- likterms[event] + log(haz)

  #Calculate hazard for which constraints are used
  if(constraint) haz.const <- link.surv$h(eta, etaD)

  #Output the negative log likelihood
  -sum( likterms ) + kappa / 2 * sum( ( haz.const[haz.const < 0] ) ^ 2 )
}
