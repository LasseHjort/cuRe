#' Fit spline-based cumulative incidence model
#'
#' The following function fits a generalized cumulative incidence function
#' using a link function formulation, i.e.,
#' \deqn{g(F_k(t|z)) \eta(t, z).
#'
#' @param formula Formula for modelling cause specific cumulative incidence function.
#' A list of formulas can be provided if different modelling is applied to the causes.
#' Reponse has to be of the form \code{Surv(time, status)}.
#' @param data Data frame in which to interpret the variables names in \code{formula}.
#' @param cens.code Numeric. The value of the event indicator which indicates censoring.
#' @param knots Knots used for the spline baseline effect of the cumulative incidence.
#' The knots are provided in a list with each entry denoting the knots corresponding to each cause.
#' @param n.knots Number of knots for the baseline splines.
#' The knots are calculated as the equidistant quantiles of the uncensored cause-specific event-times.
#' \code{n.knots} are provided in a list with each entry denoting the number of knots corresponding to each cause.
#' If \code{knots} is supplied, this argument will be ignored.
#' @param cure List of logicals. If \code{TRUE}, a cure model is fitted for the cause-specific cumulative incidence,
#' where cure is assumed after the last selected knot. Typically, cure is only assumed for one cumulative incidence function.
#' If \code{cure} is not a list, the same option is assumed for all causes.
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


fit.cif.model <- function(formula, data, cens.code = 0,
                          knots = NULL, n.knots = NULL, cure = FALSE,
                          knots.time = NULL, n.knots.time = NULL,
                          covariance = T, link = "loglog",
                          verbose = T, constr.optim = F,
                          optim.args = NULL, ortho = TRUE,
                          ini.types = c("cure", "flexpara")){

  if(!is.list(formula)) formula <- list(formula)

  #Extract relevant variables
  times <- eval(formula[[1]][[2]][[2]], envir = data)
  event <- eval(formula[[1]][[2]][[3]], envir = data)
  d.times <- times[event != cens.code]
  events <- unique(sort(event[event != cens.code]))
  n.events <- length(events)

  if(length(formula) == 1){
    formula <- rep(formula, n.events)
  }

  if(is.null(n.knots) & is.null(knots)){
    n.knots <- rep(list(3), n.events)
    names(n.knots) <- events
  }

  #Caculate placement of knots for all causes
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
    knots[[1]][length(knots[[1]])] <- knots[[1]][length(knots[[1]])] + 1e-09
    knots[[2]][length(knots[[2]])] <- knots[[2]][length(knots[[2]])] + 1e-09
  }

  #Caculate placement of knots in the time-varying covariate effects for all causes
  if(is.null(knots.time) & !is.null(n.knots.time)){
    knots.time <- vector("list", length(n.knots.time))
    for(i in 1:n.events){
      d.times <- times[event == events[i]]
      bd_knots <- range(d.times)
      knots.time[[i]] <- lapply(n.knots.time[[i]], function(nk){
        if(nk > 2){
          inner_knots <- quantile(d.times, 1 / (nk - 1)*1:(nk - 2))
          log(sort(c(bd_knots, inner_knots)))
        } else {
          log(bd_knots)
        }
      })
    }
  }

  if(!is.null(knots.time)){
    tvc.formula <- vector("list", n.events)
    for(i in 1:n.events){
      formula.string <- paste("~", paste(names(knots.time[[i]]), collapse = "+"))
      tvc.formula[[i]] <- as.formula(formula.string)
    }
  }else{
    tvc.formula <- NULL
  }


  #Extract basis functions and dbasis functions based on the cure indicator
  if(length(cure) == 1){
    cure <- rep(list(cure), n.events)
  }

  basis.func <- lapply(cure, function(cure.ind){
    if(cure.ind) basis.cure
    else basis
  })

  dbasis.func <- lapply(cure, function(cure.ind){
    if(cure.ind) dbasis.cure
    else dbasis
  })


  #Evaluate baseline model matrices
  b <- db <- R.inv <- vector("list", n.events)
  for(i in 1:n.events){
    b[[i]] <- basis.func[[i]](knots = knots[[i]], x = log(times), ortho = ortho)
    colnames(b[[i]]) <- paste0("b", 1:ncol(b[[i]]), "_c", i)
    if(ortho) R.inv[[i]] <- attr(b[[i]], "R.inv")
    db[[i]] <- dbasis.func[[i]](knots = knots[[i]], x = log(times), ortho = ortho, R.inv = R.inv[[i]])
  }

  #Construct design matrix for fixed effects and baseline splines
  for(i in 1:n.events){
    formula.rhs <- as.formula(paste0("~", rhs(formula[[i]])))
    X <- model.matrix(formula.rhs, data = data)[,-1, drop = FALSE]
    if(ncol(X)) colnames(X) <- paste0(colnames(X), "_c", i)
    b[[i]] <- cbind(b[[i]], X)
    db[[i]] <- cbind(db[[i]], matrix(0, ncol = ncol(X), nrow = nrow(X)))
  }


  #Construct desing matrix for time-varying effects
  if(!is.null(tvc.formula)){
    b.time <- db.time <- R.inv.time <- vector("list", n.events)
    for(i in 1:n.events){
      M <- model.matrix(tvc.formula[[i]], data = data)[, -1, drop = FALSE]
      b.list <- db.list <- R.inv.list <- vector("list", ncol(M))
      for(j in 1:length(knots.time[[i]])){
        b.tmp <- basis.func[[i]](knots = knots.time[[i]][[j]], x = log(times), ortho = ortho)
        b.list[[j]] <- b.tmp * M[,j]
        colnames(b.list[[j]]) <- paste0(colnames(M)[j],":b", 1:ncol(b.list[[j]]), "_c", i)
        if(ortho) R.inv.list[[j]] <- attr(b.tmp, "R.inv")
        db.tmp <- dbasis.func[[i]](knots = knots.time[[i]][[j]], x = log(times),
                                   ortho = ortho, R.inv = R.inv.list[[j]])
        db.list[[j]] <- db.tmp * M[,j]

      }
      b.time[[i]] <- do.call(cbind, b.list)
      R.inv.time[[i]] <- R.inv.list
      db.time[[i]] <- do.call(cbind, db.list)
      b[[i]] <- cbind(b[[i]], b.time[[i]])
      db[[i]] <- cbind(db[[i]], db.time[[i]])
    }

  }


  #Extract link function
  link_fun <- get.link(link)
  dlink_fun <- get.dlink(link)

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

    fit_cuminc <- cmprsk::cuminc(times, event)
    pred_cuminc <- t(cmprsk::timepoints(fit_cuminc, times = times)[[1]])
    pred_cuminc[pred_cuminc == 0] <- 0.0001
    pred_cuminc <- pred_cuminc[as.character(times),]

    pred_cuminc_trans <- get.inv.link(link)(1 - pred_cuminc)

    inivalues <- c()
    for(i in 1:length(b)){
      #finites <- is.finite(pred_cuminc_trans[, i])

      deaths <- event == i
      Dmat <- t(b[[i]][deaths,]) %*% b[[i]][deaths,]
      dvec <- t(pred_cuminc_trans[deaths, i]) %*% b[[i]][deaths,]
      eps <- rep(1e-09, length(which(deaths)))
      Amat <- t(db[[i]][deaths,])
      fit.qp <- solve.QP(Dmat = Dmat, dvec = dvec,
                         Amat = Amat, bvec = eps)

      #fit.qp$solution
      #fit.lm$coefficients

      #fit.lm <- lm(pred_cuminc_trans[finites, i] ~ -1 + b[[i]][finites,])
      #names(fit.lm$coefficients) <- colnames(b[[i]])
      names(fit.qp$solution) <- colnames(b[[i]])
      #inivalues <- c(inivalues, fit.lm$coefficients)
      inivalues <- c(inivalues, fit.qp$solution)
    }
    optim.args <- c(optim.args, list(inivalues))
  }else{
    if(verbose) cat("Initial values provided by the user... ")
  }

  #Combine parameters for optimization
  optim.pars <- c(optim.args, likelihood.pars)

  #Extract minus log likelihood function
  optim.pars$fn <- minuslog_likelihood
  #optim.pars$hessian <- TRUE

  if(verbose) cat("Completed!\nFitting the model... ")

  #Run optimization
  if(constr.optim){
    if(constr.optim & linksu != "loglog"){
      stop("Constrained optimization only works for linksu = 'loglog'")
    }
    optim.pars$ui <- cbind(matrix(0, nrow = nrow(X), ncol(X)), db)
    optim.pars$ci <- .Machine$double.eps
    optim.pars$f <- minuslog_likelihood
    optim.pars["grad"] <- list(NULL)
    res <- suppressWarnings(do.call(constrOptim, optim.pars))
  }else{
    res <- suppressWarnings(do.call(optim, optim.pars))
  }

  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(verbose) cat("Completed!\n")

  #Compute the covariance matrix matrix
  if(covariance){
    cov <- solve(pracma::hessian(minuslog_likelihood, res$par,
                                 time = times, event = event,
                                 b = b, db = db, link_fun = link_fun,
                                 dlink_fun = dlink_fun))
    #cov <- solve(res$hessian)
  }else{
    cov <- NULL
  }

  #Output the results
  L <- list(formula = formula, data = data,
            coefs = res$par, basis.func = basis.func, dbasis.func = dbasis.func,
            knots = knots, knots.time = knots.time,
            NegMaxLik = res$value, covariance = cov, ci = covariance,
            tvc.formula = tvc.formula, formula = formula,
            basis.func = basis.func, dbasis.func = dbasis.func,
            link_fun = link_fun, dlink_fun = dlink_fun,
            link = link, R.inv = R.inv, R.inv.time = R.inv.time,
            ortho = ortho, df = length(res$par) - 1, optim.pars = optim.pars,
            times = times, constr.optim = constr.optim)

  class(L) <- c("fcim", "cuRe")
  L
}

