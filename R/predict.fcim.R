#' Predict function for Flexible parametric cumulative incidence function
#'
#' Function for doing predictions for class \code{fcim}
#'
#' @param fit Object of class \code{fcim} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the data which
#' the model was fitted on.
#' @param type Type of prediction to do. Possible values are \code{cuminc} (default) for the cumulative incidence,
#' and \code{subhaz} for the sub-distribution hazard.
#' @param time Optional time points at which to compute predictions.
#' This argument is not used if type is \code{curerate}.
#' @param ci Logical. If \code{TRUE}, confidence intervals are computed.
#' @param pars Numerical vector containing the parameters values of the model.
#' In general, this argument can be ignored by the user.
#' @return A list containing the predictions of each individual in \code{newdata}.
#' @export

predict.fcim <- function(fit, newdata = NULL, type = "cuminc",
                         time = NULL, ci = T, pars = NULL){
  if(!is.null(pars)){
    fit$coefs <- pars
  }
  is_null_newdata <- is.null(newdata)
  if(is_null_newdata){
    vars <- all.vars(rhs(fit$formula[[1]]))
    if(length(vars) != 0){
      stop("'newdata' must be specified for model including covariates")
    }
    newdata <- data.frame(x = 1)
    colnames(newdata) <- "(Intercept)"
  }
  link_fun <- fit$link_fun
  dlink_fun <- fit$dlink_fun
  if(is.null(time)){
    time <- seq(min(fit$times), max(fit$times), length.out = 100)
  }

  n.causes <- length(fit$knots)
  b <- db <- vector("list", n.causes)
  for(i in 1:n.causes){
    pars <- list(knots = fit$knots[[i]], x = log(time), ortho = fit$ortho, R.inv = fit$R.inv[[i]])
    b[[i]] <- do.call(fit$basis.func[[i]], pars)
    db[[i]] <- do.call(fit$dbasis.func[[i]], pars)
    formula.rhs <- as.formula(paste0("~", rhs(fit$formula[[i]])))
  }

  if(!is.null(fit$knots.time)){
    b.time <- db.time <- vector("list", n.causes)
    for(i in 1:n.causes){
      b.time[[i]] <- lapply(1:length(fit$knots.time[[i]]), function(j){
        fit$basis.func[[i]](knots = fit$knots.time[[i]][[j]], x = log(time),
                            ortho = fit$ortho, R.inv = fit$R.inv.time[[i]][[j]])
      })
      db.time[[i]] <- lapply(1:length(fit$knots.time), function(i){
        fit$dbasis.func[[i]](knots = fit$knots.time[[i]][[j]], x = log(time),
                             ortho = fit$ortho, R.inv = fit$R.inv.time[[i]][[j]])
      })
    }
  }


  X.list <- dX.list <- vector("list", nrow(newdata))
  for(j in 1:nrow(newdata)){
    X.cause <- dX.cause <- vector("list", n.causes)
    for(i in 1:n.causes){
      M <- model.matrix(formula.rhs, data = newdata[j,, drop = F])[rep(1, nrow(b[[i]])), -1, drop = F]
      X.cause[[i]] <- cbind(b[[i]], M)
      dX.cause[[i]] <- cbind(db[[i]], matrix(0, ncol = ncol(M), nrow = nrow(M)))
      if(!is.null(fit$knots.time)){
        M_time <- model.matrix(fit$tvc.formula[[i]], newdata[j,,drop = F])[,-1, drop = F]
        b_time <- do.call(cbind, lapply(1:ncol(M_time), function(k) b.time[[i]][[k]] * M_time[,k]))
        db_time <- do.call(cbind, lapply(1:ncol(M_time), function(k) db.time[[j]][[k]] * M_time[, k]))
        X.cause[[i]] <- cbind(X.cause[[i]], b_time)
        dX.cause[[i]] <- cbind(dX.cause[[i]], db_time)
      }
    }
    X.list[[j]] <- X.cause
    dX.list[[j]] <- dX.cause
  }


  if(type == "cuminc"){
    fun <- cuminc_fun
    link.type <- "loglog"
  }else if(type == "subhaz"){
    fun <- subhaz_fun
    link.type <- "probit"
  }

  ncols <- sapply(X.list[[1]], ncol)
  par.ind <- rep(1:n.causes, ncols)
  par.ind <- split(1:length(fit$coefs), f = par.ind)

  ests <- vector("list", nrow(newdata))
  for(i in 1:nrow(newdata)){
    ests[[i]] <- vector("list", n.causes)
    for(j in 1:n.causes){
      ests[[i]][[j]] <- data.frame(Est = c(fun(fit$coefs, M = X.list[[i]][[j]], dM = dX.list[[i]][[j]],
                                               time = time, link_fun = link_fun, dlink_fun = dlink_fun,
                                               cause = j, par.ind = par.ind)))
      if(ci){
        grads <- jacobian(fun, x = fit$coefs, M = X.list[[i]][[j]], dM = dX.list[[i]][[j]],
                          time = time, link_fun = link_fun, dlink_fun = dlink_fun,
                          cause = j, par.ind = par.ind)
        ests[[i]][[j]]$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
        ci1 <- get.link(link.type)(ests[[i]][[j]]$Est - qnorm(0.975) * sqrt(ests[[i]][[j]]$var))
        ci2 <- get.link(link.type)(ests[[i]][[j]]$Est + qnorm(0.975) * sqrt(ests[[i]][[j]]$var))
        ests[[i]][[j]]$ci.lower <- pmin(ci1, ci2)
        ests[[i]][[j]]$ci.upper <- pmax(ci1, ci2)
      }
      ests[[i]][[j]]$Est <- get.link(link.type)(ests[[i]][[j]]$Est)

      if(type %in% c("cuminc")){
        if(ci){
          ests[[i]][[j]][time == 0, ] <- rep(0, 4)
        }else{
          ests[[i]][[j]][time == 0, ] <- 0
          ests[[i]][[j]]$Est[is.na(ests[[i]][[j]]$Est)] <- 0
        }
      }
    }
  }
  return(list(res = ests, time = time, type = type))
}


cuminc_fun <- function(pars, M, dM, time, link_fun, dlink_fun, cause, par.ind){
  eta <- M %*% pars[par.ind[[cause]]]
  cuminc <- 1 - link_fun(eta)
  get.inv.link("loglog")(cuminc)
}

subhaz_fun <- function(pars, M, dM, time, link_fun, dlink_fun, cause, par.ind){
  eta <- M %*% pars[par.ind[[cause]]]
  deta <- dM %*% pars[par.ind[[cause]]]
  surv <- link_fun(eta)
  dsurv <- dlink_fun(eta)
  get.inv.link("probit")(- dsurv / surv * deta / time)
}
