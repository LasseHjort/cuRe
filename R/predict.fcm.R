
#' Predict function for Flexible mixture cure model
#'
#' Function for doing predictions for class \code{fmcm}
#'
#' @param fit Object of class \code{fcm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the data which
#' the model was fitted on.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, \code{probcure} for the
#' conditional probability of being cured, and \code{survuncured} for the disease-specific survival of the uncured.
#' @param time Optional time points at which to compute predictions.
#' This argument is not used if type is \code{curerate}.
#' @param ci Logical. If \code{TRUE}, confidence intervals are computed.
#' @param pars Numerical vector containing the parameters values of the model.
#' In general, this argument can be ignored by the user.
#' @return A list containing the predictions of each individual in \code{newdata}.
#' @export

predict.fcm <- function(fit, newdata = NULL, type = c("relsurv", "ehaz", "probcure", "survuncured"),
                        time = NULL, ci = T, pars = NULL){
  type <- match(type)
  if(!is.null(pars)){
    fit$coefs <- pars[1:length(fit$coefs)]
    if(length(fit$coefs) < length(pars)){
      fit$coefs.spline <- pars[(length(fit$coefs) + 1):length(pars)]
    }else{
      fit$coefs.spline <- NULL
    }
  }
  is_null_newdata <- is.null(newdata)
  if(is_null_newdata){
    tt <- terms(fit$formula)
    formula.2 <- formula(delete.response(tt))
    vars <- c(all.vars(formula.2), all.vars(fit$formula_main), all.vars(fit$tvc.formula))
    if(length(vars) != 0){
      stop("'newdata' must be specified for model including covariates")
    }
    newdata <- data.frame(x = 1)
    colnames(newdata) <- "(Intercept)"
  }
  link_fun_pi <- fit$link_fun_pi
  link_fun_su <- fit$link_fun_su
  dlink_fun_su <- fit$dlink_fun_su
  all.coefs <- c(fit$coefs, fit$coefs.spline)
  tt <- terms(fit$formula)
  formula.2 <- formula(delete.response(tt))
  M <- model.matrix(formula.2, newdata)
  pi_fun <- function(pars, M) M %*% pars[1:ncol(M)]
  pi <- get.inv.link("logit")(link_fun_pi(pi_fun(all.coefs, M)))

  if(type == "curerate"){
    pi <- data.frame(pi = pi)
    if(ci & fit$ci){
      grads <- jacobian(pi_fun, x = all.coefs, M = M)
      pi$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      pi$ci.lower <- get.link("logit")(pi$pi - qnorm(0.975) * sqrt(pi$var))
      pi$ci.upper <- get.link("logit")(pi$pi + qnorm(0.975) * sqrt(pi$var))
    }
    pi$pi <- get.link("logit")(pi$pi)
    return(pi)
  }else{
    if(is.null(time)){
      time <- seq(min(fit$times), max(fit$times), length.out = 100)
    }
    b <- basis(knots = fit$knots, x = log(time), ortho = fit$ortho, R.inv = fit$R.inv)
    db <- dbasis(knots = fit$knots, x = log(time), ortho = fit$ortho, R.inv = fit$R.inv)
    if(!is.null(fit$knots.time)){
      tvc.b <- lapply(1:length(fit$knots.time), function(i){
        basis(knots = fit$knots.time[[i]], x = log(time), ortho = fit$ortho, R.inv = fit$R.inv_list[[i]])
      })
      tvc.db <- lapply(1:length(fit$knots.time), function(i){
        dbasis(knots = fit$knots.time[[i]], x = log(time), ortho = fit$ortho, R.inv = fit$R.inv_list[[i]])
      })
    }
    M2.list <- lapply(1:nrow(newdata), function(i){
      model.matrix(fit$formula_main, newdata[i, ,drop = F])[rep(1, nrow(b)),-1, drop = F]
    })

    M_list <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      M2 <- cbind(b, M2.list[[i]])
      dM2 <- cbind(db, M2.list[[i]])
      if(!is.null(fit$knots.time)){
        M_time <- model.matrix(fit$tvc.formula, newdata)[,-1, drop = F]
        b_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.b[[j]] * M_time[i, j]))
        db_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.db[[j]] * M_time[i, j]))
        dM2 <- cbind(dM2, db_time)
        M2 <- cbind(M2, b_time)
      }
      M_list[[i]] <- list(M2 = M2, dM2 = dM2)
    }
    if(type == "relsurv"){
      fun <- relsurv_fun
      link.type <- "loglog"
    }else if(type == "ehaz"){
      fun <- ehaz_fun
      link.type <- "identity"
    }else if(type == "probcure"){
      fun <- probcure_fun
      link.type <- "probit"
    }else if(type == "survuncured"){
      fun <- survuncured_fun
      link.type <- "loglog"
    }

    rss <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      rss[[i]] <- data.frame(Est = c(fun(all.coefs, M2 = M_list[[i]]$M2, M = M[i,, drop = FALSE],
                                         dM2 = M_list[[i]]$dM2,
                                         time = time, pi_fun = pi_fun, model = fit$type,
                                         link_fun_pi = link_fun_pi, link_fun_su = link_fun_su,
                                         dlink_fun_su = dlink_fun_su)))
      if(ci){
        grads <- jacobian(fun, x = all.coefs, M2 = M_list[[i]]$M2, M = M[i, , drop = FALSE],
                          dM2 = M_list[[i]]$dM2, time = time, pi_fun = pi_fun, model = fit$type,
                          link_fun_pi = link_fun_pi, link_fun_su = link_fun_su, dlink_fun_su = dlink_fun_su)
        rss[[i]]$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
        ci1 <- get.link(link.type)(rss[[i]]$Est - qnorm(0.975) * sqrt(rss[[i]]$var))
        ci2 <- get.link(link.type)(rss[[i]]$Est + qnorm(0.975) * sqrt(rss[[i]]$var))
        wh <- round(length(ci1) / (2 - .Machine$double.eps))
        ci_max <- rep(which.max(c(ci1[wh], ci2[wh])), length(ci1))
        rss[[i]]$ci.lower <- ifelse(ci_max == 1, ci2, ci1)
        rss[[i]]$ci.upper <- ifelse(ci_max == 1, ci1, ci2)
      }
      rss[[i]]$Est <- get.link(link.type)(rss[[i]]$Est)

      if(type %in% c("relsurv", "survuncured")){
        if(ci){
          rss[[i]][time == 0, ] <- c(1, 0, 1, 1)
        }else{
          rss[[i]][time == 0, ] <- 1
          rss[[i]]$Est[is.na(rss[[i]]$Est)] <- 1
        }
      }
    }
    return(list(res = rss, time = time, type = type))
  }
}


relsurv_fun <- function(pars, M, M2, dM2, time, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars[1:ncol(M)], M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    get.inv.link("loglog")(pi + (1 - pi) * link_fun_su(eta))
  }
  else if(model == "nmixture"){
    get.inv.link("loglog")(pi ^ (1 - link_fun_su(eta)))
  }
}

ehaz_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars, M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  deta <- dM2 %*% pars[-c(1:ncol(M))]
  ds_u <- dlink_fun_su(eta)
  if(model == "mixture"){
    s_u <- link_fun_su(eta)
    get.inv.link("identity")(-(1 - pi) * ds_u * (deta / time) / (pi + (1 - pi) * s_u))
  }else if(model == "nmixture"){
    get.inv.link("identity")(log(pi) * ds_u * deta / time)
  }
}

probcure_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars, M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    rsurv <- pi + (1 - pi) * link_fun_su(eta)
  }else if(model == "nmixture"){
    rsurv <- pi ^ (1 - link_fun_su(eta))
  }
  get.inv.link("probit")(pi / rsurv)
}

survuncured_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    get.inv.link("loglog")(link_fun_su(eta))
  }else if(model == "nmixture"){
    pi <- c(link_fun_pi(pi_fun(pars, M)))
    rsurv <- pi ^ (1 - link_fun_su(eta))
    get.inv.link("loglog")((rsurv - pi) / (1 - pi))
  }
}
