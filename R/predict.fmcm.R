get.link <- function(x, type){
  if(type %in% c("curerate", "relsurv", "crudeprob", "survuncured")){
    exp(x) / (exp(x) + 1)
  }else if (type == "probcure"){
    pnorm(x)
  }else if(type == "ehaz"){
    x
  }else if(type == "iden"){
    x
  }
}

get.inv.link <- function(x, type){
  if(type %in% c("curerate", "relsurv", "crudeprob", "survuncured")){
    log(x / (1 - x))
  }else if(type == "probcure"){
    qnorm(x)
  }else if(type == "ehaz"){
    x
  }else if(type == "iden"){
    x
  }
}


relsurv_fun <- function(pars, M, M2, dM2, time, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars[1:ncol(M)], M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    get.inv.link(pi + (1 - pi) * link_fun_su(eta), type = "relsurv")
  }
  else if(model == "nmixture"){
    get.inv.link(pi ^ (1 - link_fun_su(eta)), type = "relsurv")
  }
}

ehaz_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars, M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  deta <- dM2 %*% pars[-c(1:ncol(M))]
  ds_u <- dlink_fun_su(eta)
  if(model == "mixture"){
    s_u <- link_fun_su(eta)
    get.inv.link(-(1 - pi) * ds_u * (deta / time) / (pi + (1 - pi) * s_u), type = "ehaz")
  }else if(model == "nmixture"){
    get.inv.link(log(pi) * ds_u * deta / time, type = "ehaz")
  }
}

probcure_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(get.link(pi_fun(pars, M), type = "curerate"))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    rsurv <- pi + (1 - pi) * link_fun_su(eta)
  }else if(model == "nmixture"){
    rsurv <- pi ^ (1 - link_fun_su(eta))
  }
  get.inv.link(pi / rsurv, type = "probcure")
}

survuncured_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    get.inv.link(link_fun_su(eta), type = "survuncured")
  }else if(model == "nmixture"){
    pi <- c(get.link(pi_fun(pars, M), type = "curerate"))
    rsurv <- pi ^ (1 - link_fun_su(eta))
    get.inv.link((rsurv - pi) / (1 - pi), type = "survuncured")
  }
}



#' Predict function for Flexible mixture cure model
#'
#' Function for doing predictions for class \code{fmcm}
#'
#' @param fit Object of class \code{fmcm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the the data which
#' the model was fitted on.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, \code{probcure} for the
#' conditional probability of being cured, and \code{survuncured} for the disease-specific survival of the uncured.
#' @param time Optional time points at which to compute predictions. This argument is not used if type is \code{curerate}.
#' @param ci Logical indicating whether confidence intervals should be computed
#' @param pars Numerical vector containing the parameters values of the model.
#' In general, this argument can be ignored by the user
#' @return A list containing the predictions of each individual in \code{newdata}.
#' @export

predict.fmcm <- function(fit, newdata = NULL, type = "relsurv",
                         time = NULL, ci = T, pars = NULL){
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
  pi <- pi_fun(all.coefs, M)

  if(type == "curerate"){
    pi <- data.frame(pi = pi)
    if(ci){
      grads <- jacobian(pi_fun, x = all.coefs, M = M)
      pi$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      pi$ci.lower <- get.link(pi$pi - qnorm(0.975) * sqrt(pi$var), type = "curerate")
      pi$ci.upper <- get.link(pi$pi + qnorm(0.975) * sqrt(pi$var), type = "curerate")
    }
    pi$pi <- get.link(pi$pi, type = "curerate")
    return(pi)
  }else{
    b <- basis(knots = fit$knots, x = log(time))
    db <- dbasis(knots = fit$knots, x = log(time))
    if(!is.null(fit$knots.time)){
      tvc.b <- lapply(fit$knots.time, basis, x = log(time))
      tvc.db <- lapply(fit$knots.time, dbasis, x = log(time))
    }
    M2.list <- lapply(1:nrow(newdata), function(i){
      model.matrix(fit$formula_main, newdata[i, ,drop = F])[rep(1, nrow(b)),-1, drop = F]
    })

    M_list <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      M2 <- cbind(b, M2.list[[i]])
      dM2 <- cbind(db, M2.list[[i]])
      if(!is.null(fit$tvc.formula)){
        M_time <- model.matrix(fit$tvc.formula, newdata)[,-1, drop = F]
        b_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.b[[j]] * M_time[i, j]))
        db_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.db[[j]] * M_time[i, j]))
        dM2 <- cbind(dM2, db_time)
        M2 <- cbind(M2, b_time)
      }
      M_list[[i]] <- list(M2 = M2, dM2 = dM2)
    }
    if(type %in% c("relsurv", "ehaz", "probcure", "survuncured")){
      if(type == "relsurv"){
        fun <- relsurv_fun
      }else if(type == "ehaz"){
        fun <- ehaz_fun
      }else if(type == "probcure"){
        fun <- probcure_fun
      }else if(type == "survuncured"){
        fun <- survuncured_fun
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
                            dM2 = M_list[[i]]$dM2, time = time, pi_fun = pi_fun, model = "mixture",
                            link_fun_pi = link_fun_pi, link_fun_su = link_fun_su, dlink_fun_su = dlink_fun_su)
          rss[[i]]$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
          rss[[i]]$ci.lower <- get.link(rss[[i]]$Est - qnorm(0.975) * sqrt(rss[[i]]$var), type = type)
          rss[[i]]$ci.upper <- get.link(rss[[i]]$Est + qnorm(0.975) * sqrt(rss[[i]]$var), type = type)
        }
        rss[[i]]$Est <- get.link(rss[[i]]$Est, type = type)

        if(type %in% c("relsurv", "survuncured")){
          if(ci){
            rss[[i]][time == 0, ] <- c(1, 0, 1, 1)
          }else{
            rss[[i]][time == 0, ] <- 1
          }
        }
      }
    }
    return(list(res = rss, time = time, type = type))
  }
}
