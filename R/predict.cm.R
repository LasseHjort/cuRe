
#' Predict function for cure models
#'
#' This function is used to make predictions of the cure models.
#'
#' @param fit Object of class \code{cm} to do predictions from.
#' @param newdata Data frame from which to compute predictions. If empty, predictions are made on the data which
#' the model was fitted on.
#' @param time Optional time points at which to compute predictions.
#' @param type Type of prediction to do. Possible values are \code{relsurv} (default) for the relative survival,
#' \code{curerate} for the cure rate, \code{ehaz} for the excess hazard, \code{probcure} for the
#' conditional probability of being cured, and \code{survuncured} for the disease-specific survival of the uncured.
#' @param ci Logical. If \code{TRUE}, confidence intervals are computed.
#' @param pars Numerical vector containing the parameters values of the model.
#' @return An object of class \code{matrix} including the predictions.
#' @export
#'
predict.cm <- function(fit, newdata = NULL, type = c("relsurv", "curerate", "ehaz", "probcure", "survuncured"),
                              time = NULL, ci = T, pars = NULL){
  type <- match.arg(type)
  if(!is.null(pars)){
    groups <- factor(rep(1:4, fit$n.param.formula), 1:4, labels = c("gamma", "k1", "k2", "k3"))
    fit$coefs <- split(pars, f = groups)
  }
  is_null_newdata <- is.null(newdata)

  #Edit formulas
  formulas <- fit$formulas
  tt <- terms(formulas[[1]])
  formulas[[1]] <- formula(delete.response(tt))

  #Check if covariates are included in the model in cases where newdata is not provided
  if(is_null_newdata){
    vars <- lapply(formulas, function(formula){
      if(!is.null(formula)){
        all.vars(formula)
      }else{
        return(NULL)
      }
    })
    vars <- unlist(vars)
    if(length(vars) != 0){
      stop("'newdata' must be specified for model including covariates")
    }
    newdata <- data.frame(x = 1)
    colnames(newdata) <- "(Intercept)"
  }

  Ms <- lapply(formulas, get_design, data = newdata)
  pi_fun <- function(pars, M) M %*% pars[1:ncol(M)]
  pi <- pi_fun(unlist(fit$coefs), Ms[[1]])

  if(type == "curerate"){
    pi <- data.frame(pi = pi)
    if(ci & fit$ci){
      grads <- jacobian(pi_fun, x = unlist(fit$coefs), M = Ms[[1]])
      pi$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
      pi$ci.lower <- get.link(fit$link)(pi$pi - qnorm(0.975) * sqrt(pi$var))
      pi$ci.upper <- get.link(fit$link)(pi$pi + qnorm(0.975) * sqrt(pi$var))
    }
    pi$pi <- get.link(fit$link)(pi$pi)
    return(pi)
  }else if(type %in% c("relsurv", "ehaz", "probcure", "survuncured")){
    if(is.null(time)){
      time <- seq(min(fit$times), max(fit$times), length.out = 100)
    }

    out.fun <- switch(type,
                      relsurv = relsurv_fun_simple,
                      ehaz = ehaz_fun_simple,
                      probcure = probcure_fun_simple,
                      survuncured = survuncured_fun_simple)

    link.type <- switch(type,
                        relsurv = "logit",
                        ehaz = "identity",
                        probcure = "probit",
                        survuncured = "logit")

    rss <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      Ms_indi <- lapply(Ms, function(x) x[i,, drop = F])
      rss[[i]] <- data.frame(Est = c(out.fun(Ms_indi = Ms_indi, pars = unlist(fit$coefs), time = time,
                                             dist = fit$dist, model = fit$type, link = fit$link,
                                             surv_fun = fit$surv_fun,
                                             dens_fun = fit$dens_fun)))
      if(ci){
        grads <- jacobian(out.fun, x = unlist(fit$coefs), Ms_indi = Ms_indi, time = time,
                          dist = fit$dist, model = fit$type, link = fit$link,
                          surv_fun = fit$surv_fun,
                          dens_fun = fit$dens_fun)

        rss[[i]]$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
        rss[[i]]$ci.lower <- get.link(link.type)(rss[[i]]$Est - qnorm(0.975) * sqrt(rss[[i]]$var))
        rss[[i]]$ci.upper <- get.link(link.type)(rss[[i]]$Est + qnorm(0.975) * sqrt(rss[[i]]$var))
      }
      rss[[i]]$Est <- get.link(link.type)(rss[[i]]$Est)

      if(type %in% c("relsurv", "survuncured")){
        if(ci){
          rss[[i]][time == 0, ] <- c(1, 0, 1, 1)
        }else{
          rss[[i]][time == 0, ] <- 1
        }
      }
    }
    return(list(res = rss, time = time, type = type))
  }
}


relsurv_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "relsurv", surv_fun, dens_fun){
  surv_fun <- get.surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(link)(lps[[1]])
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link("logit")(pi.eval + (1 - pi.eval) * surv.eval)
  }else if(model == "nmixture"){
    get.inv.link("logit")(pi.eval ^ (1 - surv.eval))
  }
}

ehaz_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "ehaz", surv_fun, dens_fun){
  surv_fun <- get.surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(link)(lps[[1]])
  dens.eval <- dens_fun(time, lps = lps)
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link("identity")((1 - pi.eval) * dens.eval / (pi.eval + (1 - pi.eval) * surv.eval))
  }else if(model == "nmixture"){
    get.inv.link("identity")(-log(pi.eval) * dens.eval)
  }
}

probcure_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "probcure", surv_fun, dens_fun){
  surv_fun <- get.surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(link)(lps[[1]])
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link("probit")(pi.eval / (pi.eval + (1 - pi.eval) * surv.eval))
  }else if(model == "nmixture"){
    get.inv.link("probit")(pi.eval / (pi.eval ^ (1 - surv.eval)))
  }
}

survuncured_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "survuncured", surv_fun, dens_fun){
  surv_fun <- get.surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link("logit")(surv.eval)
  }else if(model == "nmixture"){
    pi.eval <- get.link(link)(lps[[1]])
    get.inv.link("logit")((pi.eval ^ (1 - surv.eval) - pi.eval) / (1 - pi.eval))
  }
}


