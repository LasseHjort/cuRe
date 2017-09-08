
#' Predict function for cure models
#'
#' This function is used to make predictions of the cure models.
#'
#' @param fit An object of type CureModel.
#' @param newdata A dataframe with new data to predict from.
#' @param times A numeric vector including the time points for which to make predictions.
#' @param type A character denoting the type of prediction. Possible values are \code{relsurv} (default), \code{curerate}, and \code{probcure} (see details)
#' @param ci Logical indicating whether to include confidence intervals. Default is \code{TRUE}
#' @return An object of class \code{matrix} including the predictions.
#' @details For type \code{relsurv}, relative survival predictions are made, while for \code{curerate}, \eqn{\pi} and for \code{probcure}\cr
#' the probability of cure is calculated.
#' @export
#'
predict.curemodel <- function(fit, newdata = NULL, type = "relsurv",
                              time = NULL, ci = T, pars = NULL){
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
  link_fun_pi <- fit$link_fun_pi
  link_fun_su <- fit$link_fun_su
  dlink_fun_su <- fit$dlink_fun_su

  Ms <- lapply(formulas, get_design, data = newdata)
  pi_fun <- function(pars, M) M %*% pars[1:ncol(M)]
  pi <- pi_fun(unlist(fit$coefs), Ms[[1]])

  if(type == "curerate"){
    pi <- data.frame(pi = pi)
    if(ci){
      grads <- jacobian(pi_fun, x = unlist(fit$coefs), M = Ms[[1]])
      pi$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
      pi$ci.lower <- get.link(pi$pi - qnorm(0.975) * sqrt(pi$var), type = "curerate")
      pi$ci.upper <- get.link(pi$pi + qnorm(0.975) * sqrt(pi$var), type = "curerate")
    }
    pi$pi <- get.link(pi$pi, type = "curerate")
    return(pi)
  }else if(type %in% c("relsurv", "ehaz", "probcure")){
    if(is.null(time)){
      obs.times <- eval(fit$formulas[[1]][[2]], envir = fit$data)
      time <- seq(min(obs.times), max(obs.times), length.out = 100)
    }

    out.fun <- switch(type,
                      relsurv = relsurv_fun_simple,
                      ehaz = ehaz_fun_simple,
                      probcure = probcure_fun_simple)

    rss <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      Ms_indi <- lapply(Ms, function(x) x[i,, drop = F])
      rss[[i]] <- data.frame(Est = c(out.fun(Ms_indi = Ms_indi, pars = unlist(fit$coefs), time = time,
                                             dist = fit$dist, model = fit$type, surv_fun = fit$surv_fun,
                                             dens_fun = fit$dens_fun)))
      if(ci){
        grads <- jacobian(out.fun, x = unlist(fit$coefs), Ms_indi = Ms_indi, time = time,
                          dist = fit$dist, model = fit$type, surv_fun = fit$surv_fun,
                          dens_fun = fit$dens_fun)

        rss[[i]]$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$cov %*% x)
        rss[[i]]$ci.lower <- get.link(rss[[i]]$Est - qnorm(0.975) * sqrt(rss[[i]]$var), type = type)
        rss[[i]]$ci.upper <- get.link(rss[[i]]$Est + qnorm(0.975) * sqrt(rss[[i]]$var), type = type)
      }
      rss[[i]]$Est <- get.link(rss[[i]]$Est, type = type)

      if(type == "relsurv"){
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
  surv_fun <- get_surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(lps[[1]], type = "curerate")
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link(pi.eval + (1 - pi.eval) * surv.eval, type = type)
  }else if(model == "nmixture"){
    get.inv.link(pi.eval ^ (1 - surv.eval), type = type)
  }
}

ehaz_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "ehaz", surv_fun, dens_fun){
  surv_fun <- get_surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(lps[[1]], type = "curerate")
  dens.eval <- dens_fun(time, lps = lps)
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link((1 - pi.eval) * dens.eval / (pi.eval + (1 - pi.eval) * surv.eval), type = type)
  }else if(model == "nmixture"){
    get.inv.link(-log(pi.eval) * dens.eval, type = type)
  }
}

probcure_fun_simple <- function(Ms_indi, pars, time, dist, model, link, type = "probcure", surv_fun, dens_fun){
  surv_fun <- get_surv(dist)
  lps <- calc.lps(Xs = Ms_indi, param = pars)
  lps <- lapply(lps, c)
  pi.eval <- get.link(lps[[1]], type = "curerate")
  surv.eval <- surv_fun(time, lps = lps)
  if(model == "mixture"){
    get.inv.link(pi.eval / (pi.eval + (1 - pi.eval) * surv.eval), type = type)
  }else if(model == "nmixture"){
    get.inv.link(pi.eval / (pi.eval ^ (1 - surv.eval)), type = type)
  }
}


