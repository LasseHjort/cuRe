fit_models <- function(data, bhazard, knots = NULL, n.knots = NULL, last_knot = NULL,
                       model = "relsurv", add.knots = NULL){
  if(is.null(knots)){
    ti <- data[data$status == 1, "FU_years"]
    if(is.null(last_knot)){
      last_knot <- max(ti)
    }
    knots <- c(min(ti),
               quantile(ti, 1 / (n.knots - 1)*1:(n.knots - 2)),
               last_knot)
    knots <- log(knots)
    knots <- sort(c(knots, add.knots))
  }
  if(model == "cure"){
    base_function <- basis_cure
    dbase_function <- dbasis_cure
  }else{
    base_function <- basis
    dbase_function <- dbasis
  }

  fit <- stpm2(Surv(FU_years, status == 1) ~ -1, data = data, bhazard = bhazard,
               smooth.formula = ~base_function(knots, log(FU_years)))
  L <- list(data = data, base_function = base_function, dbase_function = dbase_function,
            knots = knots, coefs = fit@coef, covariance = fit@vcov, ML = fit@min)
  class(L) <- "oldrelsurv"
  L
}

predict.oldrelsurv <- function(fit, newdata = NULL, time = NULL, pars = NULL, ci = F){
  if(!is.null(pars))
    fit$coefs <- pars
  b <- fit$base_function(knots = fit$knots, x = log(time))
  eta <- b %*% fit$coefs
  surv <- data.frame(Est = ifelse(time == 0, 1, exp(-exp(eta))))
  list(res = list(surv))
}




fit_models_time <- function(data, bhazard, knots = NULL, n.knots = NULL, last_knot = NULL,
                            model = "relsurv", add.knots = NULL, n.knots.time = NULL, knots.time = NULL){
  if(is.null(knots)){
    ti <- data[data$status == 1, "FU_years"]
    if(is.null(last_knot)){
      last_knot <- max(ti)
    }
    knots <- c(min(ti),
               quantile(ti, 1 / (n.knots - 1)*1:(n.knots - 2)),
               last_knot)
    knots <- log(knots)
    knots <- sort(c(knots, add.knots))
  }
  
  if(is.null(knots.time)){
    knots.time <- c(min(ti),
                    quantile(ti, 1 / (n.knots.time - 1)*1:(n.knots.time - 2)),
                    last_knot)
    knots.time <- log(knots.time)
  }
  
  if(model == "cure"){
    base_function <- basis_cure
    dbase_function <- dbasis_cure
  }else{
    base_function <- basis
    dbase_function <- dbasis
  }
  
  fit <- stpm2(Surv(FU_years, status == 1) ~ -1, data = data, bhazard = bhazard,
               smooth.formula = ~base_function(knots, log(FU_years)),
               tvc.formula = ~age_years:base_function(knots.time, log(FU_years)))
  L <- list(data = data, base_function = base_function, dbase_function = dbase_function,
            knots = knots, coefs = fit@coef, covariance = fit@vcov, ML = fit@min, knots.time = knots.time)
  class(L) <- "oldrelsurvtime"
  L
}

#fit <- fit_models_time(DLBCL, bhazard = DLBCL$exp_haz, n.knots = 6, model = "cure", n.knots.time = 3)
predict.oldrelsurvtime <- function(fit, newdata = NULL, time = NULL, pars = NULL, ci = F){
  if(!is.null(pars))
    fit$coefs <- pars
  b <- fit$base_function(knots = fit$knots, x = log(time))
  b.time <- fit$base_function(knots = fit$knots.time, x = log(time))
  b.time <- b.time * newdata$age_years
  b.all <- cbind(b, b.time)
  eta <- b.all %*% fit$coefs
  surv <- data.frame(Est = ifelse(time == 0, 1, exp(-exp(eta))))
  list(res = list(surv))
}



