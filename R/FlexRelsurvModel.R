FlexRelsurvModel <- function(formula, data, bhazard,
                             knots = NULL, n.knots = NULL,
                             tvc.formula = NULL, knots.time = NULL, n.knots.time = NULL,
                             hes = T, message = T, cure = F){

  #Extract relevant variables
  fu <- eval(formula[[2]][[2]], envir = data)
  status <- eval(formula[[2]][[3]], envir = data)
  death_times <- fu[status == 1]

  if(cure){
    base_function <- basis_cure
    dbase_function <- dbasis_cure
  }else{
    base_function <- basis
    dbase_function <- dbasis
  }

  #Caculate placement of knots and establish basis matrices
  if(is.null(knots)){
    bd_knots <- range(death_times)
    inner_knots <- quantile(death_times, 1 / (n.knots - 1)*1:(n.knots - 2))
    if(cure){
      inner_knots <- c(inner_knots, quantile(death_times, 0.95))
    }
    knots <- c(bd_knots, inner_knots)
    knots <- log(sort(knots))
  }

  b <- base_function(knots = knots, log(fu))
  db <- dbase_function(knots = knots, log(fu))

  if(!is.null(tvc.formula)){
    if(is.null(knots.time)){
      vars <- all.vars(tvc.formula)
      knots.time <- lapply(vars, function(x){
        bd_knots <- range(death_times)
        if(n.knots.time[[x]] > 2){
          inner_knots <- quantile(death_times, 1 / (n.knots.time[[x]] - 1)*1:(n.knots.time[[x]] - 2))
        }else{
          inner_knots <- NULL
        }
        if(cure){
          inner_knots <- c(inner_knots, quantile(death_times, 0.95))
        }
        log(sort(c(bd_knots, inner_knots)))
      })
      b_list <- lapply(knots.time, base_function, x = log(fu))
      db_list <- lapply(knots.time, dbase_function, x = log(fu))
    }
  }else{
    tvc.formula <- ~ 1
  }

  X_time <- model.matrix(tvc.formula, data = data)[,-1, drop = FALSE]
  if(ncol(X_time) > 0){
    for(i in 1:ncol(X_time)){
      b_list[[i]] <- b_list[[i]] * X_time[,i]
      db_list[[i]] <- db_list[[i]] * X_time[,i]
    }
    b_time <- do.call(cbind, b_list)
    db_time <- do.call(cbind, db_list)
  }else{
    tvc.formula <- NULL
    b_time <- NULL
    db_time <- NULL
  }


  #Construct design matrix
  X <- model.matrix(formula, data = data)
  b <- cbind(b, X[, -1], b_time)
  db <- cbind(db, matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)), db_time)

  #Generate initial values by running a mixture weibull cure model
  if(message) cat("Finding initial values...")
  #formula(as.Formula(terms(formula),formula(Formula(terms(tvc.formula)), lhs=0)), collapse=TRUE)
  if(is.null(tvc.formula)){
    formula.new <- formula
  }else{
    all_terms <- c(as.character(terms(formula)[[3]]), as.character(terms(tvc.formula)[[2]]))
    formula.new <- reformulate(termlabels = all_terms, response = terms(formula)[[2]])
  }
  fit <- coxph(formula.new, data = data)
  #fit <- coxph(Surv(FU_years, status) ~ age_years + sex, data = data[status == 1,])
  cum_base_haz <- get_basehaz(fit)

  cumhazs <- list(cum_base_haz$hazard[status == 1],
                  cum_base_haz$hazard[status == 1] - 0.1 * data[status == 1, bhazard] * fu[status == 1])

  logHs <- suppressWarnings(lapply(cumhazs, function(x) log(x * exp(fit$linear.predictors[status == 1]))))
  lm_fits <- lapply(logHs, function(logH) lm(logH ~ -1 + b[status == 1,]))

  #Make proper naming of the initial values
  ini_values <- lapply(lm_fits, function(lm_fit){
    suppressWarnings(coefs <- lm_fit$coefficients)
    nr.active.knots <- ifelse(cure, length(knots) - 1, length(knots))
    names(coefs)[1:nr.active.knots] <- paste0("gamma", 1:nr.active.knots)
    fix_var <- all.vars(formula(delete.response(terms(formula))))
    if(length(fix_var) != 0){
      if(length(coefs) > nr.active.knots){
        names(coefs)[(nr.active.knots + 1):(nr.active.knots + length(fix_var))] <- paste0("spline_", fix_var)
      }
      if(length(coefs) > (nr.active.knots + length(fix_var))){
        nr.active.knots.time <- lapply(knots.time, function(x) length(x) - ifelse(cure, 1, 0))
        var.names <- vector("list", length(nr.active.knots.time))
        for(i in 1:length(nr.active.knots.time)){
          var.names[[i]] <- paste0("gamma", 1:nr.active.knots.time[[i]], ":", names(n.knots.time)[i])
        }
        names(coefs)[(nr.active.knots + length(fix_var) + 1):length(coefs)] <- unlist(var.names)
      }
    }else{
      if(length(coefs) > nr.active.knots){
        nr.active.knots.time <- lapply(knots.time, function(x) length(x) - ifelse(cure, 1, 0))
        var.names <- vector("list", length(nr.active.knots.time))
        for(i in 1:length(nr.active.knots.time)){
          var.names[[i]] <- paste0("gamma", 1:nr.active.knots.time[[i]], ":", names(n.knots.time)[i])
        }
        names(coefs)[(nr.active.knots + 1):length(coefs)] <- unlist(var.names)
      }
    }
    coefs
  })

  if(all(ini_values[[1]] == ini_values[[2]])){
    ini_values <- ini_values[1]
  }

  if(message) cat("Completed!\nFitting the model...")

  #Fit the model

  #fit <- stpm2(Surv(FU_years, status == 1) ~ -1, data = data, bhazard = bhazard,
  #             smooth.formula = ~basis(knots, log(FU_years)))
  res_list <- lapply(ini_values, function(ini_val) optim(ini_val, fn = flexible_minuslog_likelihood,
                                                         time = fu, status = status,
                                                         b = b, db = db, bhazard = data[, bhazard],
                                                         control = list(maxit = 10000)))
  wh <- which.min(sapply(res_list, function(x) x$value))
  res <- res_list[[wh]]


  names(res$par) <- gsub("spline_", "", names(res$par))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(message) cat("Completed!\n")

  #Compute the hessian matrix
  if(hes){
    cov <- solve(pracma::hessian(flexible_minuslog_likelihood, res$par,
                                time = fu, status = status, b = b,
                                db = db, bhazard = data[, bhazard]))
    #cov <- solve(res$hessian)
  }else{
    cov <- NULL
  }

  #Output the results
  L <- list(formula = formula,
            data = data, cure = cure,
            coefs = res$par,
            base_function = base_function, dbase_function = dbase_function,
            knots = knots, knots.time = knots.time,
            ML = res$value, covariance = cov, df = length(res$par) - 1,
            tvc.formula = tvc.formula, formula = formula)

  class(L) <- "FlexRelsurvModel"
  L
}


print.FlexRelsurvModel <- function(fit){
  cat("Call:\n")
  print(fit$formula)
  if(!is.null(fit$tvc.formula)){
    cat("Call - time varying coefficients:\n")
    print(fit$tvc.formula)
  }
  cat("\nCoefficients:\n")
  print(fit$coefs)
}

summary.FlexRelsurvModel <- function(fit){
  se <- sqrt(diag(fit$cov))
  tval <- c(fit$coefs, fit$coefs.spline) / se
  coefs <- c(fit$coefs, fit$coefs.spline)
  TAB <- cbind(Estimate = fit$coefs,
                StdErr = se[1:length(fit$coefs)],
                t.value = tval[1:length(fit$coefs)],
                p.value = ifelse(is.na(tval[1:length(fit$coefs)]), rep(NA, length(fit$coefs)),
                                 2*pt(-abs(tval[1:length(fit$coefs)]), df = fit$df)))

  results <- list(coef = TAB)
  results$cure <- fit$cure
  results$ML <- fit$ML
  results$formula <- fit$formula
  results$formula.tvc <- fit$tvc.formula
  class(results) <- "summary.FlexRelsurvModel"
  results
}

print.summary.FlexRelsurvModel <- function(x)
{
  cat("Call:\n")
  print(x$formula)
  if(!is.null(x$tvc.formula)){
    cat("Call - time varying coefficients:\n")
    print(x$tvc.formula)
  }

  printCoefmat(x$coef, P.value = TRUE, has.Pvalue = T)
  cat("\n")
  cat("Cure =", x$cure, "\n")
  cat("LogLik(model) =", x$ML, "\n")
}

relsurv_fun_rs <- function(pars, M2, time, model){
  eta <- M2 %*% pars
  get.inv.link(exp(-exp(eta)), type = "relsurv")
}

predict.FlexRelsurvModel <- function(fit, newdata = NULL, time = NULL, type = "relsurv", ci = T, pars = NULL){
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
  if(type == "probcure"){
    if(!fit$cure){
      stop("Option 'probcure' only allowed for cure models")
    }
    max_knots <- max(exp(c(fit$knots, fit$knots.time)))
    last_time_check <- max(time) < max_knots
    if(last_time_check){
      time <- c(time, max_knots)
    }
  }

  all.coefs <- c(fit$coefs, fit$coefs.spline)
  b <- fit$base_function(knots = fit$knots, x = log(time))
  db <- fit$dbase_function(knots = fit$knots, x = log(time))
  if(!is.null(fit$knots.time)){
    tvc.b <- lapply(fit$knots.time, fit$base_function, x = log(time))
    tvc.db <- lapply(fit$knots.time, fit$dbase_function, x = log(time))
  }
  formula.2 <- formula(delete.response(terms(fit$formula)))
  M2.list <- lapply(1:nrow(newdata), function(i){
    model.matrix(formula.2, newdata[i, ,drop = F])[rep(1, nrow(b)),-1, drop = F]
  })

  M_list <- vector("list", nrow(newdata))
  for(i in 1:nrow(newdata)){
    M2 <- cbind(b, M2.list[[i]])
    dM2 <- cbind(db, M2.list[[i]])
    if(!is.null(fit$tvc.formula)){
      M_time <- model.matrix(fit$tvc.formula, newdata)[,-1, drop = F]
      M_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.b[[j]] * M_time[i, j]))
      db_time <- do.call(cbind, lapply(1:ncol(M_time), function(j) tvc.db[[j]] * M_time[i, j]))
      dM2 <- cbind(dM2, db_time)
      M2 <- cbind(M2, M_time)
    }
    M_list[[i]] <- list(M2 = M2, dM2 = dM2)
  }
  if(type %in% c("relsurv", "ehaz", "probcure")){
    if(type == "relsurv"){
      fun <- relsurv_fun_rs
    }else if(type == "ehaz"){
      fun <- ehaz_fun
    }else if(type == "probcure"){
      fun <- probcure_fun
    }

    rss <- vector("list", nrow(newdata))
    for(i in 1:nrow(newdata)){
      rss[[i]] <- data.frame(Est = c(fun(all.coefs, M2 = M_list[[i]]$M2, time = time)))
      if(ci){
        grads <- jacobian(fun, x = all.coefs, M2 = M_list[[i]]$M2, time = time)
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
      }else if (type == "probcure" & ci){
        l <- length(which(rss[[i]]$Est == 1))
        for(j in 1:l){
          rss[[i]][rss[[i]]$Est == 1,][j,] <- c(1, 0, 1, 1)
        }
        if(last_time_check){
          rss[[i]] <- rss[[i]][-nrow(rss[[i]]),]
        }
      }
    }
  }
  if(type == "probcure"){
    if(last_time_check){
      time <- time[-length(time)]
    }
  }
  return(list(res = rss, time = time, type = type))
}




plot.FlexRelsurvModel <- function(fit, newdata = NULL, time = NULL, ylim = c(0, 1), xlim = NULL,
                               xlab = "Time", ylab = "Relative survival", non.parametric = F,
                               type = "relsurv", col = 1, col.non.para = 2, ci = T,
                               include.knots = F, add = F, ...){

  ylab <- switch(type,
                 relsurv = "Relative survival",
                 ehaz = "Excess hazard",
                 probcure = "Conditional probability of cure",
                 crude_prob = "Probability of eventually dying from other causes than cancer",
                 lol = "Loss of lifetime")

  if(length(col) == 1 & !is.null(newdata)){
    col <- rep(col, nrow(newdata))
  }
  if(is.null(time)){
    if(is.null(xlim)){
      xlim <- c(0, max(fit$data$FU_years))
      time <- seq(xlim[1], xlim[2], length.out = 100)
    }
  }else{
    xlim <- range(time)
  }

  predict_rs <- predict(fit, newdata, time, type = type, ci = ci)
  nr.samples <- length(predict_rs$res)
  if(type == "ehaz"){
    ylim <- range(unlist(lapply(predict_rs$res, function(x) x[,-2])), na.rm = T)
  }

  for(i in 1:nr.samples){
    if(i == 1 & !add){
      plot(Est ~ time, data = predict_rs$res[[i]], type = "l", ylim = ylim, xlim = xlim,
           xlab = xlab, ylab = ylab, col = col[i], ...)
    }else{
      lines(Est ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], ...)
    }
    if(ci){
      lines(ci.upper ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], lty = 2, ...)
      lines(ci.lower ~ time, data = predict_rs$res[[i]], type = "l", col = col[i], lty = 2, ...)
    }
  }
  if(non.parametric){
    rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                     data = fit$data, ratetable = survexp.dk, method = "ederer2")
    rsfit$time <- rsfit$time / year
    lines(rsfit$surv ~ rsfit$time, type = "s", col = col.non.para, ...)
  }
  if(include.knots){
    abline(v = exp(fit$knots), lty = 2)
  }
}

