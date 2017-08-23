get_ini_values <- function(smooth.formula, tvc.formula, data, bhazard, hes, linkpi, linksu, formula, type,
                           fu, status, n.knots.time, knots, X, b, method = "mix"){
  vars <- c(all.vars(smooth.formula), all.vars(tvc.formula))
  if(method == "mix"){
    formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                             response = NULL)
    fit <- MixtureCureModel(formula, data = data, bhazard = bhazard, hes = F,
                            formula.k1 = formula.2, formula.k2 = ~ 1, type = type)
    pi_hat <- get_link("logistic")(X %*% fit$coefs[[1]])
    gpi_hat <- get_inv_link(linkpi)(pi_hat)
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    lp <- exp(model.matrix(formula.2, data = data) %*% fit$coefs[[2]])
    shat <- exp(-lp * fu ^ exp(fit$coefs[[3]]))
    gshat <- get_inv_link(linksu)(shat)
    fit_lm <- lm(gshat ~ -1 + b)
  }else if(method == "deaths"){
    formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                             response = formula[[2]])
    status2 <- 1 - status
    fit_glm <- glm(status2 ~ -1 + X, family = binomial(link = "logit"))
    pi_hat <- get_link("logistic")(predict(fit_glm))
    gpi_hat <- get_inv_link(linkpi)(pi_hat)
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    fit <- coxph(formula.2, data = data[status == 1,])
    cum_base_haz <- get_basehaz(fit)
    shat <- exp(-cum_base_haz$hazard) ^ exp(fit$linear.predictors)
    suppressWarnings(gshat <- get_inv_link(linksu)(shat))
    fit_lm <- lm(gshat ~ -1 + b[status == 1,])
  }else{
    tt <- terms(formula)
    formula.2 <- formula(delete.response(tt))
    vars <- unique(c("-1", all.vars(formula.2), all.vars(smooth.formula), all.vars(tvc.formula)))
    formula.3 <- reformulate(termlabels = vars, response = formula[[2]])
    fu_time <- all.vars(formula.3)[1]
    smooth.formula.paste <- as.formula(paste0("~basis(knots = knots, x = log(", fu_time, "))"))
    #fit <- stpm2(formula.3, data = data, smooth.formula = smooth.formula.paste, bhazard = data[, bhazard],
    #             tvc = n.knots.time)
    fit <- stpm2(formula.3, data = data, smooth.formula = smooth.formula.paste, bhazard = data[, bhazard])
    shat <- predict(fit, newdata = data, se.fit = F)
    gshat <- get_inv_link(linksu)(shat)
    data2 <- data
    data2[, fu_time] <- max(data2[, fu_time]) + 0.1
    pi_hat <- predict(fit, newdata = data2, se.fit = F)
    wh <- which(pi_hat >= shat)
    pi_hat[wh] <- shat[wh] - 0.01
    gpi_hat <- get_link(linkpi)(pi_hat)
    formula.logistic <- reformulate(termlabels = if(length(all.vars(formula.2)) == 0) "1" else all.vars(formula.2),
                                    response = "pred_pi")
    pi_fit <- lm(gpi_hat ~ -1 + X)
    ini_pi <- pi_fit$coefficients
    suhat <- (shat - pi_hat) / (1 - pi_hat)
    gsuhat <- get_inv_link(linksu)(suhat)
    fit_lm <- lm(gsuhat ~ -1 + b, data = data)
  }
  #Make proper naming of the initial values
  suppressWarnings(coefs <- summary(fit_lm)$coefficients[, "Estimate"])
  names(coefs)[1:length(knots)] <- paste0("gamma", 1:length(knots))
  fix_var <- all.vars(smooth.formula)
  if(length(fix_var) != 0){
    if(length(coefs) > length(knots)){
      names(coefs)[(length(knots) + 1):(length(knots) + length(fix_var))] <- paste0("spline_", fix_var)
    }
    if(length(coefs) > (length(knots) + length(fix_var))){
      names(coefs)[(length(knots) + length(fix_var) + 1):length(coefs)] <- paste0("gamma", 1:n.knots.time[[1]], ":",
                                                                                  all.vars(tvc.formula))
    }
  }else{
    if(length(coefs) > length(knots)){
      names(coefs)[(length(knots) + 1):length(coefs)] <- unlist(lapply(all.vars(tvc.formula),
                                                                       function(x) paste0("gamma",
                                                                                          1:n.knots.time[[x]],
                                                                                          ":", x)))
    }
  }
  ini_values <- c(ini_pi, coefs)
  names(ini_values)[1:ncol(X)] <- colnames(X)
  ini_values
}

FlexMixtureCureModel <- function(formula, data, bhazard, smooth.formula = ~ 1,
                                 knots = NULL, n.knots = NULL,
                                 tvc.formula = NULL, knots.time = NULL, n.knots.time = NULL,
                                 hes = T, message = T,
                                 type = "mixture", linkpi = "logistic", linksu = "loglog"){

  #Extract relevant variables
  fu <- eval(formula[[2]][[2]], envir = data)
  status <- eval(formula[[2]][[3]], envir = data)
  death_times <- fu[status == 1]

  #Caculate placement of knots and establish basis matrices
  if(is.null(knots)){
    bd_knots <- range(death_times)
    inner_knots <- quantile(death_times, 1 / (n.knots - 1)*1:(n.knots - 2))
    knots <- c(bd_knots, inner_knots)
    knots <- log(sort(knots))
  }

  b <- basis(knots = knots, log(fu))
  db <- dbasis(knots = knots, log(fu))

  if(!is.null(tvc.formula)){
    if(is.null(knots.time)){
      vars <- all.vars(tvc.formula)
      knots.time <- lapply(vars, function(x){
        bd_knots <- range(death_times)
        if(n.knots.time[[x]] > 2){
          inner_knots <- quantile(death_times, 1 / (n.knots.time[[x]] - 1)*1:(n.knots.time[[x]] - 2))
          log(sort(c(bd_knots, inner_knots)))
        }else{
          log(bd_knots)
        }
      })
      #names(knots.time) <- names(n.knots.time)
      b_list <- lapply(knots.time, basis, x = log(fu))
      db_list <- lapply(knots.time, dbasis, x = log(fu))
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
    b_time <- NULL
    db_time <- NULL
    tvc.formula <- NULL
  }


  #Construct design matrix
  X <- model.matrix(smooth.formula, data = data)
  b <- cbind(b, X[, -1], b_time)
  db <- cbind(db, matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)), db_time)
  X <- model.matrix(formula, data = data)

  #Generate initial values by running a mixture weibull cure model
  if(message) cat("Finding initial values...")

  #Extract link function
  link_fun_pi <- get_link(linkpi)
  link_fun_su <- get_link(linksu)
  dlink_fun_su <- get_dlink(linksu)

  types <- c("mix")
  ini_values <- lapply(types, function(x) get_ini_values(smooth.formula = smooth.formula,
                                                         tvc.formula =  tvc.formula, data = data,
                                                         bhazard = bhazard, hes = hes,
                                                         linkpi = linkpi, linksu = linksu,
                                                         formula = formula, type = type,
                                                         fu = fu, status = status,
                                                         n.knots.time = n.knots.time, knots = knots,
                                                         X = X, b = b, method = x))

  if(message) cat("Completed!\nFitting the model...")

  if(type == "mixture"){
    minusloglik <- flexible_mixture_minuslog_likelihood
  }else{
    minusloglik <- flexible_nmixture_minuslog_likelihood
  }

  # lapply(ini_values, function(inival) minusloglik(inival,
  #                                           time = fu, status = status, X = X,
  #                                           b = b, db = db, bhazard = data[, bhazard],
  #                                           link_fun_pi = link_fun_pi,
  #                                           link_fun_su = link_fun_su,
  #                                           dlink_fun_su = dlink_fun_su))
  #
  # rs_fit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
  #                   data = data, ratetable = survexp.dk, method = "ederer1")
  # rs_fit$time <- rs_fit$time / year
  # plot(rs_fit)
  #
  # pi <- exp(ini_values[[1]][1]) / (exp(ini_values[[1]][1]) + 1)
  # f <- function(t) pi + (1 - pi) * exp(-exp(basis(knots, x = log(t)) %*% ini_values[[1]][-1]))
  # curve(f, from = 0, to = 16, col = 2, add = T)
  #
  # pi <- exp(ini_values[[2]][1]) / (exp(ini_values[[2]][1]) + 1)
  # f <- function(t) pi + (1 - pi) * exp(-exp(basis(knots, x = log(t)) %*% ini_values[[2]][-1]))
  # curve(f, from = 0, to = 16, col = 3, add = T)
  #
  # pi <- exp(ini_values[[3]][1]) / (exp(ini_values[[3]][1]) + 1)
  # f <- function(t) pi + (1 - pi) * exp(-exp(basis(knots, x = log(t)) %*% ini_values[[3]][-1]))
  # curve(f, from = 0, to = 16, col = 4, add = T)

  #Fit the model
  res_list <- lapply(ini_values, function(inival) optim(par = inival,
                                                        fn = minusloglik,
                                                        time = fu, status = status, X = X,
                                                        b = b, db = db, bhazard = data[, bhazard],
                                                        link_fun_pi = link_fun_pi,
                                                        link_fun_su = link_fun_su,
                                                        dlink_fun_su = dlink_fun_su,
                                                        control = list(maxit = 10000)))
  MLs <- sapply(res_list, function(x) x$value)
  wh <- which.min(MLs)
  res <- res_list[[wh]]
  # pi <- exp(res_list[[1]]$par[1]) / (exp(res_list[[1]]$par[1]) + 1)
  # f <- function(t) pi + (1 - pi) * exp(-exp(basis(knots, x = log(t)) %*% res_list[[1]]$par[-1]))
  # curve(f, from = 0, to = 16, col = 5, add = T)
  # res_list[[1]]$value

  names(res$par) <- gsub("spline_", "", names(res$par))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(message) cat("Completed!\n")

  #Compute the hessian matrix
  if(hes){
    cov <- solve(pracma::hessian(minusloglik, res$par,
                                 time = fu, status = status,
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
            ML = res$value, covariance = cov, tvc.formula = tvc.formula, formula = formula,
            formula_main = smooth.formula, type = type,
            link_fun_pi = link_fun_pi,
            link_fun_su = link_fun_su,
            dlink_fun_su = dlink_fun_su,
            df = length(res$par) - 1, MLs = MLs)

  class(L) <- "FlexCureModel"
  L
}

print.FlexCureModel <- function(fit){
  cat("Call pi:\n")
  print(fit$formula)
  cat("Call S_u(t):\n")
  print(fit$formula_main)
  cat("\nCoefficients:\n")
  print(list(pi = fit$coefs,
             s_ut = fit$coefs.spline))
}

summary.FlexCureModel <- function(fit){
  se <- sqrt(diag(fit$cov))
  tval <- c(fit$coefs, fit$coefs.spline) / se
  coefs <- c(fit$coefs, fit$coefs.spline)
  TAB1 <- cbind(Estimate = fit$coefs,
               StdErr = se[1:length(fit$coefs)],
               t.value = tval[1:length(fit$coefs)],
               p.value = ifelse(is.na(tval[1:length(fit$coefs)]), rep(NA, length(fit$coefs)),
                                2*pt(-abs(tval[1:length(fit$coefs)]), df = fit$df)))

  TAB2 <- cbind(Estimate = fit$coefs.spline,
               StdErr = se[1:length(fit$coefs.spline)],
               t.value = tval[1:length(fit$coefs.spline)],
               p.value = ifelse(is.na(tval[1:length(fit$coefs.spline)]), rep(NA, length(fit$coefs.spline)),
                                2*pt(-abs(tval[1:length(fit$coefs.spline)]), df = fit$df)))


  results <- list(pi = TAB1, S_u = TAB2)
  results$type <- fit$type
  results$link <- fit$link
  results$ML <- fit$ML
  results$formula <- fit$formula
  results$formula.fix <- fit$formula_main
  results$formula.tvc <- fit$tvc.formula
  class(results) <- "summary.CureModel"
  results
}

print.summary.CureModel <- function(x)
{
  cat("Call - pi:\n")
  print(x$formula)
  #    cat("\n")
  printCoefmat(x$pi, P.value = TRUE, has.Pvalue = T)
  cat("\nCall - S_u - baseline: ")
  print(as.formula(deparse(x$formula.fix)))
  if(length(all.vars(x$formula.tvc))){
    cat("Call - S_u - tvc: ")
    print(deparse(x$formula.tvc))
  }
  printCoefmat(x$S_u, P.value = TRUE, has.Pvalue = T)
  cat("\n")
  cat("Type =", x$type, "\n")
  cat("Link =", x$link, "\n")
  cat("LogLik(model) =", x$ML, "\n")

}

get.link <- function(x, type){
  if(type %in% c("curerate", "relsurv", "crudeprob")){
    exp(x) / (exp(x) + 1)
  }else if (type == "probcure"){
    pnorm(x)
  } else if(type == "ehaz"){
    x
  }
}

get.inv.link <- function(x, type){
  if(type %in% c("curerate", "relsurv", "crudeprob")){
    log(x / (1 - x))
  }else if(type == "probcure"){
    qnorm(x)
  }else if(type == "ehaz"){
    x
  }
}

relsurv_fun <- function(pars, M2, M, dM2, time, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  pi <- c(link_fun_pi(pi_fun(pars, M)))
  eta <- M2 %*% pars[-c(1:ncol(M))]
  if(model == "mixture"){
    get.inv.link(pi + (1 - pi) * link_fun_su(eta), type = "relsurv")
  }
  else if(model == "nmixture"){
    get.inv.link(pi ^ link_fun_su(eta), type = "relsurv")
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
    get.inv.link(-log(pi) * ds_u * deta / time, type = "ehaz")
  }
}

probcure_fun <- function(pars, M2, M, dM2, time, pi, pi_fun, model, link_fun_pi, link_fun_su, dlink_fun_su){
  if(model == "mixture"){
    pi <- c(get.link(pi_fun(pars, M), type = "curerate"))
    eta <- M2 %*% pars[-c(1:ncol(M))]
    get.inv.link(pi / (pi + (1 - pi) * exp(-exp(eta))), type = "probcure")
  }else if(model == "nmixture"){
    eta <- M2 %*% pars
    pi <- exp(-exp(eta[length(eta)]))
    get.inv.link(pi / ifelse(time == 0, 1, exp(-exp(eta))), type = "probcure")
  }
}

predict.FlexCureModel <- function(fit, newdata = NULL, time = NULL, type = "relsurv", ci = T, pars = NULL){
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
    pi <- data.frame(pi = get.link(pi, type = "curerate"))
    if(ci){
      grads <- jacobian(pi_fun, x = all.coefs, M = M)
      pi$var <- apply(grads, MARGIN = 1, function(x) x %*% fit$covariance %*% x)
      pi$ci.lower <- get.link(pi$pi - qnorm(0.975) * sqrt(pi$var), type = type)
      pi$ci.upper <- get.link(pi$pi + qnorm(0.975) * sqrt(pi$var), type = type)
    }
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
    if(type %in% c("relsurv", "ehaz", "probcure")){
      if(type == "relsurv"){
        fun <- relsurv_fun
      }else if(type == "ehaz"){
        fun <- ehaz_fun
      }else if(type == "probcure"){
        fun <- probcure_fun
      }

      rss <- vector("list", nrow(newdata))
      for(i in 1:nrow(newdata)){
        rss[[i]] <- data.frame(Est = c(fun(all.coefs, M_list[[i]]$M2, M[i,, drop = FALSE], M_list[[i]]$dM2,
                                           time, pi_fun = pi_fun, model = "mixture",
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

        if(type == "relsurv"){
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


plot.FlexCureModel <- function(fit, newdata = NULL, time = NULL, ylim = c(0, 1), xlim = NULL,
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


