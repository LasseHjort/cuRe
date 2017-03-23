#Extract expected hazard

retrieve_background_haz <- function(formula, data){
  expected <- survexp(formula, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  exp_haz <- rep(0, nrow(data))
  D <- as.data.frame(cbind(time = c(0, expected$time, 100), surv = c(1, expected$surv, 0)))
  ll <- loess(surv ~ time, data = D)
  S_explog <- function(t) log(predict(ll, newdata = t))
  h_exp <- function(t) -numDeriv::grad(S_explog, t)
  exp_haz <- h_exp(data$FU / year)
  exp_haz
}

retrieve_background_haz_cov <- function(formula, data){
  expected <- survexp(formula, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  exp_haz <- rep(0, nrow(data))
  if(var != "1"){
    D <- as.data.frame(cbind(time = c(0, expected$time, 100), rbind(1, expected$surv, 0)))
    for(i in 1:nlevels(data[, var])){
      D_new <- D[, c("time", paste0(var, "=", levels(data[, var])[i]))]
      names(D_new) <- c("time", "surv")
      ll <- loess(surv ~ time, data = D_new)
      S_explog <- function(t) log(predict(ll, newdata = t))
      h_exp <- function(t) -numDeriv::grad(S_explog, t)
      wh <- data[, var] == levels(data[, var])[i]
      exp_haz[wh] <- h_exp(data$FU[wh] / year)
    }
  }else{
    D <- as.data.frame(cbind(time = c(0, expected$time, 100), surv = c(1, expected$surv, 0)))
    ll <- loess(surv ~ time, data = D)
    S_explog <- function(t) log(predict(ll, newdata = t))
    h_exp <- function(t) -numDeriv::grad(S_explog, t)
    exp_haz <- h_exp(data$FU / year)
  }
  exp_haz
}

extract_background <- function(t, data){
  Time <- seq(0, to = max(t), length.out = 100)
  haz <- rep(NA, length(Time))
  sex_new <- as.character(data$sex)
  for(i in 1:length(haz)){
    cat(i, "\n")
    age_new <- pmin(round((data$age + Time[i] * year) / year), 99) + 1
    diag_date_new <- format(data$diag_date + Time[i] * year, "%Y")
    wh <- which(data$OS[,1] >= Time[i])
    tmp_haz <- sapply(wh, function(x){
      survexp.dk[age_new[x], data$sex[x], diag_date_new[x]]
    })
    haz[i] <- mean(tmp_haz) * year
  }
  # plot(haz ~ Time, type = "l")
  # D <- data.frame(time = Time, haz = haz)
  ll <- loess(haz ~ time, data = D, span = 0.1)
  # a <- predict(ll, newdata = data.frame(time  = Time))
  # lines(a ~ Time)
  predict(ll, newdata = data.frame(time = data$OS[,1]))
}


extract_general <- function(data){
  haz <- rep(NA, nrow(data))
  sex_new <- as.character(data$sex)
  age_new <- pmin(round((data$age + data$FU) / year), 99)
  year_eval <- format(data$diag_date + data$FU, "%Y")
  for(i in 1:nrow(data)){
    haz[i] <- survexp.dk[age_new[i], sex_new[i], year_eval[i]]
  }
  haz <- haz * year
}


#Densities
#Weibull
wei_surv <- function(t, shape, scale, X) exp(-t ^ shape * scale)
wei_dens <- function(t, shape, scale){
  exp(-t ^ shape * scale) * scale * shape * t ^ (shape - 1)
}

slnorm <- function(t, meanlog, sdlog) 0



#Spline functions
basis_cure <- function(knots, x){
  nk <- length(knots)
  b <- matrix(nrow = length(x), ncol = nk - 1)
  knots_rev <- rev(knots)
  if (nk > 0) {
    b[, 1] <- 1
  }
  if (nk > 2) {
    for (j in 2:(nk - 1)) {
      lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
      b[, j] <- pmax(knots_rev[nk - j + 1] - x, 0)^3 - lam * pmax(knots_rev[nk] - x, 0)^3 -
        (1 - lam) * pmax(knots_rev[1] - x, 0)^3
    }
  }
  b
}

dbasis_cure <- function(knots, x){
  nk <- length(knots)
  b <- matrix(nrow = length(x), ncol = nk - 1)
  knots_rev <- rev(knots)
  if (nk > 0) {
    b[, 1] <- 0
  }
  if (nk > 2) {
    for (j in 2:(nk - 1)) {
      lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
      b[, j] <- - 3 * pmax(knots_rev[nk - j + 1] - x, 0)^2 + 3 * lam * pmax(knots_rev[nk] - x, 0)^2 +
        3 * (1 - lam) * pmax(knots_rev[1] - x, 0)^2
    }
  }
  b
}

spline_model_haz_spline_true <- function(t, lambda1, lambda2, knots, cov){
  base_functions <- sapply(1:(length(knots) - 2), base_function_cure_deriv, x = log(t), knots = knots)
  val1 <- -log(spline_model_surv_spline(t, lambda1, lambda2, knots, cov))
  val2 <- (base_functions %*% lambda1[-1] + base_functions %*% lambda2[-1] * cov) / t
  val <- val1 * val2
  c(val)
}

flex_spline_surv <- function(t, lambda1, lambda2, knots, z){
  base_functions <- basis(log(t), knots = knots)
  val <- exp(-exp(base_functions %*% lambda1 + base_functions %*% lambda2 * z))
  val[t == 0] <- 1
  c(val)
}

flex_spline_haz <- function(t, lambda1, lambda2, knots, z){
  dbase_functions <- dbasis(log(t), knots = knots)
  (dbase_functions %*% lambda1 + dbase_functions %*% lambda2 * z) / t
}


#Log-likelihood functions
minusloglike_parametricCure <- function(param, X, status, futime, surv, dens, bhazard, model){
  gamma <- param[1:ncol(X)]
  lp <- X %*% gamma
  pi <- exp(lp) / (exp(lp) + 1)
  theta <- param[(ncol(X) + 1):length(param)]
  scale <- exp(X %*% theta[1:ncol(X)])
  shape <- exp(theta[length(theta)])
  surv_eval <- surv(t = futime, shape = shape, scale = scale)
  dens_eval <- dens(t = futime, shape = shape, scale = scale)
  if(model == "mixture"){
    first_term <- status * log( bhazard + dens_eval * ( 1 - pi ) / ( pi + (1 - pi) * surv_eval ) )
    second_term <- log( pi + (1 - pi) * surv_eval )
  }else{
    first_term <- status * log( bhazard - log(pi) * dens_eval)
    second_term <- log( pi ) - log(pi) * surv_eval
  }
  -sum(first_term + second_term)
}


log_likelihood_flex_surv <- function(param, data, knots){
  wh <- which(data$status == 1)
  bs <- basis(knots = knots, x = log(data$OS[,1]))
  dbs <- dbasis(knots = knots, x = log(data$OS[wh, 1]))
  rsurv <- exp(bs %*% param)
  haz <- dbs %*% param * rsurv[wh] / data$OS[wh, 1]
  terms <- rep(0, nrow(data))
  terms[wh] <- log(haz) - rsurv[wh]
  terms[-wh] <- -rsurv[-wh]
  #cat(min(haz), "\n")
  #terms <- data$status * log( data$exp_haz +  haz) + log(surv)
  -sum(terms)
}

log_likelihood_flex <- function(param, data, knots){
  surv <- log(psurvspline(data$OS[,1], gamma = param, knots = knots, lower.tail = FALSE))
  wh <- which(data$status == 1)
  haz <- hsurvspline(data$OS[wh, 1], gamma = param, knots = knots)
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + haz)
  #cat(min(haz), "\n")
  #terms <- data$status * log( data$exp_haz +  haz) + log(surv)
  -sum(surv)
}

log_likelihood_flex2 <- function(param, data, knots){
  wh <- which(data$status == 1)
  bs <- basis(knots = knots, x = log(data$OS[,1]))
  dbs <- dbasis(knots = knots, x = log(data$OS[wh, 1]))
  rsurv <- exp(bs %*% param)
  haz <- dbs %*% param * rsurv[wh] / data$OS[wh, 1]
  surv <- -rsurv
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + haz)
  #cat(min(haz), "\n")
  #terms <- data$status * log( data$exp_haz +  haz) + log(surv)
  -sum(surv)
}

#b <- param + rnorm(5, sd = 0.01)
#log_likelihood_flex(b, data, knots)
#log_likelihood_flex2(b, data, knots)

dlog_likelihood_flex <- function(param, data, knots){
  wh <- which(data$status == 1)
  bs <- basis(knots = knots, x = log(data$OS[,1]))
  dbs <- dbasis(knots = knots, x = log(data$OS[wh, 1]))
  etas <- bs %*% param
  exp_etas <- exp(etas)
  detas <- dbs %*% param
  haz <- detas / data$OS[wh, 1] * exp_etas[wh]
  logs <- 1 / (data$exp_haz[wh] + haz)
  sapply(1:length(knots), function(x){
     dhaz1 <- exp_etas[wh] * bs[wh, x] * detas + exp_etas[wh] * dbs[, x]
     dhaz2 <- logs * dhaz1 / data$OS[wh, 1]
     survs <- -exp_etas * bs[, x]
     survs[wh] <- survs[wh] + dhaz2
     -sum(survs)
  })
}

#dlog_likelihood_flex(param, data, knots)
#numDeriv::grad(log_likelihood_flex, x = param, knots = knots, data = data)

log_likelihood_flex_spline <- function(param, data, knots, cov){
  base_coefs <- param[1:length(knots)]
  spline_coefs <- param[(length(knots) + 1):length(param)]
  surv <- log(flex_spline_surv(t = data$OS[,1],
                                   lambda1 = base_coefs,
                                   lambda2 = spline_coefs,
                                   knots = knots,
                                   z = data[, cov]))
  wh <- which(data$status == 1)
  haz <- -flex_spline_haz(t = data$OS[wh, 1],
                                 lambda1 = base_coefs,
                                 lambda2 = spline_coefs,
                                 knots = knots,
                                 z = data[wh, cov]) * surv[wh]
  #cat(min(haz), "\n")
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + haz)
  -sum(surv)
}

log_likelihood_andersson <- function(param, data, knots){
  surv <- log(spline_model_surv(t = data$OS[,1], lambda = param, knots = knots))
  wh <- which(data$status == 1)
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + spline_model_haz(t = data$OS[wh, 1], lambda = param, knots = knots))
  -sum(surv)
}

log_likelihood_andersson <- function(param, data, knots){
  wh <- which(data$status == 1)
  b <- basis_cure(knots, log(data$OS[,1]))
  db <- dbasis_cure(knots, log(data$OS[wh,1]))
  surv <- -exp(b %*% param)
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + exp(b[wh,] %*% param) * db %*% param / data$OS[wh, 1])
  -sum(surv)
}

dlog_likelihood_andersson <- function(param, data, knots){
  wh <- which(data$status == 1)
  bs <- basis_cure(knots = knots, x = log(data$OS[,1]))
  dbs <- dbasis_cure(knots = knots, x = log(data$OS[wh, 1]))
  etas <- bs %*% param
  exp_etas <- exp(etas)
  detas <- dbs %*% param
  haz <- detas / data$OS[wh, 1] * exp_etas[wh]
  logs <- 1 / (data$exp_haz[wh] + haz)
  sapply(1:(length(knots) - 1), function(x){
    dhaz1 <- exp_etas[wh] * bs[wh, x] * detas + exp_etas[wh] * dbs[, x]
    dhaz2 <- logs * dhaz1 / data$OS[wh, 1]
    survs <- -exp_etas * bs[, x]
    survs[wh] <- survs[wh] + dhaz2
    -sum(survs)
  })
}


log_likelihood_andersson_spline <- function(param, data, knots, cov){
  base_coefs <- param[1:(length(knots) - 1)]
  spline_coefs <- param[length(knots):length(param)]
  surv <- log(spline_model_surv_spline(t = data$OS[,1],
                                   lambda1 = base_coefs,
                                   lambda2 = spline_coefs,
                                   knots = knots,
                                   cov = data[, cov]))
  wh <- which(data$status == 1)
  haz <- -spline_model_haz_spline(t = data$OS[wh, 1],
                                 lambda1 = base_coefs,
                                 lambda2 = spline_coefs,
                                 knots = knots,
                                 cov = data[wh, cov]) * surv[wh]
  surv[wh] <- surv[wh] + log(data$exp_haz[wh] + haz)
  -sum(surv)
}


#Functions for extracting initial value

initialParaCure <- function(formula, data, X, surv, dens, model){
  expected <- survexp(FU ~ 1, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data, scale = year, method = "conditional")
  response <- all.vars(formula)[1]
  fit <- coxph(formula, data = data)
  base <- basehaz(fit)
  base$surv <- exp(-base$hazard)
  pi_ini <- 0.5
  if(model == "mixture"){
    survs <- (base$surv / summary(expected, times = base$time)$surv - pi_ini) / (1 - pi_ini)
  }else{
    relsurv <- base$surv / summary(expected, times = base$time)$surv
    survs <- 1 - log(relsurv) / log(pi_ini)
  }
  wei <- function(t, shape, scale) exp(-scale * t ^ shape)
  fit <- nls(survs ~ wei(base$time, shape, scale),
             start = list(shape = 1, scale = 1 ))
  coefs <- summary(fit)$coefficients[, "Estimate"]
  #f <- function(t) wei(t, shape = coefs["shape"], scale = coefs["scale"])
  #curve(f, col = 2, add = T)
  c(log(pi_ini / (1 - pi_ini)), rep(0, ncol(X) - 1),
    log(coefs["scale"]), rep(0, ncol(X) - 1),
    log(coefs["shape"]))
}

get_starting_values_flex_surv <- function(formula, data, knots){
  fit <- coxph(formula, data = data[data$OS[,2] == 1,])
  sfit <- survfit(fit)
  times <- sort(data$OS[data$OS[, 2] == 1,1])
  logH <- log(-log(summary(sfit, times)$surv))
  b <- basis(knots, log(times))
  db <- dbasis(knots, log(times))
  Dmat <- t(b) %*% b
  dvec <- t(t(logH) %*% b)
  Amat <- t(db)
  eps <- rep(1e-09, length(logH))
  solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = eps)$solution
}


get_starting_values_flex <- function(formula, data, knots){
  fit <- coxph(formula, data = data[data$OS[,2] == 1,])
  sfit <- survfit(fit)
  times <- sort(data$OS[data$OS[, 2] == 1,1])
  formula2 <- as.formula(paste0("FU ~ 1"))
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  logH <- log(-log(summary(sfit, times)$surv / summary(expected, times)$surv))
  b <- basis(knots, log(times))
  db <- dbasis(knots, log(times))
  Dmat <- t(b) %*% b
  dvec <- t(t(logH) %*% b)
  Amat <- t(db)
  eps <- rep(1e-09, length(logH))
  solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = eps)$solution
}

get_starting_values_flex_cov <- function(formula, data, X, knots){
  sfit <- survreg(formula, data = data)
  formula2 <- as.formula(paste0("FU ~ ", var))
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  times <- sort(data$OS[,1])
  survs <- summary(sfit, times)$surv / summary(expected, times)$surv
  fit <- nls(survs ~ psurvspline(times, gamma, knots = knots, lower.tail = FALSE),
             start = list(gamma = rep(0, length(knots))), control = list(warnOnly = T))
  summary(fit)$coefficients[, "Estimate"]
}


get_starting_values_flex_cov <- function(formula, data, X, knots){
  fit <- coxph(formula, data = data)
  newdata <- as.data.frame(matrix(rep(0, 1), nrow = 1))
  names(newdata) <- names(fit$coefficients)
  sfit <- survfit(fit, newdata = newdata)

  formula2 <- as.formula(paste0("FU ~ ", var))
  expected <- survexp(formula, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  times <- sort(data$OS[,1])
  survs <- summary(sfit, times)$surv / summary(expected, times)$surv
  fit <- nls(survs ~ psurvspline(times, gamma, knots = knots, lower.tail = FALSE),
             start = list(gamma = rep(0, length(knots))), control = list(warnOnly = T))
  summary(fit)$coefficients[, "Estimate"]
}


get_starting_values_flex_cov <- function(formula, data, X, knots){
  fit <- rsadd(Surv(FU, status) ~ age_years + ratetable(age = age, sex = sex, year = diag_date),
               data = DLBCL, ratetable = survexp.dk, int = 5, method = "EM")
  sm <- epa(fit, times = DLBCL$FU[DLBCL$status == 1])
  plot(sm$times, sm$lambda)
  times <- sm$times / year
  #plot(sm$times, sm$lambda, type = "l")
  logH <- log(sm$lambda)
  b <- basis(knots, log(times))
  db <- dbasis(knots, log(times))
  Dmat <- t(b) %*% b
  dvec <- t(t(logH) %*% b)
  Amat <- t(db)
  eps <- rep(1e-09, length(logH))
  solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = eps)$solution
}

get_starting_values_flex_cure <- function(formula, data, knots){
  sfit <- survfit(formula, data = data)
  formula2 <- as.formula("FU ~ 1")
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  pi_ini <- 0.5
  times <- seq(min(data$OS[,1]), max(data$OS[, 1]), length.out = 100)
  survs <- summary(sfit, times)$surv / summary(expected, times)$surv
  f <- function(times, gamma, knots) exp(-exp(basis_cure(knots, log(times)) %*% gamma))
  fit <- nls(survs ~ f(times, gamma, knots = knots),
             start = list(gamma = rep(0, length(knots) - 1)))
  summary(fit)$coefficients[, "Estimate"]
}

get_starting_values_flex_cure <- function(formula, data, knots){
  sfit <- survfit(formula, data = data)
  formula2 <- as.formula("FU ~ 1")
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data,
                      method = "conditional", ratetable = survexp.dk, scale = year)
  times <- seq(min(data$OS[,1]), max(data$OS[, 1]), length.out = 100)
  logH <- log(-log(summary(sfit, times)$surv / summary(expected, times)$surv))
  f <- function(times, gamma, knots) exp(-exp(basis_cure(knots, log(times)) %*% gamma))
  fit <- nls(logH ~ basis_cure(knots, x = log(times)) %*% gamma,
             start = list(gamma = rep(0, length(knots) - 1)))
  summary(fit)$coefficients[, "Estimate"]
}

get_starting_values_flex_cure <- function(formula, data, knots){
  fit <- coxph(formula, data = data[data$OS[,2] == 1,])
  sfit <- survfit(fit)
  times <- sort(data$OS[data$OS[, 2] == 1,1])
  formula2 <- as.formula(paste0("FU ~ 1"))
  expected <- survexp(formula2, rmap = list(age = age, sex = sex, year = diag_date),
                      data = data[data$OS[,2] == 1, ],
                      method = "conditional", ratetable = survexp.dk, scale = year)
  logH <- log(-log(summary(sfit, times)$surv / summary(expected, times)$surv))
  # b <- basis_cure(knots, log(times))
  # db <- dbasis_cure(knots, log(times))
  # Dmat <- t(b) %*% b
  # dvec <- t(t(logH) %*% b)
  # Amat <- t(db)
  # eps <- rep(1e-09, length(logH))
  # solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = eps)$solution

  fit <- nls(logH ~ basis_cure(knots, x = log(times)) %*% gamma,
             start = list(gamma = rep(0, length(knots) - 1)))
  summary(fit)$coefficients[, "Estimate"]
}
