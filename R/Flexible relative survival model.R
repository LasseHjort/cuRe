FlexModels <- function(data, ratetable = survexp.dk, model = "flex", knots = NULL, 
                       hes = T, n.knots = 5, last_knot = 10, add.knots = NULL){
  if(is.null(knots)){
    ti <- data$OS[data$OS[,2] == 1,1]
    knots <- c(log(min(ti)), quantile(log(ti), 1 / (n.knots - 1)*1:(n.knots - 2)), log(last_knot)) 
    if(!is.null(add.knots)){
      knots <- sort(c(knots, log(add.knots)))
    }
  }
  
  if(model == "andersson"){
    get_ini_val <- get_starting_values_flex_cure
    minusloglik <- log_likelihood_andersson
    dminusloglik <- dlog_likelihood_andersson
  }else{
    get_ini_val <- get_starting_values_flex
    minusloglik <- log_likelihood_flex2
    dminusloglik <- dlog_likelihood_flex
  }
  
  ini_values <- get_ini_val(formula = OS ~ 1, data = data, knots = knots)
  res <- optim(ini_values, fn = minusloglik,
               data = data, control = list(maxit = 1000),
               knots = knots)
  
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  
  if(hes){
    hessian <- solve(numDeriv::jacobian(dminusloglik, res$par, knots = knots, data = data))
    list(Data = data, var = var, model = model,
         Lambda = res$par, knots = knots, ML = res$value, Hessian = hessian)
  }else{
    list(Data = data, var = var, model = model,
         Lambda = res$par, knots = knots, ML = res$value)
  }
}

plot.FlexModels <- function(object, main = "", col = NULL, ci = TRUE, conf.int = F){
  if(!is.null(object$model)) 
    object <- list(object)
  if(is.null(names(object)))
    names(object) <- 1:length(object)
  
  rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date), 
                   data = object[[1]]$Data, ratetable = survexp.dk, method = "ederer2")
  plot(rsfit, conf.int = conf.int, xscale = year, ylab = "Relative survival", 
       xlab = "Follow-up time (years)", axes = F, main = main)
  grid()
  box()
  axis(1, at = seq(0, 16 * year, by = 2 * year), labels = seq(0, 16, by = 2))
  axis(2, at = seq(0, 1, by = 0.2))
  times <- seq(0, max(object[[1]]$Data$FU), length.out = 1000)
  if(is.null(col)) col <- 2:(length(object) + 1)
  for(i in 1:length(object)){
    fit <- object[[i]]
    if(fit$model == "andersson"){
      basis_function <- basis_cure
    }else{
      basis_function <- basis
    }
    basis_eval <- basis_function(fit$knots, x = log(times / year))
    surv_eval <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda)))
    lines(surv_eval ~ times, col = col[i])
    if(ci){
      VAR <- sapply(1:length(times), function(x){
        basis_eval[x,] %*% fit$Hessian %*% basis_eval[x,]
      })
      ci_upper <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda + sqrt(VAR) * qnorm(0.975))))
      ci_lower <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda - sqrt(VAR) * qnorm(0.975))))
      lines(ci_upper ~ times, col = col[i], lty = 2)
      lines(ci_lower ~ times, col = col[i], lty = 2)
    }
  }
}

eval.FlexModels <- function(object, times){
  if(!is.null(object$model)) 
    object <- list(object)
  if(is.null(names(object)))
    names(object) <- 1:length(object)
  
  res <- lapply(object, function(fit){
    if(fit$model == "andersson"){
      basis_function <- basis_cure
    }else{
      basis_function <- basis
    }
    basis_eval <- basis_function(fit$knots, log(times))
    VAR <- sapply(1:length(times), function(x){
      basis_eval[x,] %*% fit$Hessian %*% basis_eval[x,]
    })
    surv_eval <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda)))
    ci_upper <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda + sqrt(VAR) * qnorm(0.975))))
    ci_lower <- ifelse(times == 0, 1, exp(-exp(basis_eval %*% fit$Lambda - sqrt(VAR) * qnorm(0.975))))
    data.frame(Time = times, Surv = surv_eval, "Lower95" = ci_lower, "Upper95" = ci_upper)
  })
  if(length(res) == 1){
    return(res[[1]])
  }else{
    return(res)
  }
}

cure_fraction.FlexModels <- function(fit){
  if(fit$model != "andersson"){
    stop("Wrong type of model")
  }
  est <- exp(-exp(fit$Lambda[1]))
  var <- fit$Hessian[1,1]
  conf.lo <- exp(-exp(fit$Lambda[1] + sqrt(var) * qnorm(0.975)))
  conf.up <- exp(-exp(fit$Lambda[1] - sqrt(var) * qnorm(0.975)))
  res <- c(est, conf.lo, conf.up)
  names(res) <- c("est", "conf.lo", "conf.up")
  res
}

fit_DLBCL <- FlexModels(data = DLBCL, last_knot = 12)
fit_DLBCL2 <- FlexModels(data = DLBCL, model = "andersson", add.knots = 10, last_knot = 12)
fit_FL <- FlexModels("1", data = FL, last_knot = 12)
fit_FL2 <- FlexModels("1", data = FL, model = "andersson", add.knots = 10, last_knot = 12)
fit_PTCL <- FlexModels("1", data = PTCL, last_knot = 12)
fit_PTCL2 <- FlexModels(data = PTCL, model = "andersson", add.knots = 10, last_knot = 12)

pdf(file.path(fig.out, "RSCombined.pdf"), width = 11, height = 4)
m <- matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.9,0.1))
par(mai = mai_par, mar = c(5.1, 4.1, 4.1, 2.1))

plot.FlexModels(list(fit_DLBCL, fit_DLBCL2), main = "DLBCL", ci = F, conf.int = T, col = c("brown2", "darkolivegreen3"))
plot.FlexModels(list(fit_FL, fit_FL2), main = "FL", ci = F, conf.int = T, col = c("brown2", "darkolivegreen3"))
plot.FlexModels(list(fit_PTCL, fit_PTCL2), main = "PTCL", ci = F, conf.int = T, col = c("brown2", "darkolivegreen3"))

par(mar=c(0,0, 0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend = c("Ederer II estimate", "Nelson et al.", "Andersson et al."), 
       col = c("black", "brown2", "darkolivegreen3"), lwd=5, cex=.9, horiz = TRUE)
dev.off()

FlexModels_spline <- function(data, cov, ratetable = survexp.dk, model = "flex", 
                              knots = NULL, hes = T, n.knots = 5, last_knot = 10){
  if(is.null(knots)){
    ti <- data$OS[data$OS[,2] == 1,1]
    knots <- c(log(min(ti)), quantile(log(ti), 1 / (n.knots - 1)*1:(n.knots - 2)), log(last_knot)) 
  }
  
  if(model == "andersson"){
    get_ini_val <- get_starting_values_flex_cure
    minusloglik <- log_likelihood_andersson_spline
  }else{
    get_ini_val <- get_starting_values_flex
    minusloglik <- log_likelihood_flex_spline
  }
  
  ini_values <- get_ini_val(OS ~ 1, data = data, knots = knots)
  ini_values <- c(ini_values, rep(0, length(knots)))
  res <- optim(ini_values, fn = minusloglik,
               data = data, control = list(maxit = 3000),
               knots = knots, cov = cov)
  
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  
  if(hes){
    hessian <- solve(numDeriv::hessian(minusloglik, res$par, knots = knots, data = data, cov = cov))
    list(Data = data, model = model, convergence = res$convergence, Lambda1 = res$par[1:length(knots)], 
         Lambda2 = res$par[(length(knots) + 1):length(res$par)], 
         knots = knots, ML = res$value, Hessian = hessian) 
  }else{
    list(Data = data, model = model, convergence = res$convergence, Lambda1 = res$par[1:length(knots)], 
         Lambda2 = res$par[(length(knots) + 1):length(res$par)], 
         knots = knots, ML = res$value)   
  }
}

fit_PTCL <- FlexModels_spline(data = PTCL, cov = "age_years")
fit_PTCL <- FlexModels_spline(data = PTCL, cov = "age_years", model = "andersson")


FlexrelsurvModel_spline <- function(var, data, ratetable = survexp.dk, knots = NULL, hes = TRUE, n.knots = 5){
  if(is.null(knots)){
    ti <- data$OS[data$OS[,2] == 1,1]
    knots <- c(log(min(ti)), quantile(log(ti), 1 / (n.knots - 1)*1:(n.knots - 2)), log(10)) 
  }
  
  ini_values <- get_starting_values_flex(OS ~ 1, data, knots = knots)
  ini_values <- c(ini_values, rep(0, length(knots)))# + rnorm(6, sd = 0.0001)
  res <- optim(ini_values, fn = log_likelihood_flex_spline, 
               data = data, knots = knots, cov = var,
               control = list(maxit = 2000))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(hes){
    hessian <- solve(numDeriv::hessian(log_likelihood_flex_spline, 
                                       res$par, knots = knots, data = data, cov = var))
    list(Data = data, var = var, Lambda1 = res$par[1:length(knots)], 
         Lambda2 = res$par[(length(knots) + 1):length(res$par)], 
         knots = knots, ML = res$value, Hessian = hessian) 
  }else{
    list(Data = data, var = var, Lambda1 = res$par[1:length(knots)], 
         Lambda2 = res$par[(length(knots) + 1):length(res$par)], 
         knots = knots, ML = res$value) 
  }
}

fit <- FlexrelsurvModel_spline("age_years", DLBCL)

plot.FlexrelsurvModel_spline <- function(fit, z, add = F, col = 1){
  times <- seq(0, max(fit$Data$FU), length.out = 1000)
  f <- function(t) flex_spline_surv(t = t / year, lambda1 = fit$Lambda1, 
                                            lambda2 = fit$Lambda2, knots = fit$knots, z = z)
  if(add){
    lines(f(times) ~ times, col = col)
  }else{
    plot(f(times) ~ times, ylim = c(0, 1), type = "l", col = col) 
  }
}

plot.FlexrelsurvModel_spline(fit, 50)
plot.FlexrelsurvModel_spline(fit, 55, add = T, col = 2)
plot.FlexrelsurvModel_spline(fit, 60, add = T, col = 3)
plot.FlexrelsurvModel_spline(fit, 65, add = T, col = 4)


res1 <- calc_LOL_rel_flex_extra_cov(fit, z = 50, newdata = data.frame(age = 50 * year, sex = "male", year = 2007))
res2 <- calc_LOL_rel_flex_extra_cov(fit, z = 55, newdata = data.frame(age = 55 * year, sex = "male", year = 2007))
res3 <- calc_LOL_rel_flex_extra_cov(fit, z = 60, newdata = data.frame(age = 60 * year, sex = "male", year = 2007))
res4 <- calc_LOL_rel_flex_extra_cov(fit, z = 65, newdata = data.frame(age = 65 * year, sex = "male", year = 2007))
res5 <- calc_LOL_rel_flex_extra_cov(fit, z = 70, newdata = data.frame(age = 70 * year, sex = "male", year = 2007))

res <- rbind(res1, res2, res3, res4, res5)
res$age <- factor(rep(c(50, 55, 60, 65, 70), each = nrow(res1)))
ggplot(res, aes(x = Time, y = LossOfLifetime, colour = age, group = age)) + geom_line()


ages <- seq(45, 80, by = 1)
res <- lapply(ages, function(x){
  calc_LOL_rel_flex_extra_cov(fit, z = x, 
                              newdata = data.frame(age = x * year, sex = "male", year = 2007), 
                                                   time = seq(0, 5, by = 1))
})

LOL_base <- sapply(res, function(x) x[x$Time == 0, "LossOfLifetime"])
LOL_5 <- sapply(res, function(x) x[x$Time == 2, "LossOfLifetime"])

D <- data.frame(LOL = c(LOL_base, LOL_5), Time = factor(rep(c(0, 2), each = length(ages))), Age = rep(ages, 2))
ggplot(D, aes(x = Age, y = LOL, linetype = Time, group = Time)) + geom_line()


fit <- FlexrelsurvModel_spline("age_years", FL)
ages <- seq(45, 80, by = 1)
res <- lapply(list(DLBCL, FL, BL), function(y){
  cat("Og det var eeeeeennnn....\n")
  fit <- FlexrelsurvModel_spline("age_years", y)
  sapply(ages, function(x){
    calc_stat_cure_flex_extra_cov(fit, 2, z = x, newdata = data.frame(age = x * year, sex = "male", year = 2007))
  }) 
})

res_new <- data.frame(stat_cure = unlist(res), Age = rep(ages, 3))
res_new$disease <- rep(c("DLBCL", "FL", "BL"), each = length(ages))

p <- ggplot(res_new, aes(x = Age, y = stat_cure, group = disease, colour = disease)) + 
  geom_line() + ylab("Time to statistical cure (years)") + xlab("Age at diagnosis") + 
  theme_bw() + theme(legend.justification = c(1,1), legend.position = c(0.95,0.95),
                     legend.key = element_blank(), legend.key.width=unit(1,"cm")) + 
  scale_colour_discrete(name = NULL)

pdf(file.path(fig.out, "TOSC_LOL_disease.pdf"), width = 7, height = 5)
print(p)
dev.off()
  