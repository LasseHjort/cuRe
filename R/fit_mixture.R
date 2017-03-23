#Mixture cure models

ParaCureModel <- function(formula, data, ratetable = survexp.dk, dist = "weibull", model = "mixture",
                             bhazard = NULL, hes = FALSE){
  data <- na.omit(data, cols = all.vars(formula))
  data$FU <- data[, "OS"][,1] * year
  if(dist == "weibull"){
    surv <- wei_surv
    dens <- wei_dens
  }

  X <- model.matrix(formula, data = data)
  ini_values <- initialParaCure(formula, data, X, surv, dens, model = model)
  response <- all.vars(formula)[1]
  res <- optim(ini_values,
               fn = minusloglike_parametricCure,
               X = X,
               status = data[, response][, 2],
               futime = data[, response][, 1],
               bhazard = data$exp_haz,
               model = model,
               surv = surv,
               dens = dens,
               control = list(maxit = 2000))

  if(res$convergence != 0){
    warning("Convergence not reached")
  }

  if(hes){
    hessian <- solve(numDeriv::hessian(log_likelihood_mixtureCure, res$par, X = X,
                                       data = data, surv = surv, dens = dens))
    list(Data = data, formula = formula, dist = dist, model = model, gamma = res$par[1:ncol(X)],
         Beta = res$par[(ncol(X) + 1): (ncol(X)*2)], shape = res$par[length(res$par)],
         ML = res$value, Hessian = hessian, surv = surv)
  }else{
    list(Data = data, formula = formula, dist = dist, model = model, gamma = res$par[1:ncol(X)],
         Beta = res$par[(ncol(X) + 1): (ncol(X)*2)], shape = res$par[length(res$par)],
         ML = res$value, surv = surv)
    }
}

getCure.ParaCureModel <- function(fit, newdata = NULL){
  if(!is.null(newdata)){
    if(!"(Intercept)" %in% names(newdata)){
      newdata <- cbind("(Intercept)" = 1, newdata)
    }
  }else{
    newdata <- data.frame(1)
  }
  lp <- as.matrix(newdata) %*% fit$gamma
  as.numeric(exp(lp) / (exp(lp) + 1))
}

plot.ParaCureModel <- function(fit, newdata = NULL, from = 0, to = max(fit$Data[, all.vars(formula)[1]][,1]),
                        ylab = "Relative survival", xlab = "Follow-up time (years)", n = 100, add.ederer = F){
  if(!is.null(newdata)){
    if(!"(Intercept)" %in% names(newdata)){
      newdata <- cbind("(Intercept)" = 1, newdata)
    }
  }else{
    newdata <- data.frame("(Intercept)" = 1)
  }
  rel_survival_function <- function(t){
    scale <- exp(as.matrix(newdata) %*% fit$Beta)
    lp <- as.matrix(newdata) %*% fit$gamma
    pi <- as.numeric(exp(lp) / (exp(lp) + 1))
    shape <- exp(fit$shape)
    sapply(1:nrow(newdata), function(x){
      if(fit$model == "mixture"){
        pi[x] + (1 - pi[x]) * fit$surv(t, scale = scale[x], shape = shape)
      }else{
        pi[x] ^ (1 - fit$surv(t, scale = scale[x], shape = shape))
      }
      })
  }
  times <- seq(from, to, length.out = n)
  D <- rel_survival_function(times)
  if(add.ederer){
    rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                     data = fit$Data, ratetable = survexp.dk, method = "ederer2")
    plot(rsfit, conf.int = T, xscale = year, ylab = "Relative survival",
         xlab = "Follow-up time (years)")
  }

  for(i in 1:ncol(D)){
    if(i == 1){
      if(add.ederer){
        lines(times * year, D[, i], col = i + 1, type = "l")
      }else{
        plot(times, D[, i], col = i, type = "l", ylab = ylab, xlab = xlab, ylim = c(0, 1))
      }
    }else{
      lines(times, D[, i], col = i)
    }
  }
}

curePlot.ParaCureModel <- function(fit, newdata = NULL, from = 0, to = max(fit$Data[, all.vars(formula)[1]][,1]),
                               ylab = "Probability of cure", xlab = "Follow-up time (years)", col = NULL, n = 100){
  if(!is.null(newdata)){
    if(!"(Intercept)" %in% names(newdata)){
      newdata <- cbind("(Intercept)" = 1, newdata)
    }
  }else{
    newdata <- data.frame(1)
  }

  cure_prop <- function(t){
    scale <- exp(as.matrix(newdata) %*% fit$Beta)
    lp <- as.matrix(newdata) %*% fit$gamma
    pi <- as.numeric(exp(lp) / (exp(lp) + 1))
    shape <- exp(fit$shape)
    sapply(1:nrow(newdata), function(x){
      if(fit$model == "mixture"){
        pi[x] / (pi[x] + (1 - pi[x]) * fit$surv(t, scale = scale[x], shape = shape))
      }else{
        pi[x] / (pi[x] ^ (1 - fit$surv(t, scale = scale[x], shape = shape)))
      }
    })
  }
  times <- seq(from, to, length.out = n)
  D <- cure_prop(times)
  if(is.null(col)){
    col <- 1:nrow(newdata)
  }

  for(i in 1:ncol(D)){
    if(i == 1){
      plot(times, D[, i], col = col[i], type = "l", ylab = ylab, xlab = xlab, ylim = c(0, 1))
    }else{
      lines(times, D[, i], col = col[i])
    }
  }
}

evalCure.ParaCureModel <- function(fit, newdata = NULL, times){
  if(!is.null(newdata)){
    if(!"(Intercept)" %in% names(newdata)){
      newdata <- cbind("(Intercept)" = 1, newdata)
    }
  }else{
    newdata <- data.frame(1)
  }

  cure_prop <- function(t){
    scale <- exp(as.matrix(newdata) %*% fit$Beta)
    lp <- as.matrix(newdata) %*% fit$gamma
    pi <- as.numeric(exp(lp) / (exp(lp) + 1))
    shape <- exp(fit$shape)
    sapply(1:nrow(newdata), function(x){
      if(fit$model == "mixture"){
        pi[x] / (pi[x] + (1 - pi[x]) * fit$surv(t, scale = scale[x], shape = shape))
      }else{
        pi[x] / (pi[x] ^ (1 - fit$surv(t, scale = scale[x], shape = shape)))
      }
    })
  }
  cbind(data.frame(Time = times), cure_prop(times))
}


fit <- ParaCureModel(OS ~ age_years, data = DLBCL)
newdata <- data.frame(age_years = c(60, 70))
getCure.ParaCureModel(fit, newdata = newdata)
plot.ParaCureModel(fit, newdata = newdata)
curePlot.ParaCureModel(fit, newdata = newdata)
evalCure.ParaCureModel(fit, newdata = newdata, times = c(0, 1, 2, 3))

fit <- ParaCureModel(OS ~ 1, data = DLBCL)
plot.ParaCureModel(fit, add.ederer = T)


fit <- ParaCureModel(OS ~ 1, data = DLBCL, model = "nonmixture")
plot.ParaCureModel(fit, add.ederer = T)






  if(fit$var != "1"){
    res <- lapply(levels(fit$Data[, fit$var]), function(x){
      cat("Calculating for", fit$var, "=", x, "\n")
      data.frame(CureProb = cure_prob(times, var_val = x), Time = times)
    })
    names(res) <- levels(fit$Data[, fit$var])
  }else{
    res <- data.frame(CureProb = cure_prob(times), Time = times)
    res <- list(res)
    names(res) <- "1"
  }

  res <- lapply(1:length(res), function(x){
    res[[x]]$Group <- names(res)[x]
    res[[x]]
  })
  D <- do.call(rbind, res)

  p <- ggplot(D, aes(x = Time, y = CureProb, group = Group, colour = Group)) + geom_line() +
    xlab("Follow-up time (years)") + ylab("Probability of being cured") +
    theme_bw() + theme(legend.justification = c(1,1), legend.position = c(1,1),
                       legend.key = element_blank(), legend.key.width=unit(1,"cm")) +
    geom_hline(yintercept = 0.95, linetype = "dashed")

  if(length(unique(D$Group)) == 1){
    p <- p + scale_colour_manual(guide = F, values = "black")
  }else{
    p <- p + scale_colour_discrete(name = var_name)
  }
  p
}

fit <- MixtureCureModel("1", data = surv_data)
cure_plots_S <- function(fit, tau, n = 100, var_name = fit$var, add_lines = FALSE){
  times <- seq(0, tau, length.out = 100)
  cure_prob <- function(t, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    pi <- exp(X %*% fit$Gamma) / (exp(X %*% fit$Gamma) + 1)
    pi / (pi + (1 - pi) * exp(-scale * t ^ shape))
  }

  para_model <- function(t, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    exp(-scale * t ^ shape)
  }

  if(fit$var != "1"){
    res <- lapply(levels(fit$Data[, fit$var]), function(x){
      cat("Calculating for", fit$var, "=", x, "\n")
      data.frame(CureProb = cure_prob(times, var_val = x), Time = times)
    })
    names(res) <- levels(fit$Data[, fit$var])
  }else{
    res <- data.frame(CureProb = c(cure_prob(times), para_model(times)), Time = rep(times, 2),
                      type = rep(c("Prob of cure", "Surv of uncured"), each = length(times)))
    res <- list(res)
    names(res) <- "1"
  }

  res <- lapply(1:length(res), function(x){
    res[[x]]$Group <- names(res)[x]
    res[[x]]
  })
  D <- do.call(rbind, res)

  p <- ggplot(D, aes(x = Time, y = CureProb, group = type, colour = type)) + geom_line() +
    xlab("Follow-up time (years)") + ylab("Probability") + ylim(c(0,1)) +
    theme_bw() + theme(legend.justification = c(1,1), legend.position = c(0.3,0.3),
                       legend.key = element_blank(), legend.key.width=unit(1,"cm")) +
    geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed") + scale_colour_discrete(name = var_name)

  if(add_lines){
    f <- function(t) para_model(t) - 0.05
    uni <- uniroot(f, interval = c(0, 100))$root
    f <- function(t) cure_prob(t) - 0.95
    uni2 <- uniroot(f, interval = c(0, 100))$root

    p <- p + geom_vline(xintercept = c(uni, uni2), linetype = "dashed")
  }
  p
}

fit <- MixtureCureModel("1", data = surv_data)
cure_plots_S(fit, tau = 13, var_name = NULL)
cure_plots_S(fit, tau = 13, var_name = NULL, add_lines = T)

fit <- MixtureCureModel("age_above_50", data = surv_data)
cure_plots_S2 <- function(fit, tau, n = 100, var_name = fit$var, add_lines = FALSE){
  times <- seq(0, tau, length.out = 100)
  cure_prob <- function(t, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    pi <- exp(X %*% fit$Gamma) / (exp(X %*% fit$Gamma) + 1)
    pi / (pi + (1 - pi) * exp(-scale * t ^ shape))
  }

  para_model <- function(t, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    exp(-scale * t ^ shape)
  }

  if(fit$var != "1"){
    res <- lapply(levels(fit$Data[, fit$var]), function(x){
      cat("Calculating for", fit$var, "=", x, "\n")
      res <- data.frame(CureProb = c(cure_prob(times, var_val = x), para_model(times, var_val = x)),
                        Time = rep(times, 2),
                        type = rep(c("Prob of cure", "Surv of uncured"), each = length(times)))
    })
    names(res) <- levels(fit$Data[, fit$var])
  }else{
    res <- data.frame(CureProb = c(cure_prob(times), para_model(times)), Time = rep(times, 2),
                      type = rep(c("Prob of cure", "Surv of uncured"), each = length(times)))
    res <- list(res)
    names(res) <- "1"
  }

  res <- lapply(1:length(res), function(x){
    res[[x]]$Group <- names(res)[x]
    res[[x]]
  })
  D <- do.call(rbind, res)
  D$Group <- factor(D$Group)

  p <- ggplot(D, aes(x = Time, y = CureProb, group = Group:type, colour = Group, linetype = type)) + geom_line() +
    xlab("Follow-up time (years)") + ylab("Probability") + ylim(c(0,1)) +
    theme_bw() + theme(legend.justification = c(1,1), legend.position = c(0.25,0.4),
                       legend.key = element_blank(), legend.key.width=unit(1,"cm")) +
    geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed") + scale_linetype_discrete(name = NULL)

  if(add_lines){
    f <- function(t) para_model(t) - 0.05
    uni <- uniroot(f, interval = c(0, 100))$root
    f <- function(t) cure_prob(t) - 0.95
    uni2 <- uniroot(f, interval = c(0, 100))$root

    p <- p + geom_vline(xintercept = c(uni, uni2), linetype = "dashed")
  }
  p
}
cure_plots_S2(fit, tau = 13)
fit <- MixtureCureModel("1", data = surv_data)
cure_plots_S2(fit, tau = 13)

cure_plots(fit, tau = 12)
undebug(cure_plots)

fit <- MixtureCureModel("IPI", data = surv_data)
get_pis(fit)

fit <- MixtureCureModel("sex", data = surv_data)
get_pis(fit)


fit <- MixtureCureModel("1", data = surv_data[surv_data$IPI == "1",], n = 20)
fit <- MixtureCureModel("1", data = surv_data[surv_data$IPI == "1",], n = 20)
get_pis(fit)
plot_curves(fit)
cure_plots(fit, tau = 12)


get_data <- function(fit, tau, n = 100, var_name = fit$var){
  times <- seq(0, tau, length.out = 100)
  cure_prob <- function(t, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    pi <- exp(X %*% fit$Gamma) / (exp(X %*% fit$Gamma) + 1)
    pi / (pi + (1 - pi) * exp(-scale * t ^ shape))
  }

  if(fit$var != "1"){
    res <- lapply(levels(data[, fit$var]), function(x){
      cat("Calculating for", fit$var, "=", x, "\n")
      data.frame(CureProb = cure_prob(times, var_val = x), Time = times)
    })
    names(res) <- levels(data[, fit$var])
  }else{
    res <- data.frame(CureProb = cure_prob(times), Time = times)
    res <- list(res)
    names(res) <- "1"
  }

  res <- lapply(1:length(res), function(x){
    res[[x]]$Group <- names(res)[x]
    res[[x]]
  })
  do.call(rbind, res)
}


fit1 <- MixtureCureModel("1", data = surv_data, n = 20)
a <- get_data(fit1, tau = 12)
fit2 <- nonMixtureCureModel("1", data = surv_data, n = 20)
b <- get_data2(fit2, tau = 12)

a <- rbind(a, b)
a$Group <- rep(c("Mix", "NonMix"), each = nrow(b))

ggplot(a, aes(x = Time, y = CureProb, group = Group, colour = Group)) + geom_line() +
  xlab("Follow-up time (years)") + ylab("Probability of being cured") +
  theme_bw() + theme(legend.justification = c(1,1), legend.position = c(1,1),
                     legend.key = element_blank(), legend.key.width=unit(1,"cm")) +
  geom_hline(yintercept = 0.95, linetype = "dashed")


fit1 <- MixtureCureModel("IPI", data = surv_data)
fit2 <- MixtureCureModel_old("IPI", data = surv_data, n = 1000)
get_pis(fit1)
get_pis(fit2)
fit1$ML
fit2$ML

surv_data$sex <- factor(surv_data$sex)
fit <- MixtureCureModel("1", data = surv_data)
cure_plots(fit, tau = 12)


fit <- MixtureCureModel("1", data = surv_data)
find_cure_time <- function(fit, alpha = 0.95){
  cure_point <- function(var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    pi <- exp(X %*% fit$Gamma) / (exp(X %*% fit$Gamma) + 1)
    est <- (-log((pi - alpha * pi) / (alpha * (1 - pi))) / scale) ^ (1/shape)
    log_term <- log((pi / alpha - pi) / (1 - pi))
    pi_deriv <- scale^(1 / shape) * (-log_term ) ^ (1 / shape - 1) * 1/shape * (1 / alpha - 1) / ((1 - pi) * (pi/alpha - pi))
    scale_deriv <- -(-log_term) ^ (1 / shape) * 1/shape * scale ^ (-1/shape - 1)
    shape_deriv <- -est * (-log_term / scale) ^ (1 / shape) * 1/shape ^ 2
    jac <- matrix(c(scale_deriv, pi_deriv, shape_deriv), ncol = 1)
    VAR <- t(jac) %*% solve(-fit$Hessian) %*% jac / nrow(fit$Data)
    data.frame("Stat Cure" = est, SE = sqrt(VAR),
               Lower95 = est - qnorm(0.975) * sqrt(VAR),
               Upper95 = est + qnorm(0.975) * sqrt(VAR))
  }
  if(fit$var != "1"){
    sapply(levels(fit$Data[, fit$var]), cure_point)
  }else{
    cure_point()
  }
}


find_cure_time_boot <- function(fit, alpha = 0.95, m = 200){
  cure_point <- function(fit, var_val = NULL){
    if(fit$var != "1"){
      wh <- which(levels(fit$Data[, fit$var]) == var_val)
      X <- c(1, rep(0, nlevels(fit$Data[, fit$var]) - 1))
      X[wh] <- 1
    }else{
      X <- 1
    }
    scale <- exp(X %*% fit$Beta)
    shape <- exp(fit$Shape)
    pi <- exp(X %*% fit$Gamma) / (exp(X %*% fit$Gamma) + 1)
    (-log((pi - alpha * pi) / (alpha * (1 - pi))) / scale) ^ (1/shape)
  }
  if(fit$var != "1"){
    cure_time <- matrix(nrow = m, ncol = nlevels(fit$Data[, var]))
  }else{
    cure_time <- rep(NA, m)
  }

  for(i in 1:m){
    cat(i, "\n")
    boot <- sample(1:nrow(fit$Data), replace = T)
    fit_tmp <- MixtureCureModel(fit$var, data = fit$Data[boot,])
    get_pis(fit_tmp)
    plot_curves(fit_tmp)
    cure_plots(fit_tmp, tau = 2000)
    undebug(cure_plots)
    if(fit$var != "1"){
      cure_time[i,] <- sapply(levels(fit$Data[, fit$var]), cure_point)
    }else{
      cure_time[i] <- cure_point(fit_tmp)
    }
  }

  if(fit$var != "1"){
    est <- sapply(levels(fit$Data[, fit$var]), cure_point)
    VAR <- apply(cure_time, 2, var, na.rm = TRUE)
    D <- data.frame("Stat Cure" = est, SE = sqrt(VAR),
                    Lower95 = est - qnorm(0.975) * sqrt(VAR),
                    Upper95 = est + qnorm(0.975) * sqrt(VAR))
  }else{
    est <- cure_point(fit)
    VAR <- var(cure_time, na.rm = TRUE)
    D <- c("Stat Cure" = est, SE = sqrt(VAR),
           Lower95 = est - qnorm(0.975) * sqrt(VAR),
           Upper95 = est + qnorm(0.975) * sqrt(VAR))
  }
  D
}

find_cure_time(fit)

