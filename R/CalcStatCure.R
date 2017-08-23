quantile.calc.Crude <- function(fit, q = 0.95, newdata = NULL){
  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    f <- function(time, q) calc.Crude(fit, time = time, type = "othertime", ci = F,
                                      newdata = newdata[i,,drop = F])$prob[[1]]$prob - q
    uni <- uniroot.all(f, lower = 0, upper = 20, q = q)
    gr <- grad(f, x = uni, q = 0)
    VAR <- calc.Crude(fit, time = uni, type = "othertime")$prob[[1]]$var
    VAR2 <- gr ^ (-2) * VAR
    upper <- uni + sqrt(VAR2) * qnorm(0.975)
    lower <- uni - sqrt(VAR2) * qnorm(0.975)
    data.frame(Est = uni, var = VAR2, lower.ci = lower, upper.ci = upper)
  })

  do.call(rbind, ests)
}


quantile.calc.LL <- function(fit, q = 2, newdata = NULL){
  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    f <- function(time, q) calc.LL(fit, time = time, ci = F,
                                   newdata = newdata[i,,drop = F])$LOL[[1]]$LL - q
    uni <- uniroot.all(f, lower = 0, upper = 30, q = q)
    gr <- grad(f, x = uni, q = 0)
    VAR <- calc.LL(fit, time = uni)$LOL[[1]]$Var
    VAR2 <- gr ^ (-2) * VAR
    upper <- uni + sqrt(VAR2) * qnorm(0.975)
    lower <- uni - sqrt(VAR2) * qnorm(0.975)
    data.frame(Est = uni, var = VAR2, lower.ci = lower, upper.ci = upper)
  })
  do.call(rbind, ests)
}

