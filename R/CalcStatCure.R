quantile.calc.Crude <- function(fit, q = 0.95, newdata = NULL, max.time = 20, expected = NULL, ci = TRUE,
                                rmap, ratetable = survexp.dk, last.point = 100, reverse = FALSE){

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, last.point + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
      }else{
        data <- fit$data
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- lapply(1:nrow(newdata), function(x){
        do.call("survexp",
                list(formula = ~ 1, rmap = substitute(rmap),
                     data = newdata[x, ], ratetable = ratetable,
                     scale = ayear, times = times * ayear))
      })
    }
  }


  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    f <- function(time, q) calc.Crude(fit, time = time, type = "othertime",
                                      ci = F, newdata = newdata[i,,drop = F],
                                      expected = expected[i], reverse = reverse)$prob[[1]]$prob - q
    uni <- rootSolve::uniroot.all(f, lower = 0, upper = max.time, q = q)
    if(ci){
      gr <- grad(f, x = uni, q = 0)
      VAR <- calc.Crude(fit, time = uni, expected = expected[i], newdata = newdata[i,,drop = F],
                        type = "othertime", link = "identity", reverse = reverse)$prob[[1]]$var
      VAR2 <- gr ^ (-2) * VAR
      upper <- uni + sqrt(VAR2) * qnorm(0.975)
      lower <- uni - sqrt(VAR2) * qnorm(0.975)
      data.frame(Est = uni, var = VAR2, lower.ci = lower, upper.ci = upper)
    } else{
      data.frame(Est = uni)
    }
  })

  do.call(rbind, ests)
}


quantile.calc.LL <- function(fit, q = 2, newdata = NULL, max.time = 20, ci = TRUE,
                             expected = NULL, rmap, ratetable = survexp.dk,
                             last.point = 100, type = "ll"){

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, last.point + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
      }else{
        data <- fit$data
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- lapply(1:nrow(newdata), function(x){
        do.call("survexp",
                list(formula = ~ 1, rmap = substitute(rmap),
                     data = newdata[x, ], ratetable = ratetable,
                     scale = ayear, times = times * ayear))
      })
    }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){
    f <- function(time, q) calc.LL(fit, time = time, ci = F, expected = expected[i],
                                   newdata = newdata[i,,drop = F])$Ests[[1]][, type] - q
    uni <- rootSolve::uniroot.all(f, lower = 0, upper = max.time, q = q)
    if(ci){
      gr <- grad(f, x = uni, q = 0)
      VAR <- calc.LL(fit, time = uni, expected = expected[i],
                     newdata = newdata[i,, drop = F])$Ests[[1]]$Var
      VAR2 <- gr ^ (-2) * VAR
      upper <- uni + sqrt(VAR2) * qnorm(0.975)
      lower <- uni - sqrt(VAR2) * qnorm(0.975)
      data.frame(Est = uni, var = VAR2, lower.ci = lower, upper.ci = upper)
    } else {
      data.frame(Est = uni)
    }
  })
  do.call(rbind, ests)
}




quantile.calc.LL2 <- function(fit, q = 2, newdata = NULL, max.time = 20, ci = TRUE,
                             expected = NULL, rmap, ratetable = survexp.dk,
                             last.point = 100, type = "ll"){

  if(is.null(expected)){
    #The time points for the expected survival
    times <- seq(0, last.point + 1, by = 0.05)

    #Extract expected survival function
    if(is.null(newdata)){
      if(any(class(fit) %in% c("stpm2", "pstpm2"))){
        data <- fit@data
        newdata <- data.frame(arbritary_var = 0)
        coefs <- fit@coef
        covar <- fit@vcov
      }else{
        data <- fit$data
        coefs <- c(fit$coefs, fit$coefs.spline)
        covar <- fit$covariance
      }
      expected <- list(do.call("survexp",
                               list(formula = ~ 1, rmap = substitute(rmap),
                                    data = data, ratetable = ratetable,
                                    scale = ayear, times = times * ayear)))
    }else{
      expected <- lapply(1:nrow(newdata), function(x){
        do.call("survexp",
                list(formula = ~ 1, rmap = substitute(rmap),
                     data = newdata[x, ], ratetable = ratetable,
                     scale = ayear, times = times * ayear))
      })
    }
  }

  n.obs <- ifelse(is.null(newdata), 1, nrow(newdata))
  ests <- lapply(1:n.obs, function(i){

    g <- function(pars, q){
      f <- function(time) calc.LL(fit = fit, time = time, ci = F, expected = expected[i],
                                  newdata = newdata[i,, drop = F], pars = pars)$Ests[[1]][, type] - q
      rootSolve::uniroot.all(f, lower = 0, upper = max.time)
    }

    uni <- g(pars = NULL, q = q)
    if(ci){
      gr <- grad(g, x = coefs, q = q)
      VAR <- gr %*% covar %*% gr
      upper <- uni + sqrt(VAR) * qnorm(0.975)
      lower <- uni - sqrt(VAR) * qnorm(0.975)
      data.frame(Est = uni, var = VAR, lower.ci = lower, upper.ci = upper)
    } else {
      data.frame(Est = uni)
    }
  })
  do.call(rbind, ests)
}

