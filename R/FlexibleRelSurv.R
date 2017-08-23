
fit.flex <- function(formula, data, knots = NULL, n.knots, last_knot = NULL, bhazard,
                     model = "relsurv", tvc.formula = NULL, knots.time = NULL,
                     n.knots.time = NULL, add.knots = NULL){

  fu <- all.vars(formula)[1]
  status <- all.vars(formula)[2]

  # Find knots
  ti <- data[data[, status] == 1, fu]
  if(is.null(knots)){
    bdr_knots <- range(ti)
    inner_knots <- quantile(ti, 1 / (n.knots - 1)*1:(n.knots - 2))
    knots <- log(sort(c(bdr_knots, inner_knots)))
  }

  if(!is.null(tvc.formula) & is.null(knots.time)){
    term.names <- colnames(attr(terms(tvc.formula), "factors"))
    n.knots.time <- n.knots.time[term.names]
    knots.time <- lapply(n.knots.time, function(x){
      bdr_knots <- range(ti)
      if(x > 2){
        inner_knots <- quantile(ti, 1 / (x - 1) * 1:(x - 2))
        log(sort(c(bdr_knots, inner_knots)))
      }else{
        bdr_knots
      }
    })
  }

  # Get the base functions
  if(model == "cure"){
    base_function <- basis_cure
    dbase_function <- dbasis_cure
  }else if(model == "relsurv"){
    base_function <- basis
    dbase_function <- dbasis
  }

  # Rewrite formula
  formula.2 <- update(formula, ~ . -1)
  smooth.formula <- as.formula(paste0("~base_function(knots = knots, log(", fu, "))"))


  if(!is.null(tvc.formula)){
    b <- paste0("base_function(log(FU_years), knots = knots.time[[", 1:length(knots.time), "]])")
    tvc.formula.2 <- paste0("~ ", paste0(term.names, ":", b, collapse = " + "))
    tvc.formula.2 <- as.formula(tvc.formula.2)
    fit <- stpm2(formula.2, data = data,
                 bhazard = bhazard,
                 smooth.formula = smooth.formula,
                 tvc.formula = tvc.formula.2)
  }else{
    fit <- stpm2(formula.2, data = data,
                 bhazard = bhazard,
                 smooth.formula = smooth.formula)
  }

  L <- list(Data = data, base_function = base_function, dbase_function = dbase_function,
            knots = knots, knots.time = knots.time, coefs = fit@coef, covariance = fit@vcov,
            formula = formula, tvc.formula = tvc.formula, model = model, sum = summary(fit))
  class(L) <- "FlexSurvModel"
  L
}

#Debugging of the stpm2 code

print.FlexSurvModel <- function(fit){
  cat("Model:\n")
  if(fit$model == "relsurv"){
    cat("Parametric relative survival model\n")
  }else{
    cat("Parametric cure model\n")
  }
  cat("\nCall:\n")
  print(fit$formula)
  if(!is.null(fit$tvc.formula)){
    cat("\nCall.tvc:\n")
    print(fit$tvc.formula)
  }
  cat("\nCoefficients:\n")
  print(fit$coefs)
}

summary.FlexSurvModel <- function(fit){
  sum <- fit$sum
  coefs <- sum@coef
  rownames(coefs) <- gsub("base_function\\(knots = knots, log\\(FU_years\\)\\)", "spline", rownames(coefs))
  for(i in 1:length(fit$knots.time)){
    replace_pattern <- paste0("base_function\\(log\\(FU_years\\), knots = knots.time\\[\\[", i, "\\]\\]\\)")
    rownames(coefs) <- gsub(replace_pattern, "spline", rownames(coefs))
  }

  L <- list(coefs = coefs, ML = sum@m2logL, formula = fit$formula, tvc.formula  = fit$tvc.formula,
            cov = fit$covariance, fit$knots, fit$knots.time, model = fit$model)
  class(L) <- "summary.FlexSurvModel"
  L
}

print.summary.FlexSurvModel <- function(x){
  cat("Model:\n")
  if(x$model == "relsurv"){
    cat("Parametric relative survival model\n")
  }else{
    cat("Parametric cure model\n")
  }
  cat("\nCall:\n")
  print(x$formula)
  if(!is.null(x$tvc.formula)){
    cat("\nCall.tvc:\n")
    print(x$tvc.formula)
  }
  cat("\nCoefficients:\n")
  printCoefmat(x$coefs, P.value = TRUE, has.Pvalue = TRUE)
}

predict.FlexSurvModel <- function(fit, newdata, times, type = "relsurv"){
  if(type == "relsurv"){
    if(is.null(newdata)){
      lps <- t(matrix(unlist(fit$coefs[-1])))
      rs <- exp(-exp(fit$base_function(knots = fit$knots, log(times)) %*% fit$coefs))
      rss <- data.frame(times = times, RS = rs)
      rss
    }else{
      baseline <- fit$base_function(knots = fit$knots, log(times))
      tvc.basis <- lapply(1:length(fit$knots.time), function(x){
        fit$base_function(knots = fit$knots.time[[x]], log(times))
      })
      term.names <- colnames(attr(terms(fit$formula), "factors"))
      lp <- as.matrix(newdata[, term.names])

      splines <- matrix(ncol = nrow(newdata), nrow = length(times))
      tvc.term.names <- names(fit$knots.time)
      for(i in 1:nrow(newdata)){
        tvc.basis.i <- lapply(1:length(tvc.basis), function(x){
          tvc.basis[[x]] * newdata[i, tvc.term.names]
        })
        splines[,i] <- cbind(lp[i,], baseline, do.call(cbind, tvc.basis.i)) %*% fit$coefs
      }
      rss <- exp(-exp(splines))
      colnames(rss) <- paste0("RS", 1:ncol(rss))
      rss[times == 0,] <- 1
      rss <- cbind(times, rss)
      rss
    }
  }else if(type == "curerate"){
    if(is.null(newdata)){
      exp(-exp(fit$coefs[1]))
    }else{
      tt <- terms(fit$formula)
      formula.2 <- formula(delete.response(tt))
      design_matrices <- as.matrix(model.matrix(formula.2, data = newdata)[,-1])
      lps <- design_matrices %*% fit$coefs[[1]]
      pi <- exp(-exp(lps + fit$coefs[grepl("base_function", names(fit$coefs))][1]))
      pi
    }
  }else if(type == "probcure"){
    if(is.null(newdata)){
      pi <- exp(-exp(fit$coefs[1]))
      rs <- ifelse(times == 0, 1, exp(-exp(fit$base_function(knots = fit$knots, x = log(times)) %*% fit$coefs)))
      probcure <- data.frame(times = times, probtime = pi / rs)
      probcure
    }else{
      baseline <- fit$base_function(knots = fit$knots, log(times))
      tvc.basis <- lapply(1:length(fit$knots.time), function(x){
        fit$base_function(knots = fit$knots.time[[x]], log(times))
      })
      term.names <- colnames(attr(terms(fit$formula), "factors"))
      lp <- as.matrix(newdata[, term.names])

      prob_cure <- matrix(ncol = nrow(newdata), nrow = length(times))
      tvc.term.names <- names(fit$knots.time)
      for(i in 1:nrow(newdata)){
        tvc.basis.i <- lapply(1:length(tvc.basis), function(x){
          tvc.basis[[x]] * newdata[i, tvc.term.names]
        })
        spline_vars <- grepl("base_function", names(fit$coefs))
        pi <- exp(-exp(fit$coefs[spline_vars][1] + lp[i,] %*% fit$coefs[!spline_vars]))
        splines <- cbind(lp[i,], baseline, do.call(cbind, tvc.basis.i)) %*% fit$coefs
        prob_cure[,i] <- pi / ifelse(times == 0, 1, exp(-exp(splines)))
      }
      rss <- exp(-exp(splines))


      rss[times == 0] <- 1
      colnames(rss) <- paste0("RS", 1:ncol(rss))
      rss <- cbind(times, rss)
      rss


      tt <- terms(fit$formulas$formula.gamma)
      formula.2 <- formula(delete.response(tt))
      fit$formulas$formula.gamma <- formula.2
      design_matrices <- lapply(fit$formulas, function(x){
        if(length(x) > 0){
          return(model.matrix(x, data = newdata))
        }
      })

      lps <- vector("list", length(design_matrices))
      for(i in 1:length(lps)){
        if(length(design_matrices[[i]]) > 0){
          lps[[i]] <- design_matrices[[i]] %*% fit$coefs[[i]]
        }
      }

      #lps <- do.call(cbind, lps)
      probtime <- matrix(nrow = length(times), ncol = nrow(newdata))
      for(i in 1:nrow(newdata)){
        pi <- link_fun(lps[[1]][i,])
        probtime[, i] <- pi / (pi + (1 - pi) * surv_fun(times, lp.k1 = lps[[2]][i,], lp.k2 = lps[[3]][i,], lp.k3 = lps[[4]][i,]))
      }
      colnames(probtime) <- paste0("CP", 1:ncol(probtime))
      probcure <- cbind(times, probtime)
      probcure
    }
  }
}
