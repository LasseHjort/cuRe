

penalization.fun <- function(beta,
                             alpha,
                             lambda) lambda * (alpha * sum(beta ^ 2) + (1 - alpha) * sum(abs(beta)))

penalizedFlexRelsurvModel <- function(formula, data, bhazard,
                                      knots = NULL, n.knots = NULL,
                                      message = T, cure = F, alpha = 1, lambda = 0.1){
  cc <- complete.cases(data[,all.vars(formula)])
  data.2 <- data[cc, ]
  #Extract relevant variables
  fu <- eval(formula[[2]][[2]], envir = data.2)
  status <- eval(formula[[2]][[3]], envir = data.2)
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

  X_time <- model.matrix(tvc.formula, data = data.2)[,-1, drop = FALSE]
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
  X <- model.matrix(formula, data = data.2)
  b <- cbind(b, X[, -1], b_time)
  db <- cbind(db, matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)), db_time)

  #Generate initial values by running a mixture weibull cure model
  if(message) cat("Finding initial values...")
  #formula(as.Formula(terms(formula),formula(Formula(terms(tvc.formula)), lhs=0)), collapse=TRUE)
  response <- terms(formula)[[2]]
  formula.km <- as.formula(paste0(deparse(response), "~ 1"))
  sfit <- survfit(formula.km, data.2)
  D <- data.frame(s = summary(sfit, fu)$surv, od = order(fu))
  D <- D[order(D$od),]
  logH <- log(-log(D$s))
  lm_fit <- lm(logH[status == 1] ~ -1 + b[status == 1, 1:n.knots])

  #Make proper naming of the initial values
  suppressWarnings(coefs <- lm_fit$coefficients)
  nr.active.knots <- ifelse(cure, length(knots) - 1, length(knots))
  ini_values <- c(coefs, rep(0, ncol(X) - 1))
  names(ini_values) <- c(paste0("gamma", 1:nr.active.knots), colnames(X)[-1])

  if(message) cat("Completed!\nFitting the model...")

  #Fit the model

  #fit <- stpm2(Surv(FU_years, status == 1) ~ -1, data = data.2, bhazard = bhazard,
  #             smooth.formula = ~basis(knots, log(FU_years)))

  nr.spline <- nr.active.knots
  pena.fun.spec <- function(beta) penalization.fun(beta, alpha = alpha, lambda = lambda)
  res <- optim(par = ini_values, fn = penalized_flexible_minuslog_likelihood,
               time = fu, status = status,
               b = b, db = db, bhazard = data.2[, bhazard],
               pena.fun = pena.fun.spec, nr.spline = nr.spline,
               control = list(maxit = 10000))

  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(message) cat("Completed!\n")

  #Compute the hessian matrix

  #Output the results
  L <- list(formula = formula,
            data = data.2, cure = cure,
            coefs = res$par,
            base_function = base_function, dbase_function = dbase_function,
            knots = knots, ML = res$value, df = length(res$par) - 1,
            formula = formula)

  class(L) <- "PenFlexRelsurvModel"
  L
}
