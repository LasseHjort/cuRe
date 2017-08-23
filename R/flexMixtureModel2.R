FlexMixtureCureModel2 <- function(formula, data, bhazard, smooth.formula = ~ 1,
                                 knots = NULL, n.knots = NULL,
                                 tvc.formula = NULL, knots.time = NULL, n.knots.time = NULL,
                                 hes = T, constr = F, message = T,
                                 type = "mixture", link = "logistic", start = 0.05){

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
  }


  #Construct design matrix
  X <- model.matrix(smooth.formula, data = data)
  b <- cbind(b, X[, -1], b_time)
  db <- cbind(db, matrix(0, ncol = ncol(X) - 1, nrow = nrow(X)), db_time)
  X <- model.matrix(formula, data = data)

  #Generate initial values by running a mixture weibull cure model
  if(message) cat("Finding initial values...")
  vars <- c(all.vars(smooth.formula), all.vars(tvc.formula))
  formula.2 <- reformulate(termlabels = ifelse(length(vars) == 0, "1", vars),
                           response = formula[[2]])

  status2 <- 1 - status
  fit_glm <- glm(status2 ~ -1 + X, family = binomial(link = "logit"))
  #ini_pi <- c(log(start / (1 - start)), rep(0, ncol(X) - 1))
  ini_pi <- fit_glm$coefficients

  fit <- coxph(formula.2, data = data[status == 1,])
  cum_base_haz <- get_basehaz(fit)
  suppressWarnings(logH <- log(cum_base_haz$hazard))
  fit_lm <- lm(logH ~ -1 + b[status == 1,])

  #Make proper naming of the initial values
  suppressWarnings(coefs <- fit_lm$coefficients)
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
  if(message) cat("Completed!\nFitting the model...")

  if(type == "mixture"){
    minusloglik <- flexible_mixture_minuslog_likelihood
  }else{
    minusloglik <- flexible_nmixture_minuslog_likelihood
  }

  #Extract link function
  link_fun <- get_link(link)
  #Fit the model
  if(constr){
    res <- constrOptim(theta = ini_values, f = minusloglik, grad = NULL,
                       ui = cbind(matrix(0, ncol = ncol(X), nrow = nrow(db)), db),
                       ci = rep(0, nrow(db)),
                       time = fu, status = status, X = X,
                       b = b, db = db, bhazard = bhazard,
                       link_fun = link_fun,
                       control = list(maxit = 10000))
  }else{
    res <- optim(ini_values, fn = minusloglik,
                 time = fu, status = status, X = X,
                 b = b, db = db, bhazard = bhazard,
                 link_fun = link_fun,
                 control = list(maxit = 10000))
  }

  names(res$par) <- gsub("spline_", "", names(res$par))
  if(res$convergence != 0){
    warning("Convergence not reached")
  }
  if(message) cat("Completed!\n")

  #Compute the hessian matrix
  if(hes){
    cov <- solve(pracma::hessian(minusloglik, res$par,
                                 time = fu, status = status, X = X, b = b,
                                 db = db, bhazard = bhazard, link_fun = link_fun))
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
            formula_main = smooth.formula, type = type, link = link_fun)

  class(L) <- "FlexCureModel"
  L
}
