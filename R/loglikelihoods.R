#' Extract general population hazard
#'
#' Yearly general population hazards matched on age, gender, and calendar year is extracted from a ratetable
#'
#' @param time Either a numeric vector of follow-up times or a character indicating which variable is the follow-up time in the data
#' @param age Either a numeric vector of age or a character indicating which variable is the age in the data
#' @param sex A numeric vector including the time points for which to make predictions.
#' @param date A character denoting the type of prediction. Possible values are \code{relsurv} (default), \code{curerate}, and \code{probcure} (see details)
#' @param data The data from which to extract variables from. If \code{time}, \code{age}, \code{sex}, and \code{date} are not characters, leave this empty.
#' @param ratetable Object of class \code{ratetable} to extract background hazards from. Defaults to survexp.dk
#' @return An object of class \code{numeric} of yearly hazards.
#' @export
#' @example inst/general.haz.ex.R

general.haz <- function(time, age, sex, year, data = NULL, ratetable = survexp.dk){
  if(is.character(time)){
    time <- data[, time]
  }
  if(is.character(age)){
    age <- data[, age]
  }
  if(is.character(sex)){
    sex <- data[, sex]
  }
  if(is.character(year)){
    year <- data[, year]
  }

  dimid <- attr(ratetable, "dimid")
  od <- sapply(c("age", "sex", "year"), function(x) which(dimid == x))
  n <- length(time)

  haz <- rep(NA, n)
  sex_new <- as.character(sex)
  age_new <- pmin(round((age + time) / ayear), 99)
  year_eval <- format(year + time, "%Y")
  ryear <- range(as.numeric(dimnames(ratetable)[[od["year"]]]))
  year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
  year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)


  D <- data.frame(age = age_new, sex = sex_new, year = year_eval, stringsAsFactors = F)
  D <- D[, od]

  for(i in 1:n){
    haz[i] <- ratetable[D[i, 1], D[i, 2], D[i, 3]]
  }
  haz <- haz * ayear
  haz
}

#Global variable indicating the duration of a year
ayear <- 365.24

# Functions for extracting link functions
# Function for extracting the specified link function
get.link <- function(link){
  if(link == "logit"){
    function(x) exp(x) / (exp(x) + 1)
  }else if(link == "identity"){
    function(x) x
  }else if(link == "loglog"){
    function(x) exp(-exp(x))
  }else if(link == "probit"){
    function(x) pnorm(x)
  }else if(link == "nlogit"){
    function(x) exp(-x) / (exp(-x) + 1)
  }else if(link == "nprobit"){
    function(x) pnorm(-x)
  }else{
    stop("Link function should be either 'logit', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

get.dlink <- function(link){
  if(link == "logit"){
    function(x) exp(x) / ((exp(x) + 1) ^ 2)
  }else if(link == "identity"){
    function(x) 1
  }else if(link == "loglog"){
    function(x) -exp(x - exp(x))
  }else if(link == "probit"){
    function(x) -dnorm(x)
  }else if(link == "nlogit"){
    function(x) -exp(-x) / ((exp(-x) + 1) ^ 2)
  }else if(link == "nprobit"){
    function(x) -dnorm(-x)
  }else{
    stop("Link function should be either 'logit', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

get.inv.link <- function(link){
  if(link == "logit"){
    function(x) log(x / (1 - x))
  }else if(link == "identity"){
    function(x) x
  }else if(link == "loglog"){
    function(x) log(-log(x))
  }else if(link == "probit"){
    function(x) qnorm(x)
  }else if(link == "nlogit"){
    function(x) -log(x / (1 - x))
  }else if(link == "nprobit"){
    function(x) -qnorm(x)
  }else{
    stop("Link function should be either 'logit', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

# Function for extracting the specified survival function
get.surv <- function(dist){

  if(dist == "exponential"){

    return(function(x, lps) exp(-x * exp(lps[[2]])))

  }else if(dist == "weibull"){

    return(function(x, lps) exp(-x ^ exp(lps[[3]]) * exp(lps[[2]])))

  }else if(dist == "lognormal"){

    return(function(x, lps) 1 - pnorm((log(x) - lps[[2]]) / exp(lps[[3]])))

  }else{

    stop("Distribution should be either 'exponential', 'weibull', or 'lognormal'")

  }
}

# Function for extracting the specified density function
get.dens <- function(dist){
  if(dist == "exponential"){
    return(function(x, lps){
      scale <- exp(lps[[2]])
      scale * exp(-x * scale)
    })
  }else if(dist == "weibull"){
    return(function(x, lps){
      scale <- exp(lps[[2]])
      shape <- exp(lps[[3]])
      exp(-x ^ shape * scale) * shape * scale * x ^ (shape - 1)
    })
  }else if(dist == "lognormal"){
    return(function(x, lps){
      dnorm((log(x) - lps[[2]]) / exp(lps[[3]])) / (exp(lps[[3]]) * x)
    })
  }else{
    stop("Distribution should be either 'exponential', 'weibull', or 'lognormal'")
  }
}

#Function to calculate linear predictors for the simple parametric cure models
calc.lps <- function(Xs, param){
  lps <- vector("list", length(Xs))
  for(i in 1:length(Xs)){
    if(ncol(Xs[[i]]) != 0){
      lps[[i]] <- Xs[[i]] %*% param[1:ncol(Xs[[i]])]
      param <- param[-c(1:ncol(Xs[[i]]))]
    }
  }
  lps
}


######Likelihood functions
# Parametric Mixture cure model
mixture_minuslog_likelihood <- function(param, time, event, Xs, link_fun,
                                        surv_fun, dens_fun, bhazard){

  #Calculate linear predictors
  lps <- calc.lps(Xs, param)

  #Compute pi and the survival of the uncured
  pi <- link_fun(lps[[1]])
  surv <- surv_fun(time, lps)
  surv.term <- log(pi + (1 - pi) * surv)

  #Calculate hazard term only for uncensored patients.
  events <- which(event == 1)
  dens <- dens_fun(time[events], lapply(lps, function(lp) lp[events,]))
  pi.events <- pi[events]
  surv.events <- surv[events]
  haz.term <- log( bhazard[events] + dens * ( 1 - pi.events ) / ( pi.events + (1 - pi.events) * surv.events ))
  surv.term[events] <- surv.term[events] + haz.term

  #Output the negative log likelihood
  -sum(surv.term)
}

nmixture_minuslog_likelihood <- function(param, time, event, Xs, link_fun,
                                         surv_fun, dens_fun, bhazard){

  #Calculate linear predictors
  lps <- calc.lps(Xs, param)

  #Compute pi and the survival of the uncured
  pi <- link_fun(lps[[1]])
  surv <- surv_fun(time, lps)
  surv.term <- log(pi) - log(pi) * surv

  #Calculate hazard term only for uncensored patients.
  events <- which(event == 1)
  dens <- dens_fun(time[events], lapply(lps, function(lp) lp[events,]))
  pi.events <- pi[events]
  haz.term <- log( bhazard[events] - log(pi.events) * dens)
  surv.term[events] <- surv.term[events] + haz.term

  #Output the negative log likelihood
  -sum(surv.term)
}



# #Parametric mixture cure model 2
# mixture_minuslog_likelihood2 <- function(param, time, event, Xs, link_fun,
#                                         bhazard, surv_fun, dens_fun){
#   gamma <- param[grepl("gamma", names(param))]
#   coefs.k1 <- param[grepl("k1", names(param))]
#   coefs.k2 <- param[grepl("k2", names(param))]
#   coefs.k3 <- param[grepl("k3", names(param))]
#   lp.k1 <- Xs[[2]] %*% coefs.k1
#   if(length(coefs.k2) > 0){
#     lp.k2 <- Xs[[3]] %*% coefs.k2
#   }else{
#     lp.k2 <- NULL
#   }
#   if(length(coefs.k3) > 0){
#     lp.k3 <- Xs[[4]] %*% coefs.k3
#   }else{
#     lp.k3 <- NULL
#   }
#   deaths <- event == 1
#   pi <- link_fun(Xs[[1]] %*% gamma)
#   pi_deaths <- pi[deaths]
#   surv <- surv_fun(time, lp.k1, lp.k2, lp.k3)
#   dens <- dens_fun(time[deaths], lp.k1[deaths,], lp.k2[deaths,], lp.k3[deaths,])
#   terms <- log( pi + (1 - pi) * surv )
#   terms[deaths] <- terms[deaths] + log(bhazard[deaths] + dens * ( 1 - pi_deaths ) / ( pi_deaths + (1 - pi_deaths) * surv[deaths]))
#   -sum(terms)
# }
#
# # Parametric non-mixture cure model
# nmixture_minuslog_likelihood <- function(param, time, event, Xs, link_fun,
#                                          bhazard, surv_fun, dens_fun){
#   gamma <- param[grepl("gamma", names(param))]
#   coefs.k1 <- param[grepl("k1", names(param))]
#   coefs.k2 <- param[grepl("k2", names(param))]
#   coefs.k3 <- param[grepl("k3", names(param))]
#   lp.k1 <- Xs[[2]] %*% coefs.k1
#   if(length(coefs.k2) > 0){
#     lp.k2 <- Xs[[3]] %*% coefs.k2
#   }
#   if(length(coefs.k3) > 0){
#     lp.k3 <- Xs[[4]] %*% coefs.k3
#   }
#   pi <- link_fun(Xs[[1]] %*% gamma)
#   surv <- surv_fun(time, lp.k1, lp.k2, lp.k3)
#   dens <- dens_fun(time, lp.k1, lp.k2, lp.k3)
#   loglik <- event * log( bhazard - log(pi) * dens) + ( log(pi) - log(pi) * surv )
#   -sum(loglik)
# }

# Flexible mixture cure model
# flexible_mixture_minuslog_likelihood <- function(param, time, event, X, b, db, bhazard,
#                                                  link_fun_pi, link_fun_su, dlink_fun_su){
#   gamma <- param[1:ncol(X)]
#   beta <- param[(ncol(X) + 1):length(param)]
#   lp <- X %*% gamma
#   pi <- link_fun_pi(lp)
#   eta <- b %*% beta
#   deta <- db %*% beta / time
#   surv <- link_fun_su(eta)
#   rsurv <- pi + (1 - pi) * surv
#   dsurv <- dlink_fun_su(eta)
#   ehaz <- -(1 - pi) * dsurv * deta / rsurv
#   suppressWarnings(likterms <- event * log( bhazard + ehaz) + log(rsurv))
#   -sum(likterms)
# }

flexible_mixture_minuslog_likelihood <- function(param, time, event, X, b, db, bhazard,
                                                 link_fun_pi, link_fun_su, dlink_fun_su){
  #Get parameters
  gamma <- param[1:ncol(X)]
  beta <- param[(ncol(X) + 1):length(param)]

  #Calculate linear predictors
  lp <- X %*% gamma
  pi <- link_fun_pi(lp)
  eta <- b %*% beta
  surv <- link_fun_su(eta)
  rsurv <- pi + (1 - pi) * surv
  likterms <- log(rsurv)

  #Add the hazard term only for events
  events <- event == 1
  deta <- db[events,] %*% beta / time[events]
  dsurv <- dlink_fun_su(eta[events])
  ehaz <- -(1 - pi[events]) * dsurv * deta / rsurv[events]
  suppressWarnings(likterms[events] <- likterms[events] + log( bhazard[events] + ehaz))

  #Output the negative log likelihood
  -sum(likterms)
}


# flexible_mixture_minuslog_likelihood2 <- function(param, time, event, X2, b, db, bhazard,
#                                                  link_fun_pi, link_fun_su, dlink_fun_su){
#   #Get parameters
#   gamma <- param[1:ncol(X2)]
#   beta <- param[(ncol(X2) + 1):length(param)]
#
#   #Calculate linear predictors
#   lp <- X2 %*% gamma
#   pi <- link_fun_pi(lp)
#   eta <- b %*% beta
#   surv <- link_fun_su(eta)
#   rsurv <- pi + (1 - pi) * surv
#   likterms <- log(rsurv)
#
#   #Add the hazard term only for events
#   events <- event == 1
#   deta <- db[events,] %*% beta / time[events]
#   dsurv <- dlink_fun_su(eta[events])
#   ehaz <- -(1 - pi[events]) * dsurv * deta / rsurv[events]
#   suppressWarnings(likterms[events] <- likterms[events] + log( bhazard[events] + ehaz))
#
#   #Output the negative log likelihood
#   -sum(likterms)
# }


#Non-mixture cure
flexible_nmixture_minuslog_likelihood <- function(param, time, event, X, b, db, bhazard,
                                                  link_fun_pi, link_fun_su, dlink_fun_su){
  #Get parameters
  gamma <- param[1:ncol(X)]
  beta <- param[(ncol(X) + 1):length(param)]

  #Calculate linear predictors
  lp <- X %*% gamma
  pi <- link_fun_pi(lp)
  eta <- b %*% beta
  surv <- link_fun_su(eta)
  rsurv <- pi ^ (1 - surv)
  likterms <- log(rsurv)

  #Add the hazard term only for events
  events <- event == 1
  deta <- db[events,] %*% beta / time[events]
  ddist <- dlink_fun_su(eta[events])
  ehaz <- log(pi[events]) * ddist * deta
  suppressWarnings(likterms[events] <- likterms[events] + log( bhazard[events] + ehaz))

  #Output the negative log likelihood
  -sum(likterms)
}


# uncuredhazfun.mix <- function(pars, time, event, X, b, db, bhazard,
#                               link_fun_pi, link_fun_su, dlink_fun_su){
#   #Get parameters
#   beta <- pars[(ncol(X) + 1):length(pars)]
#
#   #Calculate linear predictors
#   eta <- b %*% beta
#   surv <- link_fun_su(eta)
#
#   #Add the hazard term only for events
#   deta <- db %*% beta / time
#   dsurv <- dlink_fun_su(eta)
#   ehaz <- - deta * dsurv / surv
#
#   c(ehaz)[1:10]
# }
#
# get_ui <- function(time, X, b, db, bhazard,
#                    link_fun_pi, link_fun_su, dlink_fun_su){
#
#   deta <- db %*% beta / time
#   dsurv <- dlink_fun_su(eta)
#   ehaz <- - deta * dsurv / surv
#
#   c(ehaz)[1:10]
# }


#Basis function
basis <- function(knots, x, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
  nx <- length(x)
  if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
  nk <- ncol(knots)
  b <- matrix(nrow=length(x), ncol=nk)
  if (nk>0){
    b[,1] <- 1
    b[,2] <- x
  }
  if (nk>2) {
    lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
    for (j in 1:(nk-2)) {
      b[,j+2] <- pmax(x - knots[,j+1], 0)^3 - lam[,j+1]*pmax(x - knots[,1], 0)^3 -
        (1 - lam[,j+1])*pmax(x - knots[,nk], 0)^3
    }
  }

  if(!intercept) b <- b[,-1]

  if(ortho){
    if(is.null(R.inv)){
      qr_decom <- qr(b)
      b <- qr.Q(qr_decom)
      R.inv <- solve(qr.R(qr_decom))
    } else{
      b <- b %*% R.inv
    }
  }
  attr(b, "R.inv") <- R.inv
  b
}

#Derivate of basis function
dbasis <- function(knots, x, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
  if(ortho & is.null(R.inv)) stop("Both 'ortho' and 'R.inv' has to be specified!")
  nx <- length(x)
  if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
  nk <- ncol(knots)
  b <- matrix(nrow=length(x), ncol=nk)
  if (nk>0){
    b[,1] <- 0
    b[,2] <- 1
  }
  if (nk>2) {
    lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
    for (j in 3:nk) {
      b[,j] <- 3*pmax(x - knots[,j-1], 0)^2 - 3*lam[,j-1]*pmax(x - knots[,1], 0)^2 -
        3*(1 - lam[,j-1])*pmax(x - knots[,nk], 0)^2
    }
  }

  if(!intercept) b <- b[, -1]
  if(ortho){
    b <- b %*% R.inv
  }

  b
}


lhs <- function(formula){
  if (length(formula) == 3) formula[[2]] else NULL
}


rhs <- function (formula)
  if (length(formula) == 3) formula[[3]] else formula[[2]]

"rhs<-" <- function (formula, value)
{
  newformula <- formula
  newformula[[length(formula)]] <- value
  newformula
}

# flexible_minuslog_likelihood <- function(param, time, event, b, db, bhazard){
#   eta <- b %*% param
#   deaths <- event == 1
#   deta <- db[deaths,] %*% param
#   exp_eta <- exp(eta)
#   inside_log <- bhazard[deaths] + deta * exp_eta[deaths] / time[deaths]
#   terms <- -exp_eta
#   suppressWarnings(terms[deaths] <- terms[deaths] + log(inside_log))
#   -sum(terms)
# }

# penalized_flexible_minuslog_likelihood <- function(param, time, event, b, db, bhazard, pena.fun, nr.spline){
#   eta <- b %*% param
#   deaths <- event == 1
#   deta <- db[deaths,] %*% param
#   exp_eta <- exp(eta)
#   inside_log <- bhazard[deaths] + deta * exp_eta[deaths] / time[deaths]
#   terms <- -exp_eta
#   suppressWarnings(terms[deaths] <- terms[deaths] + log(inside_log))
#   -sum(terms) + pena.fun(param[-(1:nr.spline)])
# }


# Cure base functions

basis.cure <- function(knots, x, ortho = TRUE, R.inv = NULL, intercept = TRUE){
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

  if(!intercept) b <- b[,-1]

  if(ortho){
    if(is.null(R.inv)){
      qr_decom <- qr(b)
      b <- qr.Q(qr_decom)
      R.inv <- solve(qr.R(qr_decom))
    } else{
      b <- b %*% R.inv
    }
  }
  attr(b, "R.inv") <- R.inv
  b

  b
}

dbasis.cure <- function(knots, x, ortho = TRUE, R.inv = NULL, intercept = TRUE){
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

  if(!intercept) b <- b[, -1]

  if(ortho){
    b <- b %*% R.inv
  }
  b
}



minuslog_likelihood <- function(param, time, event, b, db,
                                link_fun, dlink_fun){

  param_list <- split(param, rep(1:length(b), sapply(b, ncol)))

  Sks <- vector("list", length(b))
  inner_sum <- vector("list", length(b))
  for(i in 1:length(b)){
    lp <- b[[i]] %*% param_list[[i]]
    dlp <- db[[i]] %*% param_list[[i]]
    Sks[[i]] <- link_fun(lp)
    haz <- - dlink_fun(lp) / Sks[[i]] * dlp / time
    inner_sum[[i]] <- rep(0, length(time))
    suppressWarnings(inner_sum[[i]][event == i] <- log(haz[event == i] * Sks[[i]][event == i]))
  }

  sum_Fks <- length(b) - do.call("+", Sks)
  inner_sum <- do.call("+", inner_sum)

  event_logical <- as.numeric(event != 0)
  -sum(inner_sum + (1 - event_logical) * log(1 - sum_Fks))
}




