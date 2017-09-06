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
#' @param opposite Logical (default \code{FALSE}) to indicate the order in the the ratetable
#' @return An object of class \code{numeric} of yearly hazards
#' @export

# Function for extracting expected hazard
extract_general <- function(time, age, sex, date, data = NULL, ratetable = survexp.dk, opposite = F){
  haz <- rep(NA, nrow(data))
  if(is.null(data)){
    sex_new <- as.character(sex)
    age_new <- pmin(round((age + time) / year), 99)
    year_eval <- format(date + time, "%Y")
    year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
    year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)
    for(i in 1:length(time)){
      haz[i] <- survexp.dk[age_new[i], sex_new[i], year_eval[i]]
    }
  }else{
    sex_new <- as.character(data[, sex])
    age_new <- pmin(round((data[, age] + data[, time]) / year), 99)
    year_eval <- format(data[, date] + data[, time], "%Y")
    if(!opposite){
      ryear <- range(as.numeric(names(ratetable[1,1,])))
      year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
      year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)
      for(i in 1:nrow(data)){
        haz[i] <- ratetable[age_new[i], sex_new[i], year_eval[i]]
      }
    }else{
      ryear <- range(as.numeric(names(ratetable[1,,1])))
      year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
      year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)
      for(i in 1:nrow(data)){
        haz[i] <- ratetable[age_new[i], year_eval[i], sex_new[i]]
      }
    }
  }
  haz <- haz * year
  haz
}

# Functions for extracting link functions
# Function for extracting the specified link function
get_link <- function(link){
  if(link == "logistic"){
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
    stop("Link function should be either 'logistic', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

get_dlink <- function(link){
  if(link == "logistic"){
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
    stop("Link function should be either 'logistic', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

get_inv_link <- function(link){
  if(link == "logistic"){
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
    stop("Link function should be either 'logistic', 'nlogit', 'nprobit', 'identity', 'probit', or 'loglog'")
  }
}

# Function for extracting the specified survival function
get_surv <- function(dist){
  if(dist == "exponential"){
    return(function(x, lp.k1, lp.k2, lp.k3) exp(-x * exp(lp.k1)))
  }else if(dist == "weibull"){
    return(function(x, lp.k1, lp.k2, lp.k3) exp(-x ^ exp(lp.k2) * exp(lp.k1)))
  }else if(dist == "lognormal"){
    return(function(x, lp.k1, lp.k2, lp.k3) 1 - pnorm((log(x) - lp.k1) / exp(lp.k2)))
  }else{
    stop("Distribution should be either 'exponential', 'weibull', or 'lognormal'")
  }
}

# Function for extracting the specified density function
get_dens <- function(dist){
  if(dist == "exponential"){
    return(function(x, lp.k1, lp.k2, lp.k3){
      scale <- exp(lp.k1)
      scale * exp(-x * scale)
    })
  }else if(dist == "weibull"){
    return(function(x, lp.k1, lp.k2, lp.k3){
      scale <- exp(lp.k1)
      shape <- exp(lp.k2)
      exp(-x ^ shape * scale) * shape * scale * x ^ (shape - 1)
    })
  }else if(dist == "lognormal"){
    return(function(x, lp.k1, lp.k2, lp.k3){
      dnorm((log(x) - lp.k1) / exp(lp.k2)) / (exp(lp.k2) * x)
    })
  }else{
    stop("Distribution should be either 'exponential', 'weibull', or 'lognormal'")
  }
}

get_surv2 <- function(dist){
  if(dist == "exponential"){
    return(function(x, k1, k2, k3, lp1, lp2, lp3) exp(-x * exp(lp1 %*% k1)))
  }else if(dist == "weibull"){
    return(function(x, k1, k2, k3, lp1, lp2, lp3) exp(-x ^ exp(lp2 %*% k2) * exp(lp1 %*% k1)))
  }else if(dist == "lognormal"){
    return(function(x, k1, k2, k3, lp1, lp2, lp3) 1 - pnorm((log(x) - lp1 %*% k1) / exp(lp2 %*% k2)))
  }else{
    stop("Distribution should be either 'exponential', 'weibull', or 'lognormal'")
  }
}


######Likelihood functions
# Parametric Mixture cure model
mixture_minuslog_likelihood <- function(param, time, status, Xs, link_fun,
                                         bhazard, surv_fun, dens_fun){
  gamma <- param[grepl("gamma", names(param))]
  coefs.k1 <- param[grepl("k1", names(param))]
  coefs.k2 <- param[grepl("k2", names(param))]
  coefs.k3 <- param[grepl("k3", names(param))]
  lp.k1 <- Xs[[2]] %*% coefs.k1
  if(length(coefs.k2) > 0){
    lp.k2 <- Xs[[3]] %*% coefs.k2
  }else{
    lp.k2 <- NULL
  }
  if(length(coefs.k3) > 0){
    lp.k3 <- Xs[[4]] %*% coefs.k3
  }else{
    lp.k3 <- NULL
  }
  pi <- link_fun(Xs[[1]] %*% gamma)
  surv <- surv_fun(time, lp.k1, lp.k2, lp.k3)
  dens <- dens_fun(time, lp.k1, lp.k2, lp.k3)
  first_term <- status * log( bhazard + dens * ( 1 - pi ) / ( pi + (1 - pi) * surv ))
  second_term <- log( pi + (1 - pi) * surv )
  -sum(first_term + second_term)
}

#Parametric mixture cure model 2
mixture_minuslog_likelihood2 <- function(param, time, status, Xs, link_fun,
                                        bhazard, surv_fun, dens_fun){
  gamma <- param[grepl("gamma", names(param))]
  coefs.k1 <- param[grepl("k1", names(param))]
  coefs.k2 <- param[grepl("k2", names(param))]
  coefs.k3 <- param[grepl("k3", names(param))]
  lp.k1 <- Xs[[2]] %*% coefs.k1
  if(length(coefs.k2) > 0){
    lp.k2 <- Xs[[3]] %*% coefs.k2
  }else{
    lp.k2 <- NULL
  }
  if(length(coefs.k3) > 0){
    lp.k3 <- Xs[[4]] %*% coefs.k3
  }else{
    lp.k3 <- NULL
  }
  deaths <- status == 1
  pi <- link_fun(Xs[[1]] %*% gamma)
  pi_deaths <- pi[deaths]
  surv <- surv_fun(time, lp.k1, lp.k2, lp.k3)
  dens <- dens_fun(time[deaths], lp.k1[deaths,], lp.k2[deaths,], lp.k3[deaths,])
  terms <- log( pi + (1 - pi) * surv )
  terms[deaths] <- terms[deaths] + log(bhazard[deaths] + dens * ( 1 - pi_deaths ) / ( pi_deaths + (1 - pi_deaths) * surv[deaths]))
  -sum(terms)
}

# Parametric non-mixture cure model
nmixture_minuslog_likelihood <- function(param, time, status, Xs, link_fun,
                                         bhazard, surv_fun, dens_fun){
  gamma <- param[grepl("gamma", names(param))]
  coefs.k1 <- param[grepl("k1", names(param))]
  coefs.k2 <- param[grepl("k2", names(param))]
  coefs.k3 <- param[grepl("k3", names(param))]
  lp.k1 <- Xs[[2]] %*% coefs.k1
  if(length(coefs.k2) > 0){
    lp.k2 <- Xs[[3]] %*% coefs.k2
  }
  if(length(coefs.k3) > 0){
    lp.k3 <- Xs[[4]] %*% coefs.k3
  }
  pi <- link_fun(Xs[[1]] %*% gamma)
  surv <- surv_fun(time, lp.k1, lp.k2, lp.k3)
  dens <- dens_fun(time, lp.k1, lp.k2, lp.k3)
  loglik <- status * log( bhazard - log(pi) * dens) + ( log(pi) - log(pi) * surv )
  -sum(loglik)
}

# Flexible mixture cure model
flexible_mixture_minuslog_likelihood <- function(param, time, status, X, b, db, bhazard,
                                                 link_fun_pi, link_fun_su, dlink_fun_su){
  gamma <- param[1:ncol(X)]
  beta <- param[(ncol(X) + 1):length(param)]
  lp <- X %*% gamma
  pi <- link_fun_pi(lp)
  eta <- b %*% beta
  deta <- db %*% beta / time
  surv <- link_fun_su(eta)
  rsurv <- pi + (1 - pi) * surv
  dsurv <- dlink_fun_su(eta)
  ehaz <- -(1 - pi) * dsurv * deta / rsurv
  suppressWarnings(likterms <- status * log( bhazard + ehaz) + log(rsurv))
  -sum(likterms)
}

# Flexible non-mixture cure model
flexible_nmixture_minuslog_likelihood <- function(param, time, status, X, b, db, bhazard,
                                                  link_fun_pi, link_fun_su, dlink_fun_su){
  gamma <- param[1:ncol(X)]
  beta <- param[(ncol(X) + 1):length(param)]
  lp <- X %*% gamma
  pi <- link_fun_pi(lp)
  eta <- b %*% beta
  deta <- db %*% beta / time
  surv <- link_fun_su(eta)
  rsurv <- pi ^ (1 - surv)
  dsurv <- dlink_fun_su(eta)
  ehaz <- -log(pi) * dsurv * deta
  suppressWarnings(likterms <- status * log( bhazard + ehaz) + log(rsurv))
  -sum(likterms)
}

# flexible_minuslog_likelihood <- function(param, time, status, b, db, bhazard){
#   eta <- b %*% param
#   deaths <- status == 1
#   deta <- db[deaths,] %*% param
#   exp_eta <- exp(eta)
#   inside_log <- bhazard[deaths] + deta * exp_eta[deaths] / time[deaths]
#   terms <- -exp_eta
#   suppressWarnings(terms[deaths] <- terms[deaths] + log(inside_log))
#   -sum(terms)
# }

# penalized_flexible_minuslog_likelihood <- function(param, time, status, b, db, bhazard, pena.fun, nr.spline){
#   eta <- b %*% param
#   deaths <- status == 1
#   deta <- db[deaths,] %*% param
#   exp_eta <- exp(eta)
#   inside_log <- bhazard[deaths] + deta * exp_eta[deaths] / time[deaths]
#   terms <- -exp_eta
#   suppressWarnings(terms[deaths] <- terms[deaths] + log(inside_log))
#   -sum(terms) + pena.fun(param[-(1:nr.spline)])
# }


# Cure base functions

# basis_cure <- function(knots, x){
#   nk <- length(knots)
#   b <- matrix(nrow = length(x), ncol = nk - 1)
#   knots_rev <- rev(knots)
#   if (nk > 0) {
#     b[, 1] <- 1
#   }
#   if (nk > 2) {
#     for (j in 2:(nk - 1)) {
#       lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
#       b[, j] <- pmax(knots_rev[nk - j + 1] - x, 0)^3 - lam * pmax(knots_rev[nk] - x, 0)^3 -
#         (1 - lam) * pmax(knots_rev[1] - x, 0)^3
#     }
#   }
#   b
# }
#
# dbasis_cure <- function(knots, x){
#   nk <- length(knots)
#   b <- matrix(nrow = length(x), ncol = nk - 1)
#   knots_rev <- rev(knots)
#   if (nk > 0) {
#     b[, 1] <- 0
#   }
#   if (nk > 2) {
#     for (j in 2:(nk - 1)) {
#       lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
#       b[, j] <- - 3 * pmax(knots_rev[nk - j + 1] - x, 0)^2 + 3 * lam * pmax(knots_rev[nk] - x, 0)^2 +
#         3 * (1 - lam) * pmax(knots_rev[1] - x, 0)^2
#     }
#   }
#   b
# }

