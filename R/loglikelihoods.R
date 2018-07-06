#' Extract general population hazard
#'
#' Yearly general population hazards matched on age, gender, and calendar year is extracted from a ratetable.
#'
#' @param time Either a numeric vector of follow-up times or a character indicating the variable
#' containing the follow-up times in the data.
#' @param age Either a numeric vector of ages or a character indicating the variable containing the patient ages in the data.
#' @param sex Either a character vector or factor with the sex of each patient
#' or a character indicating the variable containing the patient sex in the data.
#' @param year Either a vector of class \code{Date} with the calendar time points
#' or a character indicating the variable containing the calendar times in the data.
#' @param data The data from which to extract variables from.
#' If \code{time}, \code{age}, \code{sex}, or \code{year} are not characters, this will not be used.
#' @param ratetable Object of class \code{ratetable} to extract background hazards from. Defaults to \code{survexp.dk}.
#' @return An object of class \code{numeric} containing the yearly expected hazards.
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


general.haz2 <- function(time, rmap, data = NULL, ratetable = survexp.dk, scale = ayear){
  dimid <- attr(ratetable, "dimid")
  vars <-

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

  #D2 <- D
  #D2$sex <- ifelse(D2$sex == "male", 1, 2)
  #D2$year <- as.numeric(D2$year)
  #a <- match.ratetable(as.matrix(D2), ratetable)

  dim_names <- dimnames(ratetable)
  J <- data.frame(rates = c(ratetable), rep(as.numeric(dim_names[[1]]), length(dim_names[[2]]) * length(dim_names[[3]])),
                  rep(as.numeric(dim_names[[2]]), length(dim_names[[1]]) * length(dim_names[[3]])),
                  rep(dim_names[[3]], each = length(dim_names[[1]]) *  length(dim_names[[2]])))

  names(J)[-1] <- dimid

  #head(J)

  #ratetable["38", "1835",]

  merge(D, J)$rates * scale
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
get.surv <- function(dist, link.mix = function(x) log(x / (1 - x))){

  if(dist == "exponential"){

    return(function(x, lps) exp(-x * exp(lps[[2]])))

  }else if(dist == "weibull"){

    return(function(x, lps) exp(-x ^ exp(lps[[3]]) * exp(lps[[2]])))

  }else if(dist == "lognormal"){

    return(function(x, lps) 1 - pnorm((log(x) - lps[[2]]) / exp(lps[[3]])))

  }else if(dist == "weiwei"){

    return(function(x, lps){
      p <- link.mix(lps[[2]])
      scale1 <- exp(lps[[3]])
      shape1 <- exp(lps[[4]])
      scale2 <- exp(lps[[5]])
      shape2 <- exp(lps[[6]])
      wei1 <- exp(-x ^ shape1 * scale1)
      wei2 <- exp(-x ^ shape2 * scale2)
      p * wei1 + (1 - p) * wei2
    })

  }else if(dist == "weiexp"){

    return(function(x, lps){
      p <- link.mix(lps[[2]])
      scale1 <- exp(lps[[3]])
      shape1 <- exp(lps[[4]])
      scale2 <- exp(lps[[5]])
      wei1 <- exp(-x ^ shape1 * scale1)
      exp2 <- exp(-x * scale2)
      p * wei1 + (1 - p) * exp2
    })

  }else{

    stop("Distribution should be either 'exponential', 'weibull', 'lognormal', 'weiwei', or 'weiexp'")

  }
}

# Function for extracting the specified density function
get.dens <- function(dist, link.mix = function(x) log(x / (1 - x))){
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

  }else if(dist == "weiwei"){

    return(function(x, lps){
      scale1 <- exp(lps[[3]])
      shape1 <- exp(lps[[4]])
      scale2 <- exp(lps[[5]])
      shape2 <- exp(lps[[6]])
      p <- link.mix(lps[[2]])
      f_1 <- exp(-x ^ shape1 * scale1) * shape1 * scale1 * x ^ (shape1 - 1)
      f_2 <- exp(-x ^ shape2 * scale2) * shape2 * scale2 * x ^ (shape2 - 1)
      p * f_1 + (1 - p) * f_2
    })


  }else if(dist == "weiexp"){

    return(function(x, lps){
      scale1 <- exp(lps[[3]])
      shape1 <- exp(lps[[4]])
      scale2 <- exp(lps[[5]])
      p <- link.mix(lps[[2]])
      f_1 <- exp(-x ^ shape1 * scale1) * shape1 * scale1 * x ^ (shape1 - 1)
      f_2 <- exp(-x * scale2) * scale2
      p * f_1 + (1 - p) * f_2
    })


  }else {

    stop("Distribution should be either 'exponential', 'weibull', 'lognormal', 'weiwei', or 'weiexp'")

  }
}

`%call+%` <- function(left,right) call("+",left,right)

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
minuslog_likelihoodDelayed <- function(param, time, time0, event, Xs, ind0, link.fun,
                                       surv.fun, dens.fun, bhazard, cure.type){

  #Calculate linear predictors
  lps <- calc.lps(Xs, param)

  #Compute pi and the survival of the uncured
  pi <- link.fun(lps[[1]])
  surv <- surv.fun(time, lps)
  rsurv <- cure.type$surv(pi, surv)
  likterms <- log(rsurv)
  surv0 <- surv.fun(time0, lps)
  rsurv0 <- cure.type$surv(pi, surv0)
  likterms[ind0] <- likterms[ind0] - log(rsurv0[ind0])

  #Calculate hazard term only for uncensored patients.
  #Add the hazard term only for events
  dens <- dens.fun(time, lps)
  #pi.events <- pi[events]
  #surv.events <- surv[events]
  ehaz <- cure.type$haz(pi[event], -dens[event], rsurv[event])
  haz <- bhazard[event] + ehaz
  likterms[event] <- likterms[event] + log(haz)

  #Output the negative log likelihood
  -sum(likterms)
}

grad2 <- function(func,x,...) # would shadow numDeriv::grad()
{
  h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
  temp <- x+h
  h.hi <- temp-x
  temp <- x-h
  h.lo <- x-temp
  twoeps <- h.hi+h.lo
  nx <- length(x)
  ny <- length(func(x,...))
  if (ny==0L) stop("Length of function equals 0")
  df <- if(ny==1L) rep(NA, nx) else matrix(NA, nrow=nx,ncol=ny)
  for (i in 1L:nx) {
    hi2 <- lo2 <- hi <- lo <- x
    hi[i] <- x[i] + h.hi[i]
    hi2[i] <- x[i] + 2 * h.hi[i]
    lo[i] <- x[i] - h.lo[i]
    lo2[i] <- x[i] - 2 * h.lo[i]
    if (ny==1L)
      df[i] <- (func(lo2, ...) - 8 * func(lo, ...) + 8 * func(hi, ...) - func(hi2, ...))/ (12 * h)
    else df[i,] <- (func(lo2, ...) - 8 * func(lo, ...) + 8 * func(hi, ...) - func(hi2, ...))/ (12 * h)
  }
  return(df)
}



# Likelihood for mixture cure models
# mixture_minuslog_likelihoodDelayed <- function(param, time, time0, event, Xs, link.fun,
#                                                surv.fun, dens.fun, bhazard){
#
#   #Calculate linear predictors
#   lps <- calc.lps(Xs, param)
#
#   #Compute pi and the survival of the uncured
#   pi <- link.fun(lps[[1]])
#   surv <- surv.fun(time, lps)
#   surv.term <- log(pi + (1 - pi) * surv)
#   surv0 <- surv.fun(time0, lps)
#   surv.term0 <- log(pi + (1 - pi) * surv0)
#   surv.term <- surv.term - surv.term0
#
#   #Calculate hazard term only for uncensored patients.
#   events <- which(event == 1)
#   dens <- dens.fun(time[events], lapply(lps, function(lp) lp[events,]))
#   pi.events <- pi[events]
#   surv.events <- surv[events]
#   haz.term <- log( bhazard[events] + dens * ( 1 - pi.events ) / ( pi.events + (1 - pi.events) * surv.events ))
#   surv.term[events] <- surv.term[events] + haz.term
#
#   #Output the negative log likelihood
#   -sum(surv.term)
# }
#
# # Likelihood for non-mixture cure models
# nmixture_minuslog_likelihoodDelayed <- function(param, time, time0, event, Xs, link.fun,
#                                          surv.fun, dens.fun, bhazard){
#
#   #Calculate linear predictors
#   lps <- calc.lps(Xs, param)
#
#   #Compute pi and the survival of the uncured
#   pi <- link.fun(lps[[1]])
#   surv <- surv.fun(time, lps)
#   surv.term <- log(pi) - log(pi) * surv
#   surv0 <- surv.fun(time0, lps)
#   surv.term0 <- log(pi) - log(pi) * surv0
#   surv.term <- surv.term - surv.term0
#
#   #Calculate hazard term only for uncensored patients.
#   events <- which(event == 1)
#   dens <- dens.fun(time[events], lapply(lps, function(lp) lp[events,]))
#   pi.events <- pi[events]
#   haz.term <- log( bhazard[events] - log(pi.events) * dens)
#   surv.term[events] <- surv.term[events] + haz.term
#
#   #Output the negative log likelihood
#   -sum(surv.term)
# }



# flexible_mixture_minuslog_likelihood <- function(param, time, event, X, b, db, bhazard,
#                                                  link_fun_pi, link_fun_su, dlink_fun_su){
#   #Get parameters
#   gamma <- param[1:ncol(X)]
#   beta <- param[(ncol(X) + 1):length(param)]
#
#   #Calculate linear predictors
#   lp <- X %*% gamma
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


# #Non-mixture cure
# flexible_nmixture_minuslog_likelihood <- function(param, time, event, X, b, db, bhazard,
#                                                   link_fun_pi, link_fun_su, dlink_fun_su){
#   #Get parameters
#   gamma <- param[1:ncol(X)]
#   beta <- param[(ncol(X) + 1):length(param)]
#
#   #Calculate linear predictors
#   lp <- X %*% gamma
#   pi <- link_fun_pi(lp)
#   eta <- b %*% beta
#   surv <- link_fun_su(eta)
#   rsurv <- pi ^ (1 - surv)
#   likterms <- log(rsurv)
#
#   #Add the hazard term only for events
#   events <- event == 1
#   deta <- db[events,] %*% beta / time[events]
#   ddist <- dlink_fun_su(eta[events])
#   ehaz <- log(pi[events]) * ddist * deta
#   suppressWarnings(likterms[events] <- likterms[events] + log( bhazard[events] + ehaz))
#
#   #Output the negative log likelihood
#   -sum(likterms)
# }

mix <- list(surv = function(pi, surv) pi + (1 - pi) * surv,
            haz = function(pi, gradS, surv) - (1 - pi) * gradS / surv,
            dens = function(pi, gradS, surv) - (1 - pi) * gradS)

nmix <- list(surv = function(pi, surv) pi ^ (1 - surv),
             haz = function(pi, gradS, surv) log(pi) * gradS,
             dens = function(pi, gradS, surv) log(pi) * gradS * surv)

#Likelihood for generalized mixture cure models (including delayed entry - code from rstpm2)
GenFlexMinLogLikDelayed <- function(param, event, X, XD, X.cr, X0, ind0, bhazard,
                                    link.type.cr, link.surv, kappa, constraint,
                                    cure.type){
  #Get parameters
  gamma <- param[1:ncol(X.cr)]
  beta <- param[(ncol(X.cr) + 1):length(param)]

  #Calculate linear predictors
  eta.pi <- X.cr %*% gamma
  pi <- get.link(link.type.cr)(eta.pi)
  eta <- X %*% beta
  surv <- link.surv$ilink(eta)
  rsurv <- cure.type$surv(pi, surv)
  likterms <- log(rsurv)
  eta0 <- X0 %*% beta
  surv0 <- link.surv$ilink(eta0)
  rsurv0 <- cure.type$surv(pi[ind0], surv0)
  likterms[ind0] <- likterms[ind0] - log(rsurv0)

  #Add the hazard term only for events
  etaD <- XD %*% beta
  ehaz <- cure.type$haz(pi[event], link.surv$gradS(eta[event], etaD[event]), rsurv[event])
  haz <- haz.const <- bhazard[event] + ehaz
  haz[haz <= 0] <- .Machine$double.eps
  likterms[event] <- likterms[event] + log(haz)

  #Calculate hazard for which constraints are used
  if(constraint) haz.const <- link.surv$h(eta, etaD)

  #Output the negative log likelihood
  -sum( likterms ) + kappa / 2 * sum( ( haz.const[haz.const < 0] ) ^ 2 )
}


#Likelihood for generalized mixture cure models (including delayed entry - code from rstpm2)
# GenFlexMixMinLogLikDelayed <- function(param, event, X, XD, X.cr, X0, ind0, bhazard,
#                                 link.type.cr, link.surv, kappa, constraint){
#   #Get parameters
#   gamma <- param[1:ncol(X.cr)]
#   beta <- param[(ncol(X.cr) + 1):length(param)]
#
#   #Calculate linear predictors
#   eta.pi <- X.cr %*% gamma
#   pi <- get.link(link.type.cr)(eta.pi)
#   eta <- X %*% beta
#   surv <- link.surv$ilink(eta)
#   rsurv <- pi + (1 - pi) * surv
#   likterms <- log(rsurv)
#   eta0 <- X0 %*% beta
#   surv0 <- link.surv$ilink(eta0)
#   rsurv0 <- pi[ind0] + (1 - pi[ind0]) * surv0
#   likterms[ind0] <- likterms[ind0] - log(rsurv0)
#
#   #Add the hazard term only for events
#   etaD <- XD %*% beta
#   ehaz <- - ( 1 - pi[event] ) * link.surv$gradS(eta[event], etaD[event]) / rsurv[event]
#   haz <- haz.const <- bhazard[event] + ehaz
#   haz[haz <= 0] <- .Machine$double.eps
#   likterms[event] <- likterms[event] + log(haz)
#
#   #Calculate hazard for which constraints are used
#   if(constraint) haz.const <- link.surv$h(eta, etaD)
#
#   #Output the negative log likelihood
#   -sum( likterms ) + kappa / 2 * sum( ( haz.const[haz.const < 0] ) ^ 2 )
# }
#
#
# #Likelihood for generalized non-mixture cure models (including delayed entry - code from rstpm2)
# GenFlexNmixMinLogLikDelayed <- function(param, event, X, XD, X.cr, X0, ind0,
#                                        bhazard, link.type.cr, link.surv, kappa, constraint){
#   #Get parameters
#   gamma <- param[1:ncol(X.cr)]
#   beta <- param[(ncol(X.cr) + 1):length(param)]
#
#   #Calculate linear predictors
#   eta.pi <- X.cr %*% gamma
#   pi <- get.link(link.type.cr)(eta.pi)
#   eta <- X %*% beta
#   surv <- link.surv$ilink(eta)
#   rsurv <- pi ^ (1 - surv)
#   likterms <- log(rsurv)
#   eta0 <- X0 %*% beta
#   surv0 <- link.surv$ilink(eta0)
#   rsurv0 <- pi[ind0] ^ (1 - surv0)
#   likterms[ind0] <- likterms[ind0] - log(rsurv0)
#
#   #Add the hazard term only for events
#   etaD <- XD %*% beta
#   ehaz <- log( pi[event] ) * link.surv$gradS( eta[event], etaD[event] )
#   haz <- haz.const <- bhazard[event] + ehaz
#   haz[ haz < 0 ] <- .Machine$double.eps
#   likterms[event] <- likterms[event] + log( haz )
#
#   #Calculate hazard for which constraints are used
#   if(constraint) haz.const <- link.surv$h(eta, etaD)
#
#   #Output the negative log likelihood
#   -sum(likterms) + kappa / 2 * sum( ( haz.const[haz.const < 0] ) ^ 2 )
# }


#Basis function
#' @export
cb <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
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

  if(!intercept) b <- b[,-1, drop = FALSE]

  if(ortho){
    if(is.null(R.inv)){
      qr_decom <- qr(b)
      b <- qr.Q(qr_decom)
      R.inv <- solve(qr.R(qr_decom))
    } else{
      R.inv <- matrix(R.inv, nrow = ncol(b))
      b <- b %*% R.inv
    }
  } else {
    R.inv <- diag(ncol(b))
  }

  a <- list(knots = knots[1,], ortho = ortho, R.inv = R.inv, intercept = intercept)
  attributes(b) <- c(attributes(b), a)
  class(b) <- c("cb", "matrix")
  b
}

#Predict function associated with bsx.
predict.cb <- function (object, newx, ...)
{
  if (missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "ortho",
                                              "R.inv", "intercept")])
  do.call("cb", a)
}

#Additional function needed to fix the knot location in cases where df is only specified
#' @export
makepredictcall.cb <- function (var, call)
{
  if (as.character(call)[1L] != "cb")
    return(call)
  at <- attributes(var)[c("knots", "ortho", "R.inv", "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

#Derivate of basis function
dbasis <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
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

###############KIG HER LASSSSSEEEEEE##################3
#Fix problem med cb.cure! Det virker med orthogonalisering.

# Cure base functions
#' @export
cbc <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE){
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

  if(!intercept) b <- b[,-1, drop = F]

  if(ortho){
    if(is.null(R.inv)){
      qr_decom <- qr(b)
      b <- qr.Q(qr_decom)
      R.inv <- solve(qr.R(qr_decom))
    } else{
      R.inv <- matrix(R.inv, nrow = ncol(b))
      b <- b %*% R.inv
    }
  } else {
    R.inv <- diag(ncol(b))
  }

  a <- list(knots = knots, ortho = ortho, R.inv = R.inv, intercept = intercept)
  attributes(b) <- c(attributes(b), a)
  class(b) <- c("cbc", "matrix")
  b
}

#Predict function associated with bsx.
predict.cbc <- function (object, newx, ...)
{
  if (missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "ortho",
                                              "R.inv", "intercept")])
  do.call("cbc", a)
}

#Additional function needed to fix the knot location in cases where df is only specified
#' @export
makepredictcall.cbc <- function (var, call)
{
  if (as.character(call)[1L] != "cbc")
    return(call)
  at <- attributes(var)[c("knots", "ortho", "R.inv", "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
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



#legendre.quadrature.rule.200 <- legendre.quadrature.rules(200)[[200]]
#save(legendre.quadrature.rule.200, file = "data/legendre.quadrature.rules.200.RData")


