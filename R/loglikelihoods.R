#' Extract general population hazard
#'
#' Yearly general population hazards matched on age, gender, and calendar year is extracted from a ratetable.
#'
#' @param time Either a numeric vector of follow-up times (in days) or a character indicating the variable
#' containing the follow-up times in the data.
#' @param age Either a numeric vector of ages (in days) or a character indicating the variable containing the patient ages in the data.
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

general.haz <- function(time, age, sex, year, data = NULL, ratetable = cuRe::survexp.dk){
  #Get variables if only characters are provided
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

  #Compute order of the dimnames in the ratetable - may not match with the attributes
  dimid <- attr(ratetable, "dimid")
  od <- sapply(c("age", "sex", "year"), function(x) which(dimid == x))
  n <- length(time)

  #Compute max age in the ratetable
  max_age <- as.numeric(attr(ratetable, "dimnames")[[od["age"]]])

  #Format sex
  sex_new <- as.character(sex)

  #Compute age and year after "time" days
  age_new <- pmin(round((age + time) / ayear), max_age)
  year_eval <- format(year + time, "%Y")

  #If years expand that available in the ratetable take the closest available data points
  ryear <- range(as.numeric(dimnames(ratetable)[[od["year"]]]))
  year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
  year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)


  #Create useful data frame from "ratetable" input and use dimnames from the ratetable
  df <- melt(as.matrix(ratetable), as.is = T)
  names(df)[1:length(dimid)] <- dimid

  #Create dataframe with all individuals and their characteristics
  D <- data.frame(age = age_new, sex = sex_new, year = year_eval, stringsAsFactors = F)
  #Order the dataframe to match the dimensions from the ratetable
  D <- D[, od]
  #Define order of individuals to go back after merging
  D$ord <- 1:nrow(D)

  #Merge the individual data frame with melted ratetable - only keep rows in individual data frame
  haz_df <- merge(D, df, by = c("age", "year", "sex"), all.x = T)

  #Extract ratetable values and use same order as input
  haz <- haz_df$value[order(haz_df$ord)]

  #Multiple by 365.24 and output
  haz * ayear
}

# attr <- attributes(ratetable)
# attr$dim
#
# a <- survexp.dk[,,"male"]
# b <- survexp.dk[,,"female"]
# a <- melt(as.matrix(a))
#
#
# b <- melt(as.matrix(b))
# a$sex <- "male"
# b$sex <- "female"
#
# attr$dimid
# df <- rbind(a, b)
#
# names(df)[1:2] <- attr$dimid[1:2]
#
#
# d <- merge(df, D, by = c("age", "sex", "year"), all.y = T)
# d$value


# general.haz2 <- function(time, rmap, data = NULL, ratetable = survexp.dk, scale = ayear){
#   dimid <- attr(ratetable, "dimid")
#   vars <-
#
#     od <- sapply(c("age", "sex", "year"), function(x) which(dimid == x))
#   n <- length(time)
#
#   haz <- rep(NA, n)
#   sex_new <- as.character(sex)
#   age_new <- pmin(round((age + time) / ayear), 99)
#   year_eval <- format(year + time, "%Y")
#   ryear <- range(as.numeric(dimnames(ratetable)[[od["year"]]]))
#   year_eval <- ifelse(year_eval < ryear[1], ryear[1], year_eval)
#   year_eval <- ifelse(year_eval > ryear[2], ryear[2], year_eval)
#
#
#   D <- data.frame(age = age_new, sex = sex_new, year = year_eval, stringsAsFactors = F)
#   D <- D[, od]
#
#   #D2 <- D
#   #D2$sex <- ifelse(D2$sex == "male", 1, 2)
#   #D2$year <- as.numeric(D2$year)
#   #a <- match.ratetable(as.matrix(D2), ratetable)
#
#   dim_names <- dimnames(ratetable)
#   J <- data.frame(rates = c(ratetable), rep(as.numeric(dim_names[[1]]), length(dim_names[[2]]) * length(dim_names[[3]])),
#                   rep(as.numeric(dim_names[[2]]), length(dim_names[[1]]) * length(dim_names[[3]])),
#                   rep(dim_names[[3]], each = length(dim_names[[1]]) *  length(dim_names[[2]])))
#
#   names(J)[-1] <- dimid
#
#   #head(J)
#
#   #ratetable["38", "1835",]
#
#   merge(D, J)$rates * scale
# }


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

  }else if(dist == "gmw"){

    return(function(x, lps){
      scale  <- exp(lps[[2]])
      shape1 <- exp(lps[[3]])
      frag   <- exp(lps[[4]])
      shape2 <- exp(lps[[5]])
      1 - (1 - exp(-x ^ shape1 * scale * exp(x * frag))) ^ shape2
    })

  }else{

    stop("Distribution should be either 'exponential', 'weibull', 'lognormal', 'weiwei', 'weiexp', or 'gmw'")

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

  }else if(dist == "gmw"){

    return(function(x, lps){
      scale  <- exp(lps[[2]])
      shape1 <- exp(lps[[3]])
      frag   <- exp(lps[[4]])
      shape2 <- exp(lps[[5]])
      num    <- scale * shape2 * x ^ (shape1 - 1) * (shape1 + frag * x) *
                exp(frag * x - scale * x ^ shape1 * exp(frag * x))
      denom  <- (1 - exp(-scale * x ^ shape1 * exp(frag * x))) ^ (1 - shape2)
      num / denom
    })

  }else {

    stop("Distribution should be either 'exponential', 'weibull', 'lognormal', 'weiwei', 'weiexp', or 'gmw'")

  }
}

`%call+%` <- function(left,right) call("+",left,right)

#Function to calculate linear predictors for the simple parametric cure models
calc.lps <- function(Xs, param){
  lps <- vector("list", length(Xs))
  for(i in 1:length(Xs)){
    if(ncol(Xs[[i]]) != 0){
      lps[[i]] <- Xs[[i]] %*% param[1:ncol(Xs[[i]])]
      param    <- param[-c(1:ncol(Xs[[i]]))]
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
  pi             <- link.fun(lps[[1]])
  surv           <- surv.fun(time, lps)
  rsurv          <- cure.type$surv(pi, surv)
  likterms       <- log(rsurv)

  surv0          <- surv.fun(time0, lps)
  rsurv0         <- cure.type$surv(pi, surv0)
  likterms[ind0] <- likterms[ind0] - log(rsurv0[ind0])

  #Calculate hazard term only for uncensored patients.
  #Add the hazard term only for events
  dens <- dens.fun(time, lps)
  #pi.events <- pi[events]
  #surv.events <- surv[events]
  ehaz            <- cure.type$haz(pi[event], -dens[event], rsurv[event])
  haz             <- bhazard[event] + ehaz
  likterms[event] <- likterms[event] + log(haz)

  #Output the negative log likelihood
  -sum(likterms)
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





uniroot.all <- function (f, interval, lower = min(interval), upper = max(interval),
                         tol = .Machine$double.eps^0.2, maxiter = 1000, n = 100, ...)
{
  if (!missing(interval) && length(interval) != 2)
    stop("'interval' must be a vector of length 2")
  if (!is.numeric(lower) || !is.numeric(upper) || lower >=
      upper)
    stop("lower < upper  is not fulfilled")
  xseq <- seq(lower, upper, len = n + 1)
  mod <- f(xseq, ...)
  Equi <- xseq[which(mod == 0)]
  ss <- mod[1:n] * mod[2:(n + 1)]
  ii <- which(ss < 0)
  for (i in ii) Equi <- c(Equi, uniroot(f, lower = xseq[i],
                                        upper = xseq[i + 1], ...)$root)
  return(Equi)
}
