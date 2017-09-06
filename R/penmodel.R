

minuslikelihood <- function(pars, fu, status, b, db){
  eta <- b %*% pars
  deta <- db %*% pars
  surv <- exp(-exp(eta))
  haz <- exp(eta) * deta / fu
  S <-
  q <- length(knots) + 2 ; S <- matrix(0, q, q)
  S[3:q,3:q] <- outer(knots, knots, FUN = basis, knots = knots)
  -sum(status * log(haz) + log(surv))
}


