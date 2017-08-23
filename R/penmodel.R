
event.times <- DLBCL$FU_years[DLBCL$status == 1]
knots <- log(quantile(event.times, probs = seq(0, 1, length.out = 10)))
b <- basis(x = log(DLBCL$FU_years), knots = knots)
db <- dbasis(x = log(DLBCL$FU_years), knots = knots)


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


sfit <- survfit(Surv(FU_years, status) ~ 1, data = DLBCL)
plot(sfit)

pred_surv <- sapply(DLBCL$FU_years, function(t) summary(sfit, t)$surv)
logH <- log(-log(pred_surv))
fit.lm <- lm(logH[DLBCL$status == 1] ~ -1 + b[DLBCL$status == 1,])
f <- function(t) exp(-exp(basis(x = log(t), knots = knots) %*% fit.lm$coefficients))
curve(f, from = 0, to = 16, add = T)

ini_val <- fit.lm$coefficients
optim(par = ini_val, fn = minuslikelihood, fu = DLBCL$FU_years,
      status = DLBCL$status, b = b, db = db)
