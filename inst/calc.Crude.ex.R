##Use data cleaned version of the colon disease data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 1000), ]

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

##Spline-base cure model
#Fit cure model
fit <- rstpm2::stpm2(Surv(FUyear, status) ~ 1, data = colonDC, df = 6,
                     bhazard = colonDC$bhaz, cure = TRUE)

#Compute the probability of disease-related death
res <- calc.Crude(fit, time = seq(0, 20, length.out = 50),
                  rmap = list(age = agedays, sex = sex, year = dx),
                  var.type = "n")
plot(res)

#Compute the conditional probability of dying from other causes than disease
res <- calc.Crude(fit, time = seq(0, 20, length.out = 50), type = "condother",
                  rmap = list(age = agedays, sex = sex, year = dx), var.type = "n")
plot(res)


#Simple parametric cure model
#Fit cure model
fit <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                      type = "mixture", dist = "weibull", link = "logit")

#Compute the probability of disease-related death
res <- calc.Crude(fit, time = seq(0, 20, length.out = 50),
                  rmap = list(age = agedays, sex = sex, year = dx),
                  var.type = "n")
plot(res)

#Compute the conditional probability of disease-related death
res2 <- calc.Crude(fit, time = seq(0, 20, length.out = 50), type = "condother",
                  rmap = list(age = agedays, sex = sex, year = dx), reverse = TRUE,
                  var.type = "n")
plot(res2)
