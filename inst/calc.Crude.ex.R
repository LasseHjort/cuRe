##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

#Fit flexible parametric cure model
fit <- stpm2(Surv(FUyear, status2) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz, cure = T)

#Compute the probability of cancer related death
res <- calc.Crude(fit, time = seq(0, 20, length.out = 100),
                  rmap = list(age = agedays, sex = sex, year = dx))
plot(res)

#Compute the probability of eventually dying from other causes than cancer
res <- calc.Crude(fit, time = seq(0, 20, length.out = 100), type = "othertime",
                  rmap = list(age = agedays, sex = sex, year = dx))
plot(res)


#Fit parametric cure model
fit <- fit.cure.model(Surv(FUyear, status2) ~ 1, data = colonDC, bhazard = "bhaz",
                      formula.k1 = ~ 1, formula.k2 = ~ 1,
                      type = "mixture", dist = "weibull", link = "logit")

#Compute the probability of cancer related death
res <- calc.Crude(fit, time = seq(0, 20, length.out = 100),
                  rmap = list(age = agedays, sex = sex, year = dx))
plot(res)

#Compute the probability of eventually dying from other causes than cancer
res <- calc.Crude(fit, time = seq(0, 20, length.out = 100), type = "othertime",
                  rmap = list(age = agedays, sex = sex, year = dx))
plot(res)


