##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)


###Without covariates
##Fit weibull mixture cure model
fit <- fit.cure.model(Surv(FUyear, status2) ~ 1, data = colonDC, bhazard = "bhaz",
                      formula.k1 = ~ 1, formula.k2 = ~ 1,
                      type = "mixture", dist = "weibull", link = "logit")

##Plot model
plot(fit)
plot(fit, time = seq(0, 40, length.out = 100))
plot(fit, type = "ehaz")
plot(fit, type = "survuncured")
plot(fit, type = "probcure")


###With covariates
##Fit weibull mixture cure model
fit <- fit.cure.model(Surv(FUyear, status2) ~ age, data = colonDC, bhazard = "bhaz",
                      formula.k1 = ~ age, formula.k2 = ~ 1,
                      type = "mixture", dist = "weibull", link = "logit")

##Plot model
plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), ci = F)
plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), col = 2, ci = F, add = T)

plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), ci = F, type = "ehaz")
plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), col = 2, ci = F, add = T, type = "ehaz")
