##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)


###Without covariates
##Fit weibull mixture cure model
fit.wei <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                          type = "mixture", dist = "weibull", link = "logit")

##Plot various summaries of the model (see ?predict.cm)
plot(fit.wei)
plot(fit.wei, time = seq(0, 40, length.out = 100))
plot(fit.wei, type = "hazard")
plot(fit.wei, type = "survuncured")
plot(fit.wei, type = "probcure")

#Fit a weibull-weibull mixture cure model
fit.weiwei <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                          type = "mixture", dist = "weiwei", link = "logit")

#Compare to the weibull model
plot(fit.wei, ci = F)
plot(fit.weiwei, add = T, col = 2, ci = F)

###With covariates
##Fit weibull mixture cure model with age effect for both components of the Weibull model
fit <- fit.cure.model(Surv(FUyear, status) ~ age, data = colonDC, bhazard = "bhaz",
                      formula.surv = list(~ age, ~ age),
                      type = "mixture", dist = "weibull", link = "logit")

##Plot model for age 50 and 60
plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), ci = F)
plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), ci = F, add = TRUE, col = 2)

plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), ci = F, type = "hazard")
plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), ci = F, type = "hazard", add = T, col = 2)
