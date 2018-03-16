##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)


###Without covariates
##Fit weibull mixture cure model
fit.wei <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                          type = "mixture", dist = "weibull", link = "logit")

##Plot various summaries of the model
plot(fit.wei)
plot(fit.wei, time = seq(0, 40, length.out = 100))
plot(fit.wei, type = "hazard")
plot(fit.wei, type = "survuncured")
plot(fit.wei, type = "probcure")

#Fit a weibull-weibull mixture cure model
fit.weiwei <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                          type = "mixture", dist = "weiwei", link = "logit")

#Compare to the weibull model
plot(fit.wei, var.type = "n")
plot(fit.weiwei, add = T, col = 2, var.type = "n")

###With covariates
##Fit weibull mixture cure model
fit <- fit.cure.model(Surv(FUyear, status) ~ age, data = colonDC, bhazard = "bhaz",
                      formula.surv = list(~ age, ~1),
                      type = "mixture", dist = "weibull", link = "logit")

##Plot model
plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), var.type = "n")
plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), var.type = "n", add = TRUE, col = 2)

plot(fit, newdata = data.frame(age = 60),
     time = seq(0, 15, length.out = 100), var.type = "n", type = "hazard")
plot(fit, newdata = data.frame(age = 50),
     time = seq(0, 15, length.out = 100), var.type = "n", type = "hazard", add = T, col = 2)
