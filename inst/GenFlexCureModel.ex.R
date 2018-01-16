##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

###Without covariates
##Fit mixture cure model
fit <- GenFlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC, df = 4, bhazard = "bhaz")

##Plot model
plot(fit)
plot(fit, time = seq(0.001, 40, length.out = 100))
plot(fit, type = "hazard")
plot(fit, type = "survuncured")
plot(fit, type = "probcure")

##Predict cure rate
predict(fit, type = "curerate")


##Fit non-mixture cure model
fit <- GenFlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC, df = 4,
                        bhazard = "bhaz", type = "nmixture")

##Plot relative survival
plot(fit)

##Predict cure rate
predict(fit, type = "curerate")

###With covariates
##Fit mixture cure model
fit <- GenFlexCureModel(Surv(FUyear, status) ~ sex, data = colonDC, df = 4, bhazard = "bhaz", cr.formula = ~ sex)

##Plot model
plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), ci = "n")
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), col = 2, ci = "n", add = T)


plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), ci = "n", type = "survuncured")
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), col = 2, ci = "n", add = T, type = "survuncured")

predict(fit, type = "curerate", data.frame(sex = factor("female", levels = c("male", "female"))))


##Fit mixture cure model with time-varying covariates
fit <- GenFlexCureModel(Surv(FUyear, status) ~ age, data = colonDC, df = 4, bhazard = "bhaz",
                        cr.formula = ~ age, tvc = list(age = 2))

##Plot model
plot(fit, newdata = data.frame(age = 70))
plot(fit, newdata = data.frame(age = 60), add = T, col = 2)

plot(fit, type = "hazard", newdata = data.frame(age = 70), ci = "n")
plot(fit, type = "hazard", newdata = data.frame(age = 60), add = T, col = 2, ci = "n")
