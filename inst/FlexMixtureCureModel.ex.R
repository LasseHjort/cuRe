##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

####Without covariates
##Fit relative survival model
fit <- FlexMixtureCureModel(Surv(FUyear, status2) ~ 1, data = colonDC, n.knots = 6, bhazard = "bhaz",
                            ini.types = "cure")

##Plot model
plot(fit)
plot(fit, time = seq(0, 40, length.out = 100))
plot(fit, type = "ehaz")
plot(fit, type = "survuncured")
plot(fit, type = "probcure")


####With covariates
##Fit relative survival model
fit <- FlexMixtureCureModel(Surv(FUyear, status2) ~ sex, data = colonDC, n.knots = 6, bhazard = "bhaz",
                            ini.types = "cure", smooth.formula = ~ sex)

##Plot model
plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), ci = F)
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), col = 2, ci = F, add = T)


plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), ci = F, type = "survuncured")
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), col = 2, ci = F, add = T, type = "survuncured")

#Predict cure rate for female patients
predict(fit, type = "curerate", newdata = data.frame(sex = factor("female", levels = c("male", "female"))))
