##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 1000), ]

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

###Without covariates
##Fit mixture cure model
fit <- FlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC, n.knots = 5, bhazard = "bhaz")

##Plot model
plot(fit)
plot(fit, time = seq(0, 40, length.out = 100))
plot(fit, type = "ehaz")
plot(fit, type = "survuncured")
plot(fit, type = "probcure")

##Predict cure rate
predict(fit, type = "curerate")


##Fit non-mixture cure model
fit <- FlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC, n.knots = 5, bhazard = "bhaz", type = "nmixture")

##Plot relative survival
plot(fit)

##Predict cure rate
predict(fit, type = "curerate")

###With covariates
##Fit mixture cure model
fit <- FlexCureModel(Surv(FUyear, status) ~ sex, data = colonDC, n.knots = 5, bhazard = "bhaz",
                     smooth.formula = ~ sex)

##Plot model
plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), ci = F)
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), col = 2, ci = F, add = TRUE)


plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), ci = FALSE, type = "survuncured")
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0, 15, length.out = 100), col = 2, ci = FALSE, add = TRUE, type = "survuncured")

predict(fit, type = "curerate", data.frame(sex = factor("female", levels = c("male", "female"))))


##Fit mixture cure model with time-varying covariates
fit <- FlexCureModel(Surv(FUyear, status) ~ age, data = colonDC, n.knots = 5, bhazard = "bhaz",
                     n.knots.time = list(age = 3))

##Plot model
plot(fit, newdata = data.frame(age = 70))
plot(fit, newdata = data.frame(age = 60), add = TRUE, col = 2)

plot(fit, type = "ehaz", newdata = data.frame(age = 70), ci = FALSE)
plot(fit, type = "ehaz", newdata = data.frame(age = 60), add = TRUE, col = 2, ci = FALSE)
