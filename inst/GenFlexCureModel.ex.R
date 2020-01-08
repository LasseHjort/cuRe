##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 500), ]

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

##Predict cure proportion
predict(fit, type = "curerate")


##Fit non-mixture cure model
fit <- GenFlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC, df = 4,
                        bhazard = "bhaz", type = "nmixture")

##Plot relative survival
plot(fit)

##Predict cure proportion
predict(fit, type = "curerate")

###With covariates
##Fit mixture cure model
fit <- GenFlexCureModel(Surv(FUyear, status) ~ sex, data = colonDC, df = 4,
                        bhazard = "bhaz", cr.formula = ~ sex)

##Plot model
plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), ci = FALSE)
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), col = 2, ci = FALSE, add = TRUE)


plot(fit, newdata = data.frame(sex = factor("female", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), ci = FALSE, type = "survuncured")
plot(fit, newdata = data.frame(sex = factor("male", levels = c("male", "female"))),
     time = seq(0.001, 15, length.out = 100), col = 2, ci = FALSE,
     add = TRUE, type = "survuncured")

predict(fit, type = "curerate",
        data.frame(sex = factor(c("male", "female"),
                                levels = c("male", "female"))))


##Fit mixture cure model with time-varying covariates
colonDC$gender <- as.numeric(colonDC$sex) - 1
fit <- GenFlexCureModel(Surv(FUyear, status) ~ gender, data = colonDC, df = 6,
                        bhazard = "bhaz", cr.formula = ~ gender, tvc = list(gender = 2))

##Plot model
plot(fit, newdata = data.frame(gender = 0))
plot(fit, newdata = data.frame(gender = 1), add = TRUE, col = 2)

plot(fit, type = "hazard", newdata = data.frame(gender = 0), ci = FALSE)
plot(fit, type = "hazard", newdata = data.frame(gender = 1),
     add = TRUE, col = 2, ci = FALSE)

#Predict cure proportions for a male and female patients
predict(fit, type = "curerate", newdata = data.frame(gender = 0))
predict(fit, type = "curerate", newdata = data.frame(gender = 1))
