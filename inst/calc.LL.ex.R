##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

##Spline-base cure model
#Fit cure model
fit <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz, cure = T)

#Compute and plot the loss of lifetime function
res <- calc.LL(fit, time = seq(0, 20, length.out = 50),
               rmap = list(age = agedays, sex = sex, year = dx))
plot(res)

#Compute and plot the mean residual lifetime
res <- calc.LL(fit, time = seq(0, 20, length.out = 50), type = "mrl",
               rmap = list(age = agedays, sex = sex, year = dx))
plot(res)


#Simple parametric cure model
#Fit cure model
fit <- fit.cure.model(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = "bhaz",
                      type = "mixture", dist = "weibull", link = "logit")

#Compute and plot the loss of lifetime function
res <- calc.LL(fit, time = seq(0, 20, length.out = 50),
               rmap = list(age = agedays, sex = sex, year = dx))
plot(res)

#Compute and plot the mean residual lifetime
res <- calc.LL(fit, time = seq(0, 20, length.out = 50), type = "mrl",
               rmap = list(age = agedays, sex = sex, year = dx))
plot(res)

