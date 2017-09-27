##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

##Fit flexible parametric relative survival model
fit <- stpm2(Surv(FUyear, status2) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz)

##Compute survival probabilities from 0 to 20 years
pred <- predict.lts(fit, time = seq(0, 20, length.out = 100),
                    rmap = list(age = agedays, sex = sex, year = dx))

##Plot the survival function
plot(pred)

##Fit flexible parametric cure model
fit <- FlexMixtureCureModel(Surv(FUyear, status2) ~ 1, data = colonDC, n.knots = 7, bhazard = "bhaz",
                            ini.types = "cure")

##Compute survival probabilities from 0 to 20 years
pred.cure <- predict.lts(fit, time = seq(0, 20, length.out = 100),
                         rmap = list(age = agedays, sex = sex, year = dx))

##Plot the loss of lifetime
plot(pred, ci = F)
plot(pred.cure, add = T, col = 2, ci = F)


