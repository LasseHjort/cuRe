##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

#Fit cure model and estimate cure point
fit <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz, cure = T)
cp <- calc.Crude.quantile(fit, q = 0.05, rmap = list(age = agedays, sex = sex, year = dx))

#Compare the result with the trajectory of the conditional probability of disease-related death
res <- calc.Crude(fit, type = "condother", time = seq(0, 20, length.out = 100),
                  var.type = "n", rmap = list(age = agedays, sex = sex, year = dx), reverse = T)
plot(res)
abline(h = 0.05, v = cp$Estimate)

