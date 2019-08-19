##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 1000), ]

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

##Fit flexible parametric relative survival model
fit <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz)

##Compute survival probabilities from 0 to 20 years
pred <- lts(fit, rmap = list(age = agedays, sex = sex, year = dx))

##Plot the survival function
plot(pred)


