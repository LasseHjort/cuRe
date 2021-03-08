##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 1000), ]

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", rmap = list(age = "agedays", sex = "sex", year= "dx"),
                            data = colonDC, ratetable = survexp.dk)

#Fit cure model and estimate cure point
fit <- GenFlexCureModel(Surv(FUyear, status) ~ 1, data = colonDC,
                        df = 5, bhazard = "bhaz")
calc.cure.quantile(fit, q = 0.05)
