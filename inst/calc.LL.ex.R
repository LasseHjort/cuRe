###Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")

##Extract general population hazards
colonDC$bhaz <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                            data = colonDC, ratetable = survexp.dk)

##Fit flexible parametric relative survival model
fit <- stpm2(Surv(FUyear, status2) ~ 1, data = colonDC, df = 6, bhazard = colonDC$bhaz)

#Compute loss of lifetime from 0 to 20 years
res <- calc.LL(fit, time = seq(0, 20, length.out = 100),
               rmap = list(age = agedays, sex = sex, year = dx))

#Plot the loss of lifetime
plot(res)

##Fit flexible parametric cure model
fit <- FlexCureModel(Surv(FUyear, status2) ~ 1, data = colonDC, n.knots = 7, bhazard = "bhaz",
                     ini.types = "cure")

#Compute loss of lifetime from 0 to 20 years
res.cure <- calc.LL(fit, time = seq(0, 20, length.out = 100),
                    rmap = list(age = agedays, sex = sex, year = dx))

#Plot the loss of lifetime
plot(res, ci = F, ylim = c(0, max(res$Ests[[1]]$ll)))
plot(res.cure, add = T, col = 2, ci = F)


#Compute mean residual lifetime
res.cure <- calc.LL(fit, time = seq(0, 20, length.out = 100), type = "mrl",
                    rmap = list(age = agedays, sex = sex, year = dx))
plot(res.cure)

