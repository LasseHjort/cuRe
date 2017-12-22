#Use colon cancer dataset from the rstpm2 package
library(rstpm2)
colonDC <- rstpm2::colon

#Add variables
colonDC$sex <- factor(colonDC$sex, levels = c("Female", "Male"), labels = c("female", "male"))
colnames(colonDC)[colnames(colonDC) == "status"] <- "statusCOD"
colonDC$status <- ifelse(colonDC$statusCOD %in% c("Dead: cancer", "Dead: other"), 1, 0)
colonDC$FU <- as.numeric(colonDC$exit - colonDC$dx)
colonDC$FUyear <- colonDC$FU / 365.24
colonDC$agedays <- colonDC$age * 365.24
colonDC <- subset(colonDC, select = -c(mmdx, yydx, surv_mm, surv_yy, year8594, agegrp))

#Save as .RData object
save(colonDC, file = "data/colonDC.RData")


