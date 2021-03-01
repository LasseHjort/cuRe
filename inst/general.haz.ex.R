##Use data cleaned version of the colon cancer data from the rstpm2 package
data("colonDC")
set.seed(2)
colonDC <- colonDC[sample(1:nrow(colonDC), 1000), ]

##Extract general population hazards
bhaz1 <- general.haz(time = "FU",
                     rmap = list(age = "agedays", sex = "sex", year= "dx"),
                     data = colonDC,
                     ratetable = survexp.dk)

bhaz2 <- general.haz(time = colonDC$FU,
                     rmap = list(age = colonDC$agedays, sex = colonDC$sex, year = colonDC$dx),
                     data = colonDC,
                     ratetable = survexp.dk)

all(bhaz2 == bhaz1)
