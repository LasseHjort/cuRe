bhaz1 <- general.haz(time = "FU", age = "agedays", sex = "sex", year = "dx",
                     data = colonDC, ratetable = survexp.dk)

bhaz2 <- general.haz(time = colonDC$FU, age = colonDC$agedays, sex = colonDC$sex,
                     year = colonDC$dx, ratetable = survexp.dk)
all(bhaz2 == bhaz1)
