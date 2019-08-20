
#Use the transrate.hmd function of life tables downloaded from
#the human mortality database seperately for men and women
survexp.dk <- transrate.hmd("mltper_1x1.txt",
                            "fltper_1x1.txt")

#Save as .RData object
save(survexp.dk, file = "data/survexp.dk.RData")
