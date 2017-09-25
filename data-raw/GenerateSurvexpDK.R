survexp.dk <- transrate.hmd("C:/Users/sw1y/Desktop/Phd/Projekter/Data/DanishLifeTable/mltper_1x1.txt",
                            "C:/Users/sw1y/Desktop/Phd/Projekter/Data/DanishLifeTable/fltper_1x1.txt")
save(survexp.dk, file = "data/survexp.dk.RData")
readRDS("data/survexp.dk.rds")
