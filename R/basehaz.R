
get_basehaz <- function(fit){
  id <- 1:length(fit$y[,1])
  od <- order(fit$y[,1])
  id_new <- id[od]
  fu <- fit$y[od, 1]
  status <- fit$y[od, 2]
  death_times <- fu[status == 1]
  uni_death_times <- unique(death_times)
  lp <- exp(fit$linear.predictors[od])
  h <- rep(NA, length(uni_death_times))
  for(i in 1:length(uni_death_times)){
    d_i <- length(which(death_times == uni_death_times[i]))
    sum_rj <- sum(lp[fu >= uni_death_times[i]])
    h[i] <- d_i / sum_rj
  }
  H <- c(0, cumsum(h))
  uni_death_times <- c(0, uni_death_times)
  wh <- findInterval(fu, vec = uni_death_times)
  H_all <- H[wh]
  D <- data.frame(hazard = H_all, time = fu, ID = id_new)
  D[order(D$ID), -3]
}

