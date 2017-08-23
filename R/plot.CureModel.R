plot.CureModel <- function(fit, newdata = NULL, ylim = c(0, 1), xlim = NULL,
                           xlab = "Time", ylab = "Relative survival", non_para = F,
                           type = "relsurv", col = 1, col.non.para = 2, add = F){
  if(length(col) == 1){
    col <- rep(col, nrow(newdata))
  }
  if(type == "probcure"){
    ylab <- "Conditional probability of cure"
  }
  if(is.null(xlim)){
    xlim <- c(0, max(fit$data$FU_years))
  }
  times <- seq(xlim[1], xlim[2], length.out = 100)

  predict_rs <- predict(fit, newdata, times, type = type)
  for(i in 2:ncol(predict_rs)){
    if(i == 2 & !add){
      plot(predict_rs[,i] ~ times, type = "l", ylim = c(0, 1), xlim = xlim, xlab = xlab, ylab = ylab, col = col[i - 1])
    }else{
      lines(predict_rs[,i] ~ times, type = "l", col = col[i - 1])
    }
  }
  if(non_para){
    rsfit <- rs.surv(Surv(FU, status) ~ 1 + ratetable(age = age, sex = sex, year = diag_date),
                     data = fit$data, ratetable = survexp.dk, method = "ederer2")
    rsfit$time <- rsfit$time / year
    lines(rsfit$surv ~ rsfit$time, type = "s", col = col.non.para)
  }
}
