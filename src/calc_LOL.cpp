//#include <Rcpp.h>
//using namespace Rcpp;
// [[ Rcpp :: depends ( RcppArmadillo )]]



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
// NumericVector predict_model(NumericVector coefs, NumericMatrix X){
//   return exp(-exp(X * coefs));
// }



// m_0 <- function(rel_surv, exp_function, time, tau, extra = FALSE){
//   t_new <- sort(unique(c(time, seq(0, tau, length.out = 5000))), decreasing = T)
//   df_time <- -diff(t_new)
//   vals_pop <- cumsum(c(0, rel_surv(t_new[-length(t_new)]) * exp_function(t_new[-length(t_new)]) * df_time))
//   vals_pop <- rev(vals_pop[t_new %in% time])
//   vals_exp <- cumsum(c(0, exp_function(t_new[-length(t_new)]) * df_time))
//   vals_exp <- rev(vals_exp[t_new %in% time])
//   rmrl <- vals_exp / exp_function(time) - vals_pop / (rel_surv(time) * exp_function(time))
//   if(extra){
//     return(data.frame(LossOfLifetime = rmrl, Time = time))
//   }else{
//     return(data.frame(LossOfLifetime = rmrl / (tau - time), Time = time))
//   }
// }


