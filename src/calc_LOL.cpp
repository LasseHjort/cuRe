#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo
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
NumericVector calc_area(NumericVector knots, NumericVector coefs, int n.x = 100, double lower = 0, double tau) {
  NumericVector x(n.x);
  x[0] = -1;
  int step.size = 2 / n.x;
  for(int i = 1; i < n.x; i++){
    x[i] = x[i - 1] + step.size;
  }
  NumericMatrix b = basis(knots, log((tau - lower) / 2 * x + (tau + lower) / 2));
  NumericVector lp = b coefs;
  double res = (tau - lower) / 2 * sum(lp);
  NumericMatrix b = basis(knots, log(x));
  NumericVector RS = b %*%
 return x * 2;
}

NumericVector calc_int(f)




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



// [[Rcpp::export]]
NumericMatrix basis(NumericVector knots, NumericVector x){
  double nx = x.size();
  int nk = knots.size();
  NumericMatrix b(nx, nk);
  if(nk > 0){
    NumericVector unit(nx, 1);
    b(_, 0) = unit;
    b(_, 1) = x;
  }
  if(nk > 2){
    NumericVector lam = (knots[nk - 1] - knots) / (knots[nk - 1] - knots[0]);
    for(int i = 0; i < (nk - 2); i++){
      NumericVector term1 = pow(pmax(x - knots[i + 1], 0), 3);
      NumericVector term2 = lam[i + 1] * pow(pmax(x - knots[0], 0), 3);
      NumericVector term3 = (1 - lam[i + 1]) * pow(pmax(x - knots[nk - 1], 0), 3);
      b(_,i + 2) = term1 - term2 - term3;
    }
  }
  return b;
}

