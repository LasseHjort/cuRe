gaussian_quadrature <- function(f, lower, upper, n = 1000){
  nodes <- seq(-1, 1, length.out = n)
  sum_terms <- f((upper - lower) / 2 * nodes + (upper + lower) / 2)
 (upper - lower) / 2 * sum(sum_terms)
}
