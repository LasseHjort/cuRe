########Extensions of polynomial splines to incorporate user-defined trajectories beyond the boundary knots.
#Create by: Lasse Hjort Jakobsen
#Date: 2018-05-18

#The bs function is extended to bsx by user a QR decompositon to restrict the trajectory beyond the boundary knots.
#The restriction is possible for both the first and last knot.


#' Polynomial B-splines with eXtensions
#'
#' Generate the B-spline basis matrix for a polynomial spline with derivative restrictions at the boundary knots.
#'
#' @param x the predictor variable. Missing values are allowed.
#' @param df degrees of freedom; one can specify \code{df} rather than knots; \code{bs()} then chooses
#' \code{df}-\code{degree} (minus one if there is an intercept) knots at suitable quantiles of \code{x}
#' (which will ignore missing values). The default, \code{NULL}, corresponds to no inner knots,
#' i.e., \code{degree}-\code{intercept}.
#' @param knots the internal breakpoints that define the spline. The default is \code{NULL}, which results
#' in a basis for ordinary polynomial regression. Typical values are the mean or median for one knot,
#' quantiles for more knots. See also \code{Boundary.knots}.
#' @param degree degree of the piecewise polynomial—default is \code{3} for cubic splines.
#' @param Boundary.knots boundary points at which to anchor the B-spline basis (default the range of the non-NA data).
#' If both \code{knots} and \code{Boundary.knots} are supplied, the basis parameters do not depend on \code{x}.
#' Data can extend beyond \code{Boundary.knots}.
#' @param intercept if \code{TRUE}, an intercept is included in the basis; default is \code{FALSE}.
#' @param deriv an integer vector of length 2 with values between 0 and \code{degree + 1} giving the
#' derivative constraint order at the left and right boundary knots;
#' an order of 2 constrains the second derivative to zero (f”(x)=0);
#' an order of 1 constrains the first and second derivatives to zero (f'(x)=f”(x)=0);
#' an order of 0 constrains the zero, first and second derivatives to zero (f(x)=f'(x)=f”(x)=0)
#' An order of \code{degree + 1} computes the basis matrix similarly to \code{bs}.
#' @return A matrix with containing the basis functions evaluated in \code{x}.
bsx <- function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                 Boundary.knots = range(x), deriv = NULL)
{
  if(is.null(deriv)){
    deriv <- rep(degree + 1, 2)
  }
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1)
    stop("'degree' must be integer >= 1")
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  } else FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept) + if(all(ord == deriv)) 0 else sum(ord - deriv - c(0, 1))
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                         2L)[-c(1L, nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    #warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree -
                          1L))
    e <- 1/4
    if (any(ol)) {
      #k.pivot <- (1 - e) * Boundary.knots[1L] + e * Aknots[ord +
      #                                                       1]
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      #k.pivot <- (1 - e) * Boundary.knots[2L] + e * Aknots[length(Aknots) -
      #                                                       ord]
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
                           "^"))
      #xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree - 1,
      #                     "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs)
      basis[or, ] <- xr %*% (tt/scalef)
      #basis[or, ] <- splineDesign(Aknots, rep(k.pivot, length(x[or])), ord,
      #                            derivs = 0)
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
                                      ord)
  } else basis <- splineDesign(Aknots, x, ord)

  if(!intercept){
    basis <- basis[, -1L, drop = FALSE]
  }

  x.tmp <- rep(Boundary.knots, if(all(ord == deriv)) c(0,0) else ord - deriv - c(0, 1))
  deriv1 <- if(deriv[1] < ord) deriv[1]:(ord - 1) else NULL
  deriv2 <- if(deriv[2] < ord - 1) deriv[2]:(ord - 2) else NULL
  deriv.tmp <- c(deriv1, deriv2)

  if(!is.null(deriv.tmp)){
    const <- splineDesign(knots = Aknots, x = x.tmp, ord = ord, deriv = deriv.tmp)

    if (!intercept){
      const <- const[, -1L, drop = FALSE]
    }
    qr.const <- qr(t(const))
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:nrow(const)), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), drop = FALSE])
  } else q.const <- NULL

  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
            Boundary.knots = Boundary.knots, intercept = intercept, deriv = deriv,
            q.const = q.const)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bsx", "basis", "matrix")
  basis
}

#Predict function associated with bsx.
predict.bsx <- function (object, newx, ...)
{
  if (missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots",
                                              "intercept", "deriv", "degree")])
  do.call("bsx", a)
}

#Additional function needed to fix the knot location in cases where df is only specified
makepredictcall.bsx <- function (var, call)
{
  if (as.character(call)[1L] != "bsx")
    return(call)
  at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                          "deriv", "degree")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}


#Prespecified arguments for testing
library(splines)
x <- -20:20
df <- NULL
knots <- c(2, 5,7)
Boundary.knots <- c(1, 10)
intercept <- FALSE
deriv <- NULL
degree <- 3


# #Lets test if it works for cubic splines. Define x values and boundary knots
x <- seq(-1, 14, length.out = 100)
Boundary.knots <- c(1, 10)
#Compute basis
b <- bsx(x = x, df = 5, degree = 3, deriv = c(4,0), Boundary.knots = Boundary.knots)
#Pick random coefficients
coefs <- rnorm(ncol(b))
#Plot the trajectory
plot(x, b %*% coefs, type = "l")
abline(v = Boundary.knots)
#
#
# #Lets try degree 4 and 5
# b <- bsx(x = x, df = 5, degree = 4, deriv = c(3,1), Boundary.knots = Boundary.knots)
# plot(time, b %*% coefs, type = "l")
# abline(v = Boundary.knots)
#
# b <- bsx(x = x, df = 5, degree = 5, deriv = c(2,1), Boundary.knots = Boundary.knots)
# plot(time, b %*% coefs, type = "l")
# abline(v = Boundary.knots)
#
# #Lets try degree 3 and no restrictions beyond the first knot
# b <- bsx(x = x, df = 5, degree = 3, deriv = c(4,1), Boundary.knots = Boundary.knots)
# plot(time, b %*% coefs, type = "l")
# abline(v = Boundary.knots)
#
#
# #So seems to work. Lets try it in rstpm2::stpm2.
# library(rstpm2)
# load(file = "../TestScripts/colonDC.RData")
#
# fit <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz, df = 5, cure = T)
# plot(fit, newdata = data.frame(age = 50), ylim = c(0, 1))
#
# fit2 <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 5, deriv = c(2, 1)))
# plot(fit2, newdata = data.frame(age = 50), add = T, line.col = 2)
#
# fit3 <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 5, deriv = c(2, 1), degree = 4))
# plot(fit3, newdata = data.frame(age = 50), add = T, line.col = 3)
#
# fit4 <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 5, deriv = c(2, 1), degree = 5))
# plot(fit4, newdata = data.frame(age = 50), add = T, line.col = 4)
# #It is a bit weird that the cure rate seems to increase as I increase the polynomial degree.
# #However, I cannot immediately explain this
#
# #Lets try to adjust the derivative constraints
# fit <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 3, deriv = c(2, 1)))
# plot(fit, newdata = data.frame(age = 50), ylim = c(0, 1))
#
# fit2 <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 3, deriv = c(3,3)))
# plot(fit2, newdata = data.frame(age = 50), add = T, line.col = 2)
#
# fit3 <- stpm2(Surv(FUyear, status) ~ 1, data = colonDC, bhazard = bhaz,
#               smooth.formula = ~bsx(x = log(FUyear), df = 5))
# plot(fit3, newdata = data.frame(age = 50), add = T, line.col = 3)

#Does seem to provide decent results

