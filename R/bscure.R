x <- -20:20
df <- NULL
knots <- c(2, 5,7)
Boundary.knots <- c(1, 10)
intercept <- FALSE
deriv <- c(2,1)
degree <- 3

bsx <- function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                 Boundary.knots = range(x), deriv = c(2, 2))
{
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
  }
  else FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
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
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree -
                          1L))
    e <- 1/4
    if (any(ol)) {
      k.pivot <- (1 - e) * Boundary.knots[1L] + e * Aknots[ord +
                                                             1]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- (1 - e) * Boundary.knots[2L] + e * Aknots[length(Aknots) -
                                                             ord]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
                         derivs)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
                                      ord)
  }
  else basis <- splineDesign(Aknots, x, ord)
  #basis <- splineDesign(Aknots, x, ord, outer.ok = T)
  const <- splineDesign(knots = Aknots, x = rep(Boundary.knots, 3 - deriv),
                        ord = ord, deriv = c(deriv[1]:(ord - 2), deriv[2]:(ord - 2)))
  if (!intercept){
    basis <- basis[, -1L, drop = FALSE]
    const <- const[, -1L, drop = FALSE]
  }
  qr.const <- qr(t(const))
  q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), drop = FALSE])

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


time <- -20:20
a <- bsx(x = time, knots = c(2,5,7), Boundary.knots = c(1, 10), deriv = c(2,1))

plot(time, a %*% 1:4, type = "l")

debug(nsx)
a <- nsx(x = time, knots = c(2, 5,7), Boundary.knots = c(1, 10), deriv = c(2,1))
plot(time, a %*% 1:3, type = "l")

knots <- log(seq(0, 20, length.out = 7))

func <- function(time){
  b <- bs(x = time, knots = exp(knots[-c(1, length(knots))]), Boundary.knots = exp(range(knots)))
  exp(b %*% c(5, -5, -5, -5, 5, -5, -15, 10))
}

gausswd <- statmod::gauss.quad(300)

surv_fun <- function(time){
  res <- rep(NA, length(time))
  for(i in 1:length(time)){
    res[i] <-   exp(-time[i] / 2 * sum(gausswd$weights * func(time[i] / 2 * gausswd$nodes + time[i] / 2)))
  }
  res
}


curve(surv_fun, from = 0, to = 21)
surv_fun(15)

func(20)


debug(nsx)
nsx(x = log(times), knots = knots[-c(1, length(knots))], Boundary.knots = range(knots), cure = T)
b <- bs(x = log(times), knots = knots[-c(1, length(knots))], Boundary.knots = range(knots))

plot(times, exp(b %*% c(0.1, 0.1, 0.1, 0.1, -0.1, 0.1, -0.1, 0.1)))


a <- function(x) splineDesign(knots = c(1,2,3,4,5, 6, 7, 8), ord = 3, x = x) %*% c(0.5,0.5, 0.1, 0.5, 0.5)
curve(a, from = 3, to = 6, n = 100)




