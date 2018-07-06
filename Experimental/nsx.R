nsx <- function (x, df = NULL, knots = NULL, intercept = FALSE,
                 Boundary.knots = range(x),
                 derivs = if (cure) c(2,1) else c(2,2),
                 log=FALSE, # deprecated: only used in rstpm2:::stpm2Old
                 centre = FALSE, cure = FALSE, stata.stpm2.compatible=FALSE)
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  if (!missing(df) && missing(knots)) {
    nIknots <- df - 1 - intercept + 4 - sum(derivs)
    if (nIknots < 0) {
      nIknots <- 0
      warning("'df' was too small; have used ", 1 + intercept)
    }
    knots <- if (nIknots > 0) {
      knots <- if (!cure)
        seq.int(0, 1, length.out = nIknots + 2L)[-c(1L,
                                                    nIknots + 2L)]
      else c(seq.int(0, 1, length.out = nIknots + 1L)[-c(1L,
                                                         nIknots + 1L)], 0.95)
      if (!stata.stpm2.compatible)
        stats::quantile(x[!outside], knots)
      else stats::quantile(x[!outside], round(knots,2), type=2)
    }
  }
  else nIknots <- length(knots)
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0,
                                                      1))$design
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0,
                                                      1))$design
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside))
      basis[inside, ] <- spline.des(Aknots, x[inside],
                                    4)$design
  }
  else basis <- spline.des(Aknots, x, 4)$design
  const <- splineDesign(Aknots, rep(Boundary.knots, 3-derivs), 4, c(derivs[1]:2, derivs[2]:2))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  if (centre) {
    centreBasis <- nsx(centre,
                       knots=if (is.null(knots)) numeric(0) else knots,
                       Boundary.knots=Boundary.knots,
                       intercept=intercept, derivs=derivs, centre=FALSE, log=log)
    oldAttributes <- attributes(basis)
    basis <- t(apply(basis,1,function(x) x-centreBasis))
    attributes(basis) <- oldAttributes
  }
  a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots,
            Boundary.knots = Boundary.knots, intercept = intercept, derivs=derivs,
            centre=centre, log=log, q.const=q.const)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("nsx", "basis", "matrix")
  basis
}
