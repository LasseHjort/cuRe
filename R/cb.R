#' Restricted cubic splines
#'
#' Function for computing the basis matrix for restricted cubic splines
#'
#' @param x Values to evaluate the basis functions in.
#' @param knots Chosen knots for the spline.
#' @param ortho Logical. If \code{TRUE} (default) orthogonalization of the basis matrix is carried out.
#' @param R.inv Matrix or vector containing the values of the R matrix from the QR decomposition of the basis matrix.
#' This is used for making new predictions based on the initial orthogonalization.
#' Therefore the default is \code{NULL}.
#' @param intercept Logical. If \code{FALSE}, the intercept of the restricted cubic spline is removed.
#' @return A matrix with containing the basis functions evaluated in \code{x}.
#' @references Royston P. and Parmar M.K. (2002) Flexible parametric proportional-hazards and proportional-odds
#' models for censored survival data, with application to prognostic modelling and
#' estimation of treatment effects. \emph{Statistics in Medicine}, 21:15.
#' @export

cb <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
  nx <- length(x)
  if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
  nk <- ncol(knots)
  b <- matrix(nrow=length(x), ncol=nk)
  if (nk>0){
    b[,1] <- 1
    b[,2] <- x
  }
  if (nk>2) {
    lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
    for (j in 1:(nk-2)) {
      b[,j+2] <- pmax(x - knots[,j+1], 0)^3 - lam[,j+1]*pmax(x - knots[,1], 0)^3 -
        (1 - lam[,j+1])*pmax(x - knots[,nk], 0)^3
    }
  }

  if(!intercept) b <- b[,-1, drop = FALSE]

  if(ortho){
    if(is.null(R.inv)){
      qr_decom <- qr(b)
      b <- qr.Q(qr_decom)
      R.inv <- solve(qr.R(qr_decom))
    } else{
      R.inv <- matrix(R.inv, nrow = ncol(b))
      b <- b %*% R.inv
    }
  } else {
    R.inv <- diag(ncol(b))
  }

  a <- list(knots = knots[1,], ortho = ortho, R.inv = R.inv, intercept = intercept)
  attributes(b) <- c(attributes(b), a)
  class(b) <- c("cb", "matrix")
  b
}

#Predict function associated with bsx.
predict.cb <- function (object, newx, ...)
{
  if (missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "ortho",
                                              "R.inv", "intercept")])
  do.call("cb", a)
}

#Additional function needed to fix the knot location in cases where df is only specified
#' @export
makepredictcall.cb <- function (var, call)
{
  if (as.character(call)[1L] != "cb")
    return(call)
  at <- attributes(var)[c("knots", "ortho", "R.inv", "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

#Derivate of basis function
dbasis <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE) {
  if(ortho & is.null(R.inv)) stop("Both 'ortho' and 'R.inv' has to be specified!")
  nx <- length(x)
  if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
  nk <- ncol(knots)
  b <- matrix(nrow=length(x), ncol=nk)
  if (nk>0){
    b[,1] <- 0
    b[,2] <- 1
  }
  if (nk>2) {
    lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
    for (j in 3:nk) {
      b[,j] <- 3*pmax(x - knots[,j-1], 0)^2 - 3*lam[,j-1]*pmax(x - knots[,1], 0)^2 -
        3*(1 - lam[,j-1])*pmax(x - knots[,nk], 0)^2
    }
  }

  if(!intercept) b <- b[, -1]
  if(ortho){
    b <- b %*% R.inv
  }

  b
}
