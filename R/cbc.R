###############KIG HER LASSSSSEEEEEE##################3
#Fix problem med cb.cure! Det virker med orthogonalisering.

#' Restricted cubic splines with cure
#'
#' Function for computing the basis matrix for restricted cubic splines which are constant beyond the last knot
#'
#' @param x Values to evaluate the basis functions in.
#' @param knots Chosen knots for the spline.
#' @param ortho Logical. If \code{TRUE} (default) orthogonalization of the basis matrix is carried out.
#' @param R.inv Matrix or vector containing the values of the R matrix from the QR decomposition of the basis matrix.
#' This is used for making new predictions based on the initial orthogonalization.
#' Therefore the default is \code{NULL}.
#' @param intercept Logical. If \code{FALSE}, the intercept of the restricted cubic spline is removed.
#' @return A matrix with containing the basis functions evaluated in \code{x}.
#' @references Andersson T.M.-L., et al. (2011) Estimating and modelling cure in population-based cancer
#' studies within the framework of flexible parametric survival models.
#' \emph{BMC Medical Research Methodology}, 11:96.
#' @export
cbc <- function(x, knots, ortho = TRUE, R.inv = NULL, intercept = TRUE){
  nk <- length(knots)
  b <- matrix(nrow = length(x), ncol = nk - 1)
  knots_rev <- rev(knots)
  if (nk > 0) {
    b[, 1] <- 1
  }
  if (nk > 2) {
    for (j in 2:(nk - 1)) {
      lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
      b[, j] <- pmax(knots_rev[nk - j + 1] - x, 0)^3 - lam * pmax(knots_rev[nk] - x, 0)^3 -
        (1 - lam) * pmax(knots_rev[1] - x, 0)^3
    }
  }

  if(!intercept) b <- b[,-1, drop = F]

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

  a <- list(knots = knots, ortho = ortho, R.inv = R.inv, intercept = intercept)
  attributes(b) <- c(attributes(b), a)
  class(b) <- c("cbc", "matrix")
  b
}

#Predict function associated with bsx.
predict.cbc <- function (object, newx, ...)
{
  if (missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "ortho",
                                              "R.inv", "intercept")])
  do.call("cbc", a)
}

#Additional function needed to fix the knot location in cases where df is only specified
#' @export
makepredictcall.cbc <- function (var, call)
{
  if (as.character(call)[1L] != "cbc")
    return(call)
  at <- attributes(var)[c("knots", "ortho", "R.inv", "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

dbasis.cure <- function(knots, x, ortho = TRUE, R.inv = NULL, intercept = TRUE){
  nk <- length(knots)
  b <- matrix(nrow = length(x), ncol = nk - 1)
  knots_rev <- rev(knots)
  if (nk > 0) {
    b[, 1] <- 0
  }
  if (nk > 2) {
    for (j in 2:(nk - 1)) {
      lam <- (knots_rev[nk - j + 1] - knots_rev[1])/(knots_rev[nk] - knots_rev[1])
      b[, j] <- - 3 * pmax(knots_rev[nk - j + 1] - x, 0)^2 + 3 * lam * pmax(knots_rev[nk] - x, 0)^2 +
        3 * (1 - lam) * pmax(knots_rev[1] - x, 0)^2
    }
  }

  if(!intercept) b <- b[, -1]

  if(ortho){
    b <- b %*% R.inv
  }
  b
}
