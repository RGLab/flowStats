#' The version of fdPar from fda 2.4.0 because the new API changes the output.
#' (specifically res$fd$coefs) and thus breaks the landmarkreg call.
#' @name fdPar
#' @md
#' @param fdobj functional data object, functional basis object, a functional parameter object or a matrix. If it a matrix, it is replaced by fd(fdobj). If class(fdobj) == 'basisfd', it is converted to an object of class fd with a coefficient matrix consisting of a single column of zeros.
#' @param Lfdobj either a nonnegative integer or a linear differential operator object.
#' If NULL, Lfdobj depends on `fdobj[['basis']][['type']]`
#' - bspline `Lfdobj <- int2Lfd(max(0, norder-2))`, where `norder = norder(fdobj)`
#' - fourier Lfdobj = a harmonic acceleration operator:
#'  `Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)`
#'  where rng = `fdobj[['basis']][['rangeval']]`.
#' - anything else `Lfdobj <- int2Lfd(0)`
#' @param lambda a nonnegative real number specifying the amount of smoothing to be applied to the estimated functional parameter.
#' @param estimate not currently used.
#' @param penmat a roughness penalty matrix. Including this can eliminate the need to compute this matrix over and over again in some types of calculations.
#' @export
fdPar <- function (fdobj = NULL, Lfdobj = NULL, lambda = 0, estimate = TRUE, 
    penmat = NULL) 
{
  if (!inherits(fdobj, "fd")) {
    if (is.null(fdobj)) {
      fdobj = fd()
    }
    else {
      if (inherits(fdobj, "basisfd")) {
        nbasis <- fdobj$nbasis
        dropind <- fdobj$dropind
        coefs <- matrix(0, nbasis - length(dropind), 
            1)
        fdnames <- list("time", "reps 1", "values")
        if (!is.null(fdobj$names)) {
          basisnames <- {
            if (length(dropind) > 0) 
              fdobj$names[-dropind]
            else fdobj$names
          }
          dimnames(coefs) <- list(basisnames, NULL)
          fdnames[[1]] <- basisnames
        }
        fdobj <- fd(coefs, fdobj, fdnames)
      }
      else if (is.numeric(fdobj)) 
        fdobj <- fd(fdobj)
      else stop("First argument is neither a functional data object ", 
            "nor a basis object.")
    }
  }
  {
    if (is.null(Lfdobj)) {
      if (fdobj$basis$type == "fourier") {
        rng <- fdobj$basis$rangeval
        Lfdobj <- vec2Lfd(c(0, (2 * pi/diff(rng))^2, 
                0), rng)
      }
      else {
        norder <- {
          if (fdobj$basis$type == "bspline") 
            norder.bspline(fdobj$basis)
          else 2
        }
        Lfdobj <- int2Lfd(max(0, norder - 2))
      }
    }
    else Lfdobj <- int2Lfd(Lfdobj)
  }
  if (!inherits(Lfdobj, "Lfd")) 
    stop("'Lfdobj' is not a linear differential operator object.")
  if (!is.numeric(lambda)) 
    stop("Class of LAMBDA is not numeric.")
  if (lambda < 0) 
    stop("LAMBDA is negative.")
  if (!is.logical(estimate)) 
    stop("Class of ESTIMATE is not logical.")
  if (!is.null(penmat)) {
    if (!is.numeric(penmat)) 
      stop("PENMAT is not numeric.")
    penmatsize <- dim(penmat)
    if (any(penmatsize != nbasis)) 
      stop("Dimensions of PENMAT are not correct.")
  }
  fdParobj <- list(fd = fdobj, Lfd = Lfdobj, lambda = lambda, 
      estimate = estimate, penmat = penmat)
  oldClass(fdParobj) <- "fdPar"
  fdParobj
}