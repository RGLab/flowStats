\name{iProcrustes}

\alias{iProcrustes}

\title{Procrustes analysis. Using singular value decomposition (SVD) to
       determine a linear transformation of a matrix X to align the
       points in X to the points in a reference matrix Y.
}

\description{
  As a part of generalized Procrustes analysis, this function determines a
  linear transformation (rotation/reflection and scalling) of the points in
  matrix \code{x} to align them to their reference points in
  \code{xbar}. The alignemnt is carried out by
  minimizing the distance between the points in \code{x} and \code{xbar}.
}

\usage{
iProcrustes(x, xbar, rotation.only=FALSE, scalling=TRUE, translate=FALSE)
}

\arguments{
  \item{x}{A numerical matrix. The columns represents the coordinates
  of the points. The matrices \code{x} and \code{xbar} must have the
  same dimensions.} 

  \item{xbar}{A numerical, reference matrix to which points in matrix \code{x}
  are to be aligned.}
 
  \item{rotation.only}{Logical. When \code{rotaion.only} is TRUE, it allows the
  function to lose reflection component of the linear transformation. Although
  it might not give the best-fitting aligenment, in the case 
  of dealing with flow cytometry data, reflection transformation is
  not desirable. When \code{rotaion.only} is FALSE, it allows the function to
  retain the reflection component.}

  \item{scalling}{Logical. When \code{scalling} is FALSE, it allows the
    function to calculate the linear transformation without a scalling factor.
    That is, the returning scalling factor is set to \eqn{1}.} 

  \item{translate}{Logical. Set it to FALSE when the points in matrices
    x and xbar are already centralized prior to applying this function.
    When \code{translate} is TRUE, it allows the
    function to translate the centroid the points in matrices \code{x} to that
    of points in \code{xbar}.}
}

\details{
  Suppose the points in matrix \eqn{X} and \eqn{\bar{X}} are centralized
  (meaning their centroids are at the origin). The
  linear transformation of \eqn{X} for aligning \eqn{X} to its reference
  matrix \eqn{\bar{X}}., i.e., min \eqn{||sXQ - \bar{X}||_F}, is given by:
    \deqn{Q = VU^T,}
  and
    \deqn{s = trace(\bar{X}^TXQ) / trace(X^T X),}
  where V and U are the sigular value vectors of \eqn{\bar{X}^T X} (that is, 
  \eqn{\bar{X}^T X = U \Sigma V^T}), and \eqn{s} is the scalling factor. 

}

\value{
  A list of the linear tranformation with items
  \item{Q}{An orthogonal, rotation/reflection matrix.}
  \item{scal}{A scalling factor}
  \item{T}{(optional) A translation vector used to shift the centroid of the
  points in matrix \code{x} to the origin. Returned when \code{translate}
  is TRUE.} 
  \item{T.xbar}{(optional) Centered \code{xbar} (that is, the centroid
  of the points in \code{xbar} is translated to the origin). Returned when
  \code{translate} is TRUE.}

  Note that the return values of this function do not include the transformed
  matrix \eqn{scal* x* Q} or \eqn{scal*(x-IT)*Q}, \eqn{T} is the translation
  vector and \eqn{I} is an \eqn{n-by-1} vector with elements \eqn{1}. 
}

\author{C. J. Wong \email{cwon2@fhcrc.org}}

\seealso{\code{\link{gpaSet}}}

\examples{
## Example 1 
x <- matrix(runif(20), nrow=10, ncol=2)+ 1.4
s <- matrix(c(cos(60), -sin(60), sin(60), cos(60)), 
            nrow=2, ncol=2, byrow=TRUE)
xbar <- 2.2 *(x \%*\% s) - 0.1

lt <- iProcrustes(x, xbar, translate=TRUE) ## return linear transformation
lt

## showing result
I <- matrix(1, nrow=nrow(x), ncol=1)
tx <- x - I \%*\% lt$T
## get the transformed matrix xnew
xnew <- lt$scal * (tx \%*\% lt$Q)

if (require(lattice)) {
   xyplot(V1 ~ V2, 
          do.call(make.groups, lapply(list(x=x, xbar=xbar, T.xbar=lt$T.xbar,
                  xnew=xnew),as.data.frame)),  
          group=which, aspect=c(0.7), pch=c(1,3,2,4), col.symbol="black",
	  main=("Align the points in x to xbar"),
          key=list(points=list(pch=c(1,3,2,4), col="black"), space="right",
                   text=list(c("x", "xbar", "T.xbar", "xnew"))))
}

## Example 2. centralized x and xbar prior to using iProcrustes
x <- matrix(runif(10), nrow=5, ncol=2)
s <- matrix(c(cos(60), -sin(60), sin(60), cos(60)), 
            nrow=2, ncol=2, byrow=TRUE)
xbar <- 1.2 *(x \%*\% s) - 2
I <- matrix(1, nrow=nrow(x), ncol=1)
x <- x-(I \%*\% colMeans(x)) ## shift the centroid of points in x to the origin
xbar <- xbar - (I \%*\% colMeans(xbar)) ## shift centroid to the origin
lt <- iProcrustes(x, xbar, translate=FALSE) ## return linear transformation
## only return the rotation/reflection matrix and scalling factor
lt

xnew=lt$scal *(x \%*\% lt$Q) ## transformed matrix aligned to centralized xbar
if (require(lattice)) {
    xyplot(V1 ~ V2,
           do.call(make.groups, lapply(list(x=x,xbar=xbar, 
                   xnew=xnew), as.data.frame)), 
           group=which, auto.key=list(space="right"))
}

}
