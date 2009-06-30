\name{gpaSet}
\alias{gpaSet}

\title{Multi-dimensional normalization of flow cytometry data}

\description{

  This function performs a multi-dimensional normalization of flow cytometry
  data (\code{flowSets}) using a generalized Procrustes analysis (GPA) method,
  coupled with feature registration based on backgating. 

}

\usage{
gpaSet(x, params, register="backgating", bgChannels=NULL, ref.method="mean", 
       rotation.only=TRUE, iter.max=10, thres=1e-3, plot=FALSE, ...)
}

\arguments{
  \item{x}{A \code{\link[flowCore:flowSet-class]{flowSet}}.}

  \item{params}{A character vector of length 2 describing the channels of
    interest.}

  \item{register}{A character indicating the method to be used for identifying
    features. Avaliable methods include \dQuote{backgating} and
    \dQuote{curve1Filter}. } 

  \item{bgChannels}{A character vector indicating the channels used for
    backgating.}

  \item{ref.method}{A character indicating the method to be used for
    calculating the reference feature. Up to this moment, the only avaliable
    method is \dQuote{mean}.}

  \item{rotation.only}{Logical for coarsing a reflection matrix to a
    rotation matrix.}

  \item{iter.max}{An integer indicating the maximum number of iterations for
    GPA.}

  \item{thres}{A numeric indicating the threshold of the convergence test for
    GPA.}

  \item{fig}{Logical. Indicating whether to plot the resulting
    figures of feature identification and GPA methods. Mainly used for
    testing.}

  \item{...}{Further arguments.}
}

\details{
   Normalization is acheived by first identifying features for each
   \code{\link[flowCore:flowFrame-class]{flowFrame}} in the \code{flowSet} for
   designated channels using backgating, sebsequently leballing features, and
   finally aligning the feautres to a reference feature in the sense of
   minimizing minimizing the Frobenus norm of 
   \deqn{||sFQ - \bar{F}||,}
   where \eqn{s} is a scaler, \eqn{Q} a rotational matrix, \eqn{F} the matrix of
   features, and \eqn{\bar{F}} the reference feature. Both  \eqn{s} and
   \eqn{Q} are solved by using signular value decomposition (SVD).  

}

\value{
  The normalized \code{flowSet} with a set of "GPA" attribute.
}

\references{in progress}

\author{C. J. Wong \email{cwon2@fhcrc.org}}

\examples{
## Example 1:
data(ITN)
ITN <- ITN[1:8, ]
wf <- workFlow(ITN)
tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
add(wf, tl)
x <- Data(wf[["asinh"]])
xy = c("FSC", "SSC")
s <- gpaSet(x, params=xy, bgChannels="CD8", ref.method="mean")
if(require(flowViz)) {
   d1 <- densityplot(~., s, channels=c("FSC", "SSC"), 
                     layout=c(2,1), main="After GPA using bg")
   d2 <- xyplot(FSC ~ SSC, as(s, "flowFrame"), 
                channels=c("FSC", "SSC"), main="All flowFrames")
   plot(d1)
   plot(d2, newpage=FALSE)
}

## view "GPA" attribute
attr(s, "GPA")

## Example 2: using work flow
xy = c("FSC", "SSC")
norm <- normalization(normFun=function(x, parameters, ...)
        gpaSet(x, parameters, ...),
        parameters = c("FSC", "SSC"), 
	arguments=list(bgChannels=bgChannels, register="backgating"),
	normalizationId="Procrustes")

add(wf, norm2, parent="asinh")
s <- Data(wf[["Procrustes"]])
if(require(flowViz)) {
   d1 <- densityplot(~., s, channels=c("FSC", "SSC"), 
                     layout=c(2,1), main="After GPA using bg")
   d2 <- xyplot(FSC ~ SSC, as(s, "flowFrame"), 
                channels=c("FSC", "SSC"), main="All flowFrames")
   plot(d1)
   plot(d2)
}
}
