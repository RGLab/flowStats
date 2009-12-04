\name{gpaSet}
\alias{gpaSet}

\title{Multi-dimensional normalization of flow cytometry data}

\description{

  This function performs a multi-dimensional normalization of flow cytometry
  data (\code{flowSets}) using a generalized Procrustes analysis (GPA) method.
   
}

\usage{
gpaSet(x, params, register="backgating", bgChannels=NULL,
       bg=NULL,  rotation.only=TRUE,
       downweight.missingFeatures=FALSE, thres.sigma=2.5,
       show.workflow=FALSE,
       ask=names(dev.cur())!="pdf")
}

\arguments{
  \item{x}{A \code{\link[flowCore:flowSet-class]{flowSet}}.}

  \item{params}{A character vector of length 2 describing the channels of
    interest.}

  \item{register}{A character indicating the method to be used for identifying
    features. Available methods include \dQuote{backgating} and
    \dQuote{curve1Filter}. } 
 
  \item{bgChannels}{A character vector indicating the channels used for
    backgating. If NULL, \code{backGating} will find the appropriate
    backgating channels.}
  
  \item{bg}{A data frame as the returning value of the
    \code{backGating} function. If not NULL, \code{gpaSet} will skip the
    \code{backGating} process and use the given data frame to extract
    potential features.}

  \item{rotation.only}{Logical for coarsing a reflection matrix to a
    rotation matrix.}

  \item{downweight.missingFeatures}{Logical. If TRUE, the missing
  features, labeled as bogus features, are downweighted to zero. See details. }

  \item{thres.sigma}{A numerical value indicating the threshold of where
  to cut the tree, e.g., as resulting from \code{diana}, into several clusters.
  It is default to 2.5 sigma of the distribution of the
  heights of the cluster points.}

  \item{show.workflow}{Logical. If TRUE, the workflow of 
    \code{gpaSet} will be displayed.}

  \item{ask}{Logical. If TRUE, the display operates in interactive
    mode.}
    
}

\details{
   Normalization is achieved by first identifying features for each
   \code{\link[flowCore:flowFrame-class]{flowFrame}} in the \code{flowSet} for
   designated channels using backgating, subsequently labeling features, and
   finally aligning the features to a reference feature in the sense of
   minimizing  the Frobenus norm of 
   \deqn{||sFQ - \bar{F}||,}
   where \eqn{s} is a scalar, \eqn{Q} a rotational matrix, \eqn{F} the matrix of
   features, and \eqn{\bar{F}} the reference feature. Both  \eqn{s} and
   \eqn{Q} are solved by using singular value decomposition (SVD).

   Note that if feature \eqn{F_{ij}} is missing, it is given a bogus
   value as \eqn{\bar{F}_{ij}}.
   
   If \code{downweight.missingFeatures} is TRUE, the cost function
   becomes
   \deqn{||s W_0 FQ - W_0 \bar{F}||,}
   where the weighting function \eqn{W_0} is zero if the corresponding
   feature is bogus.
   
}

\value{
  The normalized \code{flowSet} with "GPA" attribute.
}

\references{in progress}

\author{C. J. Wong \email{cwon2@fhcrc.org}}

\examples{
## Example 1: calling up gpaSet directly
data(ITN)
data(BackGating)

tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
dat <- transform(ITN, tl)

xy = c("FSC", "SSC")
bgChannels = c("CD8", "CD4", "CD3")
## bg <- flowStats:::backGating(x, xy=xy, channels=bgChannels)
s <- gpaSet(dat, params=xy, bgChannels=bgChannels, bg=BackGating)

if(require(flowViz)) {
   d1 <- densityplot(~., s, channels=c("FSC", "SSC"), 
                     layout=c(2,1), main="After GPA using bg")
   d2 <- xyplot(FSC ~ SSC, as(s, "flowFrame"), 
                channels=c("FSC", "SSC"), main="All flowFrames")
   plot(d1)
   plot(d2)
}

## view "GPA" attribute
attr(s, "GPA")

\dontrun{
## Example 2: using work flow and normalization objects
data(ITN)
ITN <- ITN[1:8, ]
wf <- workFlow(ITN)
tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
add(wf, tl)
x <- Data(wf[["asinh"]])
## normalize 'FSC' and 'SSC' channels
norm <- normalization(normFun=function(x, parameters, ...)
        gpaSet(x, parameters, ...),
        parameters = c("FSC", "SSC"), 
	arguments=list(bgChannels=c("CD8", "CD3"),
                       register="backgating"),
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
} ## end of dontrun

}
