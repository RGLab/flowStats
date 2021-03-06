\name{quadrantGate}
\alias{quadrantGate}

\title{Automated quad gating}

\description{
 This function tries to find the most likely separation of
 two-dimensional flow cytometry in four quadrants.  
}

\usage{
quadrantGate(x, stains, alpha=c("min", "min"), sd=c(2, 2), plot=FALSE,
    filterId="defaultQuadGate", refLine.1=NULL, refLine.2=NULL
				 		,rare=c(FALSE,FALSE)
				 		,sig=c(NULL,NULL)
						,...)
}

\arguments{
  
  \item{x}{ A \code{\link[flowCore:flowSet-class]{flowSet}} or
    \code{\link[flowCore:flowSet-class]{flowFrame}}. }
  
  \item{stains}{ A character vector of length two giving the two flow
    parameters for which the quad gate is to be computed. }
  
  \item{alpha, sd}{ Tuning factors to control the computation of the
    gate boundaries. See \code{\link{rangeGate}} for details. }
  
  \item{plot}{ Logical. Produce plots of intermediate results. }

  \item{filterId}{ Character, the name assigned to the resulting
    filter. }

  \item{refLine.1}{Either \code{NULL} or a numeric of lenth 1. If
    \code{NULL}, this parameter is ignored. When it is set to a numeric, the
    minor sub-population (if any) below this reference line in the first
    stain channel will be igored while determining the separator between
    positive and negative.}

  \item{refLine.2}{Either \code{NULL} or a numeric of lenth 1. If
    \code{NULL}, this parameter is ignored. When it is set to a numeric, the
    minor sub-population (if any) below this reference line in the second
    stain channel will be igored while determining the separator between
    positive and negative.}
	
	\item{rare}{logical flags for two channels, Refer to \code{\link{density1d}} for more details.
	}
	
	\item{sig}{parameters for two channels. Refer to \code{\link{density1d}} for more details.
	}
  \item{\dots}{ Additional arguments }
}

\details{

  The most likely separation between postitive and negative stains for
  two-dimensional data is computed based on density
  estimates. Essentially, the gate parameters are first fitted
  separately for the two parameters and later combined. See the
  documentation for \code{\link{rangeGate}} for details. There is a
  certain amount of heuristics involved in this process. The algorithm
  can be slightly tweaked using the \code{alpha} and \code{sd}
  arguments. Their values will be recycled for the two dimensions unless
  explicitely given as vectors of length 2.
  
}

\value{
  
 An object of class \code{\link[flowCore]{quadGate}}. 

} 

\author{Florian Hahne } 

\seealso{

  \code{\link[flowCore]{quadGate}},
  \code{\link{rangeGate}}

}

\examples{
\dontrun{
library(flowCore)
data(GvHD)
dat <- GvHD[pData(GvHD)$Patient==10]
dat <- transform(dat, "FL4-H"=asinh(`FL4-H`), "FL2-H"=asinh(`FL2-H`))
qg <- quadrantGate(dat, c("FL2-H", "FL4-H"))
qg

if(require(flowViz))
xyplot(`FL2-H`~`FL4-H`, dat, filter=qg)

qg <- quadrantGate(dat, c("FL2-H", "FL4-H"), alpha=c(0.1, 0.9), plot=TRUE)
qg
split(dat, qg)
}
}
