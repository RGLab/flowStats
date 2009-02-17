\name{warpSet}
\alias{warpSet}

\title{Normalization based on landmark registration}

\description{

  This function will perform a normalization of flow cytometry
  data based on warping functions computed on high-density region
  landmarks for individual flow channels.
  
}

\usage{
warpSet(x, stains, grouping = NULL, monwrd = TRUE, subsample=NULL, ...)
}


\arguments{
  
  \item{x}{ A \code{\link[flowCore:flowSet-class]{flowSet}}. }
  
  \item{stains}{ A character vector of flow parameters in \code{x} to be
    normalized. }
  
  \item{grouping}{ A character indicating one of the phenotypic
    variables in the \code{phenoData} slot of \code{x} used as a grouping
    factor. The within-group and between-group variance is computed and
    a warning is issued in case the latter is bigger than the former,
    indicating the likely removal of signal by the normalization
    procedure.  }

  \item{monwrd}{ Logical. Compute strictly monotone warping
    functions. This gets directly passed on to
    \code{\link[fda]{landmarkreg}}. }

  \item{subsample}{ Numeric. Reduce the number of events in each \code{flowSet}
    by sub sampling for all density estimation steps and the calculation
    of the warping functions. This can increase computation time for
    large data sets, however it might reduce the accuracy of the density
    estimates. To be used with care. }

  \item{\dots}{ Further arguments that are passed on to
    \code{landmarkreg}. }
  
}

\details{

  Normalization is archived by first identifying high-density regions
  (landmarks) for each \code{\link[flowCore:flowFrame-class]{flowFrame}}
  in the \code{flowSet} for a single channel and subsequently by
  computing warping functions for each \code{flowFrame} that best align
  these landmarks. This is based on the algorithm implemented in the
  \code{landmarkreg} function in the \code{\link[fda:fda-package]{fda}}
  package. An intermediate step classifies the high-density regions, see
  \code{\link{landmarkMatrix}} for details.

  Please note that this normalization is on a channel-by-channel
  basis. Multiple channels are normalized in a loop.  
  
  
}


\value{

  The normalized \code{flowSet}.
  
}

\references{
  
 J.O. Ramsay and B.W. Silverman: Applied Functional Data Analysis,
 Springer 2002

}


\author{ Florian Hahne }

\note{

  We currently use a patched fda version.
  
}

\seealso{

  \code{\link[flowCore:curv1Filter-class]{curv1Filter}},
  \code{\link{landmarkMatrix}}
  
}

\examples{

data(ITN)
dat <- transform(ITN, "CD4"=asinh(CD4), "CD3"=asinh(CD3), "CD8"=asinh(CD8))
lg <- lymphGate(dat, channels=c("CD3", "SSC"),
preselection="CD4",scale=1.5)
dat <- Subset(dat, lg$n2gate)
datr <- warpSet(dat, "CD8", grouping="GroupID", monwrd=TRUE)
if(require(flowViz)){
d1 <- densityplot(~CD8, dat, main="original", filter=curv1Filter("CD8"))
d2 <- densityplot(~CD8, datr, main="normalized", filter=curv1Filter("CD8"))
plot(d1, split=c(1,1,2,1))
plot(d2, split=c(2,1,2,1), newpage=FALSE)
}

}