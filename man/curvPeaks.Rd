\name{curvPeaks}
\alias{curvPeaks}


\title{ Parse curv1Filter output }


\description{

  Parse the output of
  \code{\link[flowStats:curv1Filter-class]{curv1Filter}} and find modes
  and midpoints of the high-density regions. This function is considered
  to be internal.
  
}

\usage{

curvPeaks(x, dat, borderQuant = 0.01, n = 201, from, to, densities=NULL)

}


\arguments{
  \item{x}{ A
    \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResult}}
    produced by a \code{curv1Filter} operation. }

  \item{dat}{ The corresponding
    \code{\link[flowCore:flowFrame-class]{flowFrame}}. }
  
  \item{borderQuant}{ A numeric in \code{[0,1]} giving the extreme
    quantiles for which high-density regions are ignored. }
  
  \item{n, from, to}{ Arguments are passed on to
    \code{\link[stats]{density}}. }

  \item{densities}{The optional y values of the density estimate
    computed for the respective data.}
  
}


\value{
  
  A list with items
  
  \item{peaks }{x and y locations of the modes of the regions in the
    density estimates. }
  
  \item{regions }{the left and right margins of the regions. }

  \item{midpoints }{the mean of \code{regions}. }
  
  \item{regPoints }{x and y locations of the outline of the significant
    density regions. }

   \item{densFuns }{an approximation function of the density estimate }
  
}

\author{ Florian Hahne}

\seealso{ \code{\link{landmarkMatrix}} }

\examples{

data(GvHD)
tmp <- filter(GvHD[[10]], curv1Filter("FSC-H"))
res <-  flowStats:::curvPeaks(tmp, exprs(GvHD[[10]])[, "FSC-H"])

}