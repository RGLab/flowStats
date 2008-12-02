\name{autoGate}
\alias{autoGate}

\title{ Automated gating of single populations in 2D}
\description{
  This function tries to fit a single norm2Filter based on a rough
  preselection of the data. This function is considered internal. Please
  use the API provided by \code{\link{lymphGate}}.
}

\usage{
autoGate(x, \dots, scale = 2.5)
}

\arguments{
  
  \item{x}{ An object of class \code{\link[flowCore]{flowSet}}}
  
  \item{\dots}{ Named arguments or a list of the ranges used for the
    initial rough preselection. This gets passed on to
    \code{\link[flowCore]{rectangleGate}}, see it's documentation for
    details. }

  \item{scale}{ The \code{scale} parameter that gets passed on to
    \code{\link[flowCore]{norm2Filter}}. }
}

\details{

  The \code{flowSet} is first filtered using a \code{rectangleGate} and
  the \code{norm2Filter} is subsequently fitted to the remaining subset.

}

\value{

  A list with items:
  
  \item{x }{The filtered \code{flowSet}. }
  
  \item{n2gate }{The \code{norm2Filter} object. }

  \item{n2gateResults }{The \code{\link[flowCore]{filterResult}} after
    applying the \code{norm2Filter} on the \code{flowSet}. }
  
}


\author{ Florian Hahne }

\seealso{ \code{\link{lymphGate}}, \code{\link[flowCore]{norm2Filter}} }

\examples{

data(GvHD)
flowStats:::autoGate(GvHD[10:15], "FSC-H"=c(100,500), "SSC-H"=c(0, 400))

}

\keyword{misc}

