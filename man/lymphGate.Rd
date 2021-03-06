\name{lymphFilter-class}
\docType{class}
\alias{lymphFilter-class}
\alias{\%in\%,flowFrame,lymphFilter-method}
\alias{lymphGate}
\alias{lymphFilter}

\title{ Automated gating of elliptical cell populations in 2D. }

\description{

  Cell populations of roughly elliptical shape in two-dimensional projections
  are of huge interest in many flow cytometry applications. This
  function identifies a single such population, potentially from a mixture
  of multiple populations. 
  
}

\section{Extends}{
  
  Class \code{\link[flowCore:parameterFilter-class]{parameterFilter}},
  directly.
  
  Class \code{\link[flowCore:concreteFilter-class]{concreteFilter}},
  by class "parameterFilter", distance 2.

  Class \code{\link[flowCore:filter-class]{filter}}, by class
  "parameterFilter", distance 3.

}

\section{Slots}{

  See \code{Arguments} section for details.
  
  \describe{
    
    \item{\code{preselection}:}{Object of class \code{character}, the
      name of the flow parameter used for preselection. }
    
    \item{\code{rectDef}:}{Object of class \code{list}, the initial
      rectangular selection. }
    
    \item{\code{scale}:}{Object of class \code{numeric}.}

    \item{\code{bwFac}:}{Object of class \code{numeric}.}
    
    \item{\code{parameters}:}{Object of class \code{parameters}, the flow
      parameters to operate on. }
    
    \item{\code{filterId}:}{Object of class \code{"character"}, the
      filter identifier. }
    
  }
}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("lymphFilter",
    parameters, ...)} or using the constructor \code{lymphFilter}. The
  constructor is the recommended way of object instantiation.

}

\usage{
lymphGate(x, channels, preselection=NULL, scale=2.5,  bwFac=1.3,
          filterId="defaultLymphGate", plot=FALSE, \dots)
}

\arguments{
  
  \item{x}{ An object of class \code{\link[flowCore:flowSet-class]{flowSet}}. }

  \item{channels}{ A character vector of length 2 of valid flow
    parameters in \code{x}. }
  
  \item{preselection}{ Either \code{NULL}, in which case this boils down
    to fitting a regular
    \code{\link[flowStats]{norm2Filter}}, a 
    character scalar giving one of the flow parameters in \code{x}, or a
    named list of numerics specifying the initial rough
    preselection. The latter gets passed on to
    \code{\link[flowCore]{rectangleGate}}, see it's documentation for
    details. }

  \item{scale}{ The \code{scaleFactor} parameter that gets passed on to
    \code{norm2Filter}. }

  \item{bwFac}{ The bandwidth factor that gets passed on to
    \code{\link[flowStats:curv1Filter-class]{curv1Filter}}. }

  \item{filterId}{ A character used as filterId. }

  \item{plot}{Logical. Produce plots of filter results}
  
  \item{\dots}{ Additional arguments. }
  
}

\details{

  This algorithm does not apply real mixture modelling, however it is
  able to identify a single elliptical cell population from a mixture of
  multiple such populations. The idea is to first define a rough
  rectangular preselection and, in a second step, fit a bivariate normal
  distribution to this subset only.

  Depending on the value of \code{preselection}, the initial rough
  selection is either

  \describe{

    \item{NULL:}{ No preselection at all}

    \item{character scalar}{ Preselection based on cells that are
      positive for a single marker only. This allows for back-gating,
      for instances by selecting \code{CD4+} T-cells and using this
      information to back-gate lymphocytes in \code{FSC} and
      \code{SSC}. Positive cells are identified using a
      \code{curv1Filter}. }

    \item{a named list of numerics:}{ Preselection by a rectangular
      gate. The items of the list have to be numerics of length one
      giving the gate boundaries in the respective dimensions. }


  }

}

\value{

  An \code{\link[flowCore]{ellipsoidGate}} or list 
  of \code{ellipsoidGate} objects
  
}


\author{ Florian Hahne }

\seealso{ \code{\link[flowStats]{norm2Filter}},
  \code{\link[flowStats:curv1Filter-class]{curv1Filter}} } 

\examples{

library(flowCore)
data(GvHD)
dat <- GvHD[pData(GvHD)$Patient==10]
dat <- transform(dat, "FL4-H"=asinh(`FL4-H`))
lg <- lymphGate(dat, channels=c("FSC-H", "SSC-H"), preselection="FL4-H",scale=1.5)

if(require(flowViz))
xyplot(`SSC-H`~`FSC-H`, dat, filter=lg)


}

\keyword{misc}
\keyword{classes}

