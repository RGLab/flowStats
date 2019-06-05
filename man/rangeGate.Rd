\name{rangeGate}

\alias{rangeGate}
\alias{rangeFilter}
\alias{rangeFilter-class}
\alias{\%in\%,flowFrame,rangeFilter-method}
\alias{oneDGate}

\title{ Find most likely separation between positive and negative
  populations in 1D }


\description{
  
 The function tries to find a reasonable split point between the two
 hypothetical cell populations "positive" and "negative".
}

\usage{
rangeGate(x, stain, alpha="min", sd=2, plot=FALSE, borderQuant=0.1,
absolute=TRUE, filterId="defaultRectangleGate", positive=TRUE,
refLine=NULL, simple = FALSE,...)

rangeFilter(stain, alpha="min", sd=2, borderQuant=0.1, 
filterId="defaultRangeFilter")
}

\arguments{
  
 \item{x}{ A \code{\link[flowCore:flowSet-class]{flowSet}} or
    \code{\link[flowCore:flowFrame-class]{flowFrame}}. }
  
  \item{stain}{ A character scalar giving the flow parameter for which
    to compute the separation. }
  
  \item{alpha}{ A tuning parameter that controls the location of the
    split point between the two populations. This has to be a numeric in
    the range \code{[0,1]}, where values closer to 0 will shift the
    split point closer to the negative population and values closer to 1
    will shift towards the positive population. Additionally, the value
    of \code{alpha} can be \code{"min"}, in which case the split point
    will be selected as the area of lowest local density between the two
    populations. }
  
  \item{sd}{For the case where there is only a single population, the
    algorithm falls back to esitmating the mode of this population and a
    robust measure of the variance of it distribution. The \code{sd}
    tuning parameter controls how far away from the mode the split point
    is set. }
  
  \item{plot}{ Create a plot of the results of the computation. }
  
  \item{borderQuant}{Usualy the instrument is set up in a way that the
    positive population is somewhere on the high end of the measurement
    range and the negative population is on the low end. This parameter
    allows to disregard populations with mean values in the extreme
    quantiles of the data range. It's value should be in the range
    \code{[0,1]}. }
  
  \item{absolute}{Logical controling whether to classify a population
    (positive or negative) relative to the theoretical measurment range
    of the instrument or the actual range of the data. This can be set
    to \code{TRUE} if the alignment of the measurment range is not
    optimal and the bulk of the data is on one end of the theoretical
    range.}

  \item{filterId}{ Character, the name assigned to the resulting
    filter. }

  \item{positive}{ Create a range gate that includes the positive
    (\code{TRUE}) or the negative (\code{FALSE}) population. }

  \item{refLine}{Either \code{NULL} or a numeric of lenth 1. If
    \code{NULL}, this parameter is ignored. When it is set to a numeric, the
    minor sub-population (if any) below this reference line 
    will be igored while determining the separator between positive and
    negative. }  
  \item{simple}{ \code{logical} scalar indicating whether to use a simple peak finding version of density1d algorithm.
   }
  \item{\dots}{ Further arguments. }
}

\details{

  The algorithm first tries to identify high density regions in the
  data. If the input is a \code{flowSet}, density regions will be
  computed on the collapsed data, hence it should have been normalized
  before (see \code{\link{warpSet}} for one possible normalization
  technique). The high density regions are then clasified as positive
  and negative populations, based on their mean value in the theoretical
  (or absolute if argument \code{absolute=TRUE}) measurement range. In
  case there are only two high-density regions the lower one is usually
  clasified as the negative populations, however the heuristics in the
  algorithm will force the classification towards a positive population
  if the mean value is already very high. The \code{absolute} and
  \code{borderQuant} arguments can be used to control this
  behaviour. The split point between populations will be drawn at the
  value of mimimum local density between the two populations, or, if the
  \code{alpha} argument is used, somewhere between the two populations
  where the value of alpha forces the point to be closer to the negative
  (\code{0 - 0.5}) or closer to the positive population (\code{0.5 -
  1}).

  If there is only a single high-density region, the algorithm will fall
  back to estimating the mode of the distribution
  (\code{\link[MASS]{hubers}}) and a robust measure of it's variance
  and, in combination with the \code{sd} argument, set the split point
  somewhere in the right or left tail, depending on the classification
  of the region.

  For more than two populations, the algorithm will still classify each
  population into positive and negative and compute the split point
  between those clusteres, similar to the two population case.
  
  
  The \code{rangeFilter} class and constructor provide the means to
  treat \code{rangeGate} as regular \code{flowCore} filters.
  
}

\section{Methods}{
  
  \describe{
  
     \item{\%in\%}{\code{signature(x =
	"flowFrame", table = "rangeFilter")}: the work horse for doing the
      actual filtering. Internally, this simply calls the \code{rangeGate}
      function. }
  }
}

\value{

  A range gate, more explicitely an object of class
  \code{\link[flowCore]{rectangleGate}}.  

}


\author{ Florian Hahne, Kyongryun Lee }

\seealso{
  
  \code{\link{warpSet}},  \code{\link{rangeGate}},
  \code{\link[flowCore]{rectangleGate}}
    
}
  
\examples{
library(flowCore)
data(GvHD)
dat <- GvHD[pData(GvHD)$Patient==10]
dat <- transform(dat, "FL4-H"=asinh(`FL4-H`), "FL3-H"=asinh(`FL3-H`))
rg <- rangeGate(dat, "FL4-H", plot=TRUE)
rg
split(dat, rg)

## Test rangeGate when settting refLine=0; it does not do anything since
## there is no sub-population below zero.
rangeGate(dat, "FL4-H", plot=FALSE, refLine=0)

rf <- rangeFilter("FL4-H")
filter(dat, rf)
}
