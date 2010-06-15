\name{density1d}

\alias{density1d}

\title{ Find most likely separation between positive and negative
  populations in 1D }


\description{
  
 The function tries to find a reasonable split point between the two
 hypothetical cell populations "positive" and "negative". This function
 is considered internal, please use the API provided by
 \code{\link{rangeGate}}.

}

\usage{
density1d(x, stain, alpha = "min", sd = 2, plot = FALSE, borderQuant =
0.1, absolute = TRUE, inBetween = FALSE, refLine=NULL, ...)
}

\arguments{
  
 \item{x}{ A \code{\link[flowCore:flowSet-class]{flowSet}} or
    \code{\link[flowCore:flowSet-class]{flowFrame}}. }
  
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

  \item{inBetween}{Force the algorithm to put the separator in between
    two peaks. If there are more than two peaks, this argument is
    ignored.}

  \item{refLine}{Either \code{NULL} or a numeric of lenth 1. If
    \code{NULL}, this parameter is ignored. When it is set to a numeric, the
    minor sub-population (if any) below this reference line 
    will be igored while determining the separator between positive and
    negative.}
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
  
}


\value{


  A numeric indicating the split point between positive and negative populations.
}


\author{ Florian Hahne }

\seealso{
  
  \code{\link{warpSet}},  \code{\link{rangeGate}}
    
}
  
\examples{

data(GvHD)
dat <- GvHD[pData(GvHD)$Patient==10]
dat <- transform(dat, "FL4-H"=asinh(`FL4-H`), "FL3-H"=asinh(`FL3-H`))
d <- flowStats:::density1d(dat, "FL4-H", plot=TRUE)
if(require(flowViz))
densityplot(~`FL4-H`, dat, refline=d)

## tweaking the location
flowStats:::density1d(dat, "FL4-H", plot=TRUE, alpha=0.8)

## only a single population
flowStats:::density1d(dat, "FL3-H", plot=TRUE)
flowStats:::density1d(dat, "FL3-H", plot=TRUE, sd=2)

}
