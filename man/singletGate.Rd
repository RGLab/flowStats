\name{singletGate}
\alias{singletGate}

\title{Creates a singlet polygon gate using the prediction bands from a robust linear model}


\description{
  
 The function constructs a singlet \code{polygonGate} using the prediction bands
 from a robust linear model trained on the area and height of one flow
 parameter (usually forward scatter). 
	
}

\usage{
singletGate (x, area, height, sidescatter = NULL, maxit = 100,
                        prediction_level = 0.99, filterId = "singlet")
}

\arguments{
  
 \item{x}{ A \code{\link[flowCore:flowFrame-class]{flowFrame}}. }
  
  \item{area}{ A character giving the flow parameter that records the signal intensity as peak area. }
  \item{height}{ A character giving the flow parameter that records the signal intensity as peak height. }

  \item{sidescatter}{ A character giving an optional flow prameter name for the side scatter signal. }
  
  \item{maxit}{ An integer specifying the maximum number of iterations used to train the robust linear model  }
  
  \item{prediction_level}{ A numeric value between 0 and 1 that determines the width of the prediction bands.}

  \item{filterId}{ Character, the name assigned to the resulting filter. }
  
}

\details{

  The singlet \code{polygonGate} is constructed using the prediction
  bands at the minimum and maximum observations in the \code{area}
  channel. The width of the prediction bands is determined by the
  \code{prediction_level} percentile: the larger the value, the wider
  the bands.

  The robust linear model is trained using the \code{\link[MASS]{rlm}}
  function from the \code{MASS} package.
  
}


\value{
  
  \code{\link[flowCore]{polygonGate}}.  

}


\author{ Greg Finak, John Ramey }

\seealso{
  
  \code{\link{rangeGate}},
  \code{\link[flowCore]{polygonGate}}
    
}
  
\examples{
\dontrun{
# fr is a flowFrame
sg <- singletGate(fr, area = "FSC-A", height = "FSC-H")
sg
# plot the gate 
xyplot(`FSC-H` ~ `FSC-A`, fr, filter = sg)
}

}
