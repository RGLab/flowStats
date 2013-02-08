\name{singletGate}
\alias{singletGate}

\title{a polygon gate to select singlet cell populations}


\description{
  
 The function tries to create a singlet polygon gate based area and height of one flow parameter (usually forward scatter). 
	
}

\usage{
singletGate <- function(x, area, height, sidescatter = NULL
						, lower = NULL, upper = NULL
                        ,maxit = 100, nsd = 5, prediction_level = 0.99,
                        filter_id = "singlet") 
}

\arguments{
  
 \item{x}{ A \code{\link[flowCore:flowFrame-class]{flowFrame}}. }
  
  \item{area}{ A character scalar giving the flow parameter that records the signal intensity as peak area. }
  \item{height}{ A character scalar giving the flow parameter that records the signal intensity as peak height. }

  \item{sidescatter}{ A character scalar giving an optional flow prameter name for side scatter signal. }
  
  \item{lower,upper}{}
  
  \item{maxit}{ }
  
  \item{nsd}{}
  \item{prediction_level}{}

  \item{filterId}{ Character, the name assigned to the resulting
    filter. }

    
  
}

\details{

  The algorithm first tries to use the prediction bands from a robust linear model to 
  create a polygon gate to select singlet population.
  
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
\donotrun{
#fr is a flowFrame
sg<-singletGate(fr,area = "FSC-A", height = "FSC-H")
rg
#plot the gate 
xyplot(`FSC-H`~`FSC-A`,fr,filter=sg)
}


}
