\name{gaussNorm}
\alias{gaussNorm} 

\title{Per-channel normalization based on landmark registration}

\description{
  
  This funciton normalizes a set of flow cytometry data samples by
  identifying and aligning the high density regions (landmarks or peaks) for each channel.
  The data of each channel is shifted in such a way that the identified high
  density regions are moved to fixed locations called base landmarks.  
}

\usage{
  gaussNorm (flowset, channel.names, max.lms=2, base.lms=NULL, debug=FALSE, fname='')

}


\arguments{
  
  \item{flowset}{ A \code{\link[flowCore:flowSet-class]{flowSet}}. }
  
  \item{channel.names}{ A character vector of flow parameters in \code{flowset} to be
    normalized. }
  
  \item{max.lms}{ A numeric vector of the maximum number of base
  landmarks to be used for normalizing each channel. If it has only one
  value that will be used as the maximum number of base landmarks for
  all the channels.  }

  \item{base.lms}{ A list of vector for each channel that contains the
  base landmarks for normalizing that channel. If not specified the base
  landmarks are computed from the set of extracted landmarks.}

  \item{debug}{ Logical. Forces the function to draw before and
    after normalization plots for each sample. The plot of the i-th sample is stored in \code{paste(fname, i)} file.}
  
  \item{fname}{ The pre- and post- normalization plots of the i-th sample  is
  stored in \code{paste(fname, i)} file if debug is set to
  \code{TRUE}. If default value is used the plots are drawn on separate X11 windows
  for each sample. In this case, the function waits for a user input to draw the plots
  for the next sample.}

}

\details{

  Normalization is archived in three phases: (i) identifying high-density regions
  (landmarks) for each \code{\link[flowCore:flowFrame-class]{flowFrame}}
  in the \code{flowSet} for a single channel; (ii) computing the best
  matching between the landmarks and a set of fixed reference
  landmarks for each channel called base landmarks; (iii) manipulating
  the data of each channel  in such a way that each landmark is moved to
  its matching base landmark. 
  Please note that this normalization is on a channel-by-channel
  basis. Multiple channels are normalized in a loop.  
  
  
}


\value{

  The normalized \code{flowSet}.
  
}


\author{ Alireza Hadj Khodabakhshi }

\examples{

data(ITN)
dat <- transform(ITN, "CD4"=asinh(CD4), "CD3"=asinh(CD3), "CD8"=asinh(CD8))
lg <- lymphGate(dat, channels=c("CD3", "SSC"),
preselection="CD4",scale=1.5)
dat <- Subset(dat, lg$n2gate)
datr <- gaussNorm(dat, "CD8")
if(require(flowViz)){
d1 <- densityplot(~CD8, dat, main="original", filter=curv1Filter("CD8"))
d2 <- densityplot(~CD8, datr, main="normalized", filter=curv1Filter("CD8"))
plot(d1, split=c(1,1,2,1))
plot(d2, split=c(2,1,2,1), newpage=FALSE)
}

}