\name{normQA}
\alias{normQA}

\title{Normalization quality assessment}

\description{Create QA plots for a flow cytometry normalization process.}

\usage{
normQA(data, morph = c("^fsc", "^ssc"), channels, odat = NULL, ask = names(dev.cur()) != "pdf", grouping = NULL, tag.outliers = FALSE, peaksOnly = TRUE)
}

\arguments{
  
  \item{data}{a normalized \code{\link[flowCore:flowSet-class]{flowSet}}.}
  
  \item{morph}{A character vector of channel names to use for the
    backgating into the morphological channels.}
  
  \item{channels}{The channels for which to create plots. Defaults to
    all normalized channels.}
  
  \item{odat}{The original data set, always needed if there are no
    warping functions available.}
  
  \item{ask}{Ask before creating a new plot.}
  
  \item{grouping}{A grouping variable in \code{data}'s \code{phenoData} slot.}
  
  \item{tag.outliers}{Logical. Add sample name to outliers in the
    plots.}
  
  \item{peaksOnly}{Logical. Only use data when a peak was detected in a
    particular sample. If set to \code{FALSE}, a average peak location
    is estimated.}
  
}

\details{

  This function assumes that the necessary information has been added as
  attributes to \code{data} during the normalization procedure. Depending
  on the available information, a set of QA plots is generated.
  Available plots are:

    Amount of peak adjustment
    
	Warping functions

    Landmark classification confidence
    
Backgating of peak events in morphological channels

}
\value{
This function is called for its side effect of generating plots.
}

\author{
Florian Hahne
}
