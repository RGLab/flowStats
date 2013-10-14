\name{BackGating}
\alias{BackGating}
\docType{data}
\title{Sample backgating results}
\description{
  A data frame containing the sub-populations of ITN dataset corresponding
  to the high-density areas on "FSC" and "SSC" channels. This dataset is
  yielded by \code{backGating} on channel CD3,   CD8, and CD4 of the
  \code{ITN} sample data. 
}

\usage{data(BackGating)}

\source{ Results from executing the following code:
  
  data(ITN)
  
  flowStats:::backGating(ITN, xy=c("FSC", "SSC"),
      channels=c("CD3", "CD8", "CD4"))
}

\keyword{datasets}
