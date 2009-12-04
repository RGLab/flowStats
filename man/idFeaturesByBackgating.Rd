\name{idFeaturesByBackgating}
\alias{idFeatures}
\alias{idFeaturesByBackgating}

\title{(Internal use only) Identify features of flow cytometry data using
  backgating} 

\description{Identify and labeling significant features using divisive
  clustering method such as \code{\link[cluster:diana]{diana}}.
}

\usage{
idFeaturesByBackgating(bg, nDim, thres.sigma=2.5, lambda=1/1.4,
                       reference.method="median",
                       plot.workflow=FALSE, ask=names(dev.cur())!="pdf")
}

\arguments{
  \item{bg}{A data frame containing subpopulations on channels of
    interests. Must be a returning result from \code{flowStats:::backGating}}
  \item{nDim}{An integer indicating the length of channels of interest.}
  \item{thres.sigma}{An numerical value indicating the threshold at
    which to cut tree, e.g., as resulting from 'diana', into several clusters.}
  \item{lambda}{a numerical value of the threshold for deciding outlier
    clusters. The default value, 1/1.4, is derived using "robust
  statistics". See 'details'.}
  \item{reference.method}{A character vector indicating the method for
    computing the reference features. If \code{median}, the reference
    feature is defined by the medain of clusters of features. Valid
    methods include \code{median} and \code{mean}. }
  \item{plot.workflow}{Logical. If TURE, display the workflow of feature
    identification.}
  \item{ask}{Logical. If TRUE, the display operates in interactive mode.}
}

\details{
  Using the resulting data frame from \code{backGating} as potential
  features, the algorithm follows four major steps: (i) centering 
  the potential features, which yields the returning value
  \code{TransMatrix}, (ii) using \code{diana} to compute a clustering of
  the potential features, (iii) cutting the tree into several clusters,
  and (iv) accessing outliers and rendering the final registered
  features with labels.

  In step three, the threshold for cutting the tree is computed by
  \deqn{sd * thres.sigma,}
  where \eqn{sd} is the standard deviation of the distribution of the
  height between entities computed by \code{diana}.

  A cluster is determined as an outlier if the number of its members is
  less than the median of the numbers of all clusters' members times 'lambda'.
}

\value{ 
  \item{register}{A list containing registered features for each sample.}
}
  
\author{Chao-Jen Wong}
\seealso{
  \code{\link[cluster:diana]{diana}},
  \code{\link{BackGating}}, \code{\link{gpaSet}}.
}

\examples{
\dontrun{
data(ITN)
wf <- workFlow(ITN)
tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
dat <- trnasformList(ITN, tl)
bg <- backGating(dat, xy=c("FSC", "SSC"), channels="CD3")
}

data(BackGating)
results <- flowStats:::idFeaturesByBackgating(bg=BackGating, nDim=2,
                       plot.workflow=TRUE, ask=TRUE)

}
  
\keyword{misc}
