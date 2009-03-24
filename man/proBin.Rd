\name{proBin}
\alias{proBin}

\title{ Probability binning - a metric for evaluating multivariate differences}

\description{
  
  This function divides the flowframe events into bins such that each
  bin contains the same number of events. The number of events falling
  into each bin can then be compared across the control and test samples
  using statistical methods such as the Chi-squared test.

}

\usage{
proBin(m, minEvents)
}

\arguments{
  
  \item{m}{ An object of class \code{\link[flowCore]{flowFrame}}}
  
  \item{minEvents}{ The \code{minEvents} The minimum number of events in
    each bin. (i.e. the termination criterion for the probability
    binning algorithm)}
  
}

\details{

  The \code{flowSet} is first filtered using a \code{rectangleGate} and
  the \code{norm2Filter} is subsequently fitted to the remaining subset.

}

\value{

  A list with items:
  
  \item{table }{A \code{data.frame} that stores information regarding
    each node of the tree generated during each stage of the probability
    binning algorithm.  Each row in the table represents a node, the
    first row representing the original \code{flowFrame} matrix.

    The dataIndx column provides indexes for retrieving the matrices
    during each stage of the binning process from the enviroment
    \code{data} .

    The parent field indicates the row number in the table that holds
    the parent information for the corresponding node.

    The left and right columns indicates the row numbers in the table
    that stores information regarding the children of that particular
    node. The leaf nodes that hold the binned data can be identified by
    the nodes with the left of right values of zero( ie. no children
    nodes)

    The visited column is used internally by the algorithm to check if a
    particular node has been visited during the computation process.
    
  }
  
  \item{data }{An enviroment in which the matrices generated during each
    stage of the probability binning process is stored. The matrices
    stored at the leaf nodes represent the binned events obtained after
    the stop criterion of minEvents has been achieved. These can be
    identified by the corresponding dataIndx fields provided by the rows
    in the table with the left or right column values of zero.  }

  \item{limits }{A list containing the the boundaries of each hyperplane
    generated during probability binning}

  \item{splitPars }{A \code{data.frame} containing two columns splitCol
    - indicates the column number of the \code{flowFrame} , the split
    was performed.

    splitMed - The median value which was used as the threshold for
    splitting the \code{flowFrame}

    The splitCol and splitMed parameters are utilized by the plotBins
    and shadeBins functions in visualizing the differences between
    control and test sample cases.
    
  }
  
}


\author{ Nishant Gopalakrishnan }

\seealso{ \code{\link{plotBins}}, \code{\link{binByRef}}}

\examples{

data(GvHD)
res<-proBin(GvHD[[1]][,1:7],200)  # Time column in the flowFrame is removed

}

\keyword{misc}