\name{normalize-methods}
\docType{methods}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize,GatingSetInternal,missing-method}
\alias{normalize,GatingSet,missing-method}
\title{normalize a GatingSet imported with flowWorkspace, using sequential normalization on the manual gates in the GatingHierarchy.}
\description{
The method will step through the gating hierarchy in a breadth first search manner and normalize each dimension and gate not explicitly excluded in skipdims,or skipgates. The normalization approach is based on warpSet, but uses sequential normalization to alternately normalize then perform gating of the cell populations.  This often helps with feature registration of popualtions lower in the gating hierarchy. FSC and SSC, as well as time are generally excluded by default. The rule of thumb, is to only normalize a channel in a gate if it is absolutely warranted.
}
\usage{
normalize(data,x,...)
}
\arguments{

	\item{data}{
	The \code{GatingSet} to be normalized.
}
\item{x}{\code{missing}. Not used in here.}
\item{\dots}{
		Arguments passed to downstream functions.
		
		target:
		The target sample to normalize the other samples in the gating set to. A \code{character} vector. Must match a sample name in \code{x}, otherwise \code{NULL} will use the mean (average) of the peaks identified in all samples

 		populations:
		A \code{character} vector of population names that are to be normalized.  
		

  		dims:
		A \code{character} vector of parameter names to be normalized.  
		
		
		chunksize:
		For a memory-efficient implementation of normalization, set the \code{chunksize},(an \code{integer}), which will peform normalization on chunks of the data of size \code{chunksize}. 
		
		nPeaks:
		A \code{list} of \code{integer} or an \code{integer} vector that specifies the expected number of peaks for each sample. Can be omitted to keep all peaks.

		bwFac:
		The bandwidth for density estimation, a \code{numeric}. Affects the sensitivity for smoothing and detecting distinct peaks.
}
		

    
}
\details{
This function implements sequential normalization using a GatingSet and a set of manual gates. For each gate in the gating hierarchy, the algorithm checks if the gate should be normalized, and which dimensions in the gate should be normalized. If normalization is warranted, this is performed prior to gating. After gating, the counts for the gate in the GatingSet are updated, and the next gate is processed. This is useful in the application of template gates to data that has staining variability in one or more channels.
}
\value{
Returns a \code{GatingSet} of normalized data.
}
\author{
Greg Finak <gfinak@fhcrc.org>
}

\seealso{
See also \code{\link[flowWorkspace]{GatingSet-class}}, \code{\link[flowWorkspace]{GatingHierarchy-class}},  \code{\link[ncdfFlow]{ncdfFlowSet}}
}
\examples{
\dontrun{
#gs is a GatingSet
gs_norm <- normalize(gs
                , target = "M+T panel_903997-25.fcs"
                , populations = "cd27gate"
                , dims = "<Violet A 610/20-A>"
                , minCountThreshold = 100
                , nPeaks = list('cd27gate' = 2)
                , chunksize = 10
                , bwFac = 2
                )



#show the population statistics for before and after normalization				
getPopStats(gs)
getPopStats(gs_norm)		

#plot the gate to see the effects of normalization
grid.arrange(
		plotGate(gs, "cd27gate", type = "densityplot", stack = T)
		,plotGate(gs_norm, "cd27gate", type = "densityplot", stack = T)
		)
}
}
\keyword{methods}




 


