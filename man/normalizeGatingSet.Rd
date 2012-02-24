\name{normalizeGatingSet}
\alias{normalizeGatingSet}
\title{
normalize a GatingSet imported with flowWorkspace, using sequential normalization on the manual gates in the GatingHierarchy.
}
\description{
The function will step through the gating hierarchy in a breadth first search manner and normalize each dimension and gate not explicitly excluded in skipdims,or skipgates. The normalization approach is based on warpSet, but uses sequential normalization to alternately normalize then perform gating of the cell populations.  This often helps with feature registration of popualtions lower in the gating hierarchy. FSC and SSC, as well as time are generally excluded by default. The rule of thumb, is to only normalize a channel in a gate if it is absolutely warranted.
}
\usage{
normalizeGatingSet(x, target = NULL, skipgates = NULL, skipdims = c("FSC-A", "SSC-A", "FSC-H", "SSC-H", "Time"), subsample = NULL, chunksize = 10, nPeaks = list(), bwFac = 2, ...)
}
\arguments{
  \item{x}{
	The \code{GatingSet} to be normalized.
}
  \item{target}{
The target sample to normalize the other samples in the gating set to. A \code{character} vector. Must match a sample name in \code{x}, otherwise \code{NULL} will use the mean (average) of the peaks identified in all samples
}
  \item{skipgates}{
	A \code{numeric} vector of gate indices that are to be skipped (i.e. not normalized). The index calculation is a bit convoluted. We provide a mapping function \code{\link[flowWorkspace]{includedGate2ExcludedGate}} to map from flowJo gate path names to indices.
}
  \item{skipdims}{
	A \code{character} vector of parameter names to be omitted during sequential normalization. For convenience, a mapping function \code{\link[flowWorkspace]{includedChannel2ExcludedChannel}} performs the mapping from explicitly including dimensions for normalization, to a list of excluded dimensions for the function call.
}
  \item{subsample}{
	
}
  \item{chunksize}{i
	For a memory-efficient implementation of normalization, set the \code{chunksize},(an \code{integer}), which will peform normalization on chunks of the data of size \code{chunksize}. 
}
  \item{nPeaks}{
A \code{list} of \code{integer} or an \code{integer} vector that specifies the expected number of peaks for each sample. Can be omitted to keep all peaks.
}
  \item{bwFac}{
The bandwidth for density estimation, a \code{numeric}. Affects the sensitivity for smoothing and detecting distinct peaks.
}
  \item{\dots}{
Additional arguments passed to downstream functions.
}
}
\details{
This function implements sequential normalization using a GatingSet and a set of manual gates. For each gate in the gating hierarchy, the algorithm checks if the gate should be normalized, and which dimensions in the gate should be normalized. If normalization is warranted, this is performed prior to gating. After gating, the counts for the gate in the GatingSet are updated, and the next gate is processed. This is useful in the application of template gates to data that has staining variability in one or more channels. 
}
\value{
Returns a \code{GatingSet} of normalized data.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link[flowWorkspace]{GatingSet-class}}, \code{\link[flowWorkspace]{GatingHierarchy-class}}, \code{\link[flowWorkspace]{includedGate2ExcludedGate}}, \code{\link[flowWorkspace]{includedChannel2ExcludedChannel}}, \code{\link[ncdfFlow]{ncdfFlowSet}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
