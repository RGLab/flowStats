## Align data in a flowSet by estimating high density regions and using this
## information as landmarks. This works separately on each parameter.


warpSetGS <- function(x,stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
		peakNr=NULL, clipRange=0.01, nbreaks=11, fres, bwFac=2,
		warpFuns=FALSE,gate=NULL,target=NULL,chunksize=10,
		...){
	#Perform normalization on the subset
	#return the normalized data
	if(!inherits(x,"GatingSet"))
		stop("x must be of class GatingSet")
	if(!flowWorkspace:::isNcdf(x)){
		#TODO code the regular flowSet (not ncdfFlowSet) normalization code.
		message("Gating Set not gated using netcdf. We'll use the regular warpSet function");
		if(is.null(gate)){
			flowset<-getData(x);
		}else{
			if(gate<=length(flowWorkspace:::getNodes(x[[1]], showHidden = TRUE))){
				flowset<-getData(x,gate);
			}else{
				stop("gate ",gate," out of range");
			}
		}
		colnames(flowset)<-colnames(getData(x[[1]]))
		warped<-warpSet(x=flowset,stains=stains,grouping=grouping,monwrd=monwrd,subsample=subsample,peakNr=peakNr,clipRange=clipRange,nbreaks=nbreaks,fres=fres,bwFac=bwFac,warpFuns=warpFuns,target=target,...);
	}else{
		if(is.null(gate)){
			#ncflowset<-graph:::nodeData(x[[1]]@tree,x[[1]]@nodes[1],"data")[[1]][["data"]]$ncfs
			ncflowset<-flowWorkspace:::flowData(x)
			#subset for the correct samples
			ncflowset<-ncflowset[flowWorkspace:::sampleNames(x)]
			#ncflowset<-ncflowset[setdiff(sampleNames(ncflowset),setdiff(sampleNames(ncflowset),sampleNames(x)))]
			
		}else{
			#check if the gate is in range
			if(gate<=length(flowWorkspace:::getNodes(x[[1]], showHidden = TRUE))){
#				browser()
              ncflowset <- getData(x,gate)
				#subset for the correct samples
				#ncflowset<-ncflowset[setdiff(sampleNames(ncflowset),setdiff(sampleNames(ncflowset),sampleNames(x)))]
				ncflowset <- ncflowset[sampleNames(x)]
			}else{
				stop("gate ",gate," out of range");
			}
		}
		#Update the column names on the ncflowset so that they are consistent with the column names on the GatingSet
		colnames(ncflowset)<-colnames(getData(x[[1]]))
		warped<-warpSetNCDF(x=ncflowset,stains=stains,grouping=grouping,monwrd=monwrd,subsample=subsample,peakNr=peakNr,clipRange=clipRange,nbreaks=nbreaks,fres=fres,bwFac=bwFac,warpFuns=warpFuns,target=target,chunksize=chunksize,...);
		#Should the warped ncdfflowset be assigned to the GatingSet data environment?
		#Do we now need to clear the indices after calling Subset?
		return(warped);
	}
}
warpSetNCDFLowMem <- function(x, stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
		peakNr=NULL, clipRange=0.01, nbreaks=11, fres=NULL, bwFac=2,
		warpFuns=FALSE,target=NULL,chunksize=10,isNew=FALSE,newNcFile=NULL,
		...)
{
    .Deprecated("warpSetNCDF")
	## Some type-checking first
	flowCore:::checkClass(x, "flowSet")
	flowCore:::checkClass(stains, "character")
	mt <- stains %in% colnames(x)
	if(!all(mt))
		stop("Invalid stain(s) not matching the flowSet:\n    ",
				paste(stains[!mt], collapse=", "))
	#expData <- as(x, "list")
	#expData should now be x...
#	browser()
	if(isNew){
		expData<-ncdfFlow:::clone.ncdfFlowSet(x,isNew=TRUE,isEmpty=FALSE,ncdfFile=newNcFile)
	}else{
		expData<-x;
	}
	samples<-sampleNames(expData)
	ranges<-vector("list",length(x));
	names(ranges)<-sampleNames(x)
	for(i in 1:length(x)){ranges[[i]]<-range(x[[i]])}
	#ranges <- fsApply(x, range) 
	if(!is.null(grouping))
		flowCore:::checkClass(grouping, "character", 1)
	if(!is.null(subsample))
	{
		flowCore:::checkClass(subsample, "numeric", 1)
		## subsample data set before all density estimation steps
		#Does Subset clobber anything in a permanent way for an ncdfFlowSet? Check this, otherwise use subsample.
		x <- Subset(x, sampleFilter(size=subsample))
	}
	flowCore:::checkClass(monwrd, "logical", 1)
	flowCore:::checkClass(bwFac, "numeric", 1)
	
	##Chunking for filterResult
	chunkgroups<-(1:length(expData))%/%chunksize
	chunksamples<-split(samples,chunkgroups)
	chunks.size<-table(chunkgroups)
	nchunks<-length(chunks.size)
	densY<-vector("list",nchunks)
	peaksandregions<-vector("list",nchunks)
	#Goal here is to compute filterResults in chunks, then to compute and register landmarks, which take up MUCH less space.
	#chunk the data and compute filterResults
	## find landmarks
	for(k in 1:nchunks){
		if(is.null(formals()$fres))
		{
			fres<-vector("list",nchunks)
			#fres <- list()
			for(p in stains){
				cat("\rEstimating landmarks for channel", p, " Sample chunk ",k,"...")
				#fres[[p]] <- filter(x, curv1Filter(p, bwFac=bwFac))
				fres[[k]][[p]]<-filter(x[chunksamples[[k]]],curv1Filter(p,bwFac))
			}
			cat("\n")    
		}
		else
		{
			if(!is.list(fres[[k]]) || !all(stains %in% names(fres[[k]])))
				stop("The supplied list of filter results does not match ",
						"the channels to warp.")
			if(!all(sapply(fres, is, "filterResult")))
				stop("The argument 'fres' has to be a list of filterResultLists.")
			cat("Extracting landmarks from user supplied filterResults\n")
		}
		
		## define some variables
		nb <- 1001
		lm <- list()
		z <- NULL 
		
		## iterate over stains
		eps <- .Machine$double.eps
		for(p in stains)
		{
			## set up fda parameters
			extend <- 0.15
			from <- min(sapply(ranges, function(z) z[1,p]-diff(z[,p])*extend), na.rm=TRUE)
			to <- max(sapply(ranges, function(z) z[2,p]+diff(z[,p])*extend), na.rm=TRUE)
			wbasis <- create.bspline.basis(rangeval=c(from, to),
					norder=4, breaks=seq(from, to, len=nbreaks))
			
			WfdPar <- fdPar(wbasis, 1, 1e-4)
			##append to the previously computed density estimates of the last chunk and now all stains
			densY[[k]][[p]] <- t(fsApply(x[chunksamples[[k]]], function(z){
								r <- range(z)[,p]
								z <- exprs(z)
								z <- z[z[,p]>r[1]+eps & z[,p]<r[2]-eps, p]
								density(z, from=from, to=to, n=nb, na.rm=TRUE)$y
							}))
			argvals <- seq(from, to, len=nb) 
			##Can't call landmarkMatrix until we have ALL the filterResults.. Need to go into the landmarkMatrix function and factor our the code that deals with filterResults so that we pull out the relevant information chunk by chunk.
			peaksandregions[[k]][[p]]<-flowStats:::getPeakRegions(x[chunksamples[[k]]],fres[[k]],p,border=clipRange,peakNr=peakNr,densities=densY[[k]][[p]],n=nb)
			
			## create matrix of landmarks from curv1Filter peaks.. append to the chunked landmarks.. and all stains.
			#cat("Registering curves for parameter", p, "...\n")
			#landmarks[[k]][[p]] <- flowStats:::landmarkMatrix(x[chunksamples[[k]]], fres[[k]], p, border=clipRange, peakNr=peakNr,densities=densY[[k]][[p]], n=nb)
			fres[[k]][[p]]<-list();gc();
		}
	}
	## create matrix of landmarks from peaks.. append to the chunked landmarks.. and all stains.
	
	dns<-list()
	landmarks<-list()
	for(p in stains){
		dns[[p]]<-do.call(cbind,lapply(densY,function(q)q[[p]]))
		cat("Registering curves for parameter", p, "...\n")
		landmarks[[p]] <- flowStats:::landmarkMatrixWithoutFilterResult(x, peaksandregions, p, border=clipRange, peakNr=peakNr)
	}
	#TODO check all landmark matrices to ensure they are valid.
	densY<-dns;
	rm(dns);gc()
	
	for(p in stains){
		#TODO this should be computed outside the loop once we have all the density estimates.
		fdobj   <- Data2fd(argvals,densY[[p]], wbasis)
#                           argnames = c("x", "samples", "density"))
		densY[[p]]<-NULL;gc();
		## check if we remove signal between groups
		sig <- 0.05
		if(!is.null(grouping)){
			if(!grouping %in% names(pData(x)))
				stop("'", grouping, "' is not a phenoData variable.")
			grps <- as.factor(pData(x)[,grouping])
			anv <- numeric(ncol(landmarks[[p]]))
			for(i in seq_len(ncol(landmarks[[p]])))
				anv[i] <- anova(lm(landmarks[[p]][,i] ~ grps))$Pr[1]
			if(any(anv < sig))
				message("Within-group variances are smaller than ",
						"across-group variances for stain ", p,
						".\nWarping might have removed signal!")
		}
		## fill NAs with column medians
		regions <- attr(landmarks[[p]], "regions")
		dists <- attr(landmarks[[p]], "cdists") 
		attr(landmarks[[p]], "regions") <- NULL
		attr(landmarks[[p]], "cdists") <- NULL
		hasPeak <- !is.na(landmarks[[p]])
		for(n in 1:ncol(landmarks[[p]]))
		{
			nar <- is.na(landmarks[[p]][,n])
			landmarks[[p]][nar,n] <- mean(landmarks[[p]][,n], na.rm=TRUE)
			reg<-sapply(regions[!nar],function(r) if(is.null(dim(r))) r else r[n,])
			if(class(reg)=="list"){
				reg<-do.call(cbind,reg)
			}
			m <- matrix(apply(reg,1, mean, na.rm=TRUE), ncol=2)
			for(i in names(which(nar)))
				regions[[i]][n,] <- m
		}
		
		## register the densities
		if(ncol(landmarks[[p]])==1){ ## only one peak: offset
			if(is.null(target)){
				offsets <- landmarks[[p]]-median(landmarks[[p]])
				names(offsets)<-sampleNames(expData)
			}else{				
				offsets<-landmarks[[p]]-landmarks[[p]][sampleNames(expData)%in%target]
				names(offsets)<-sampleNames(expData)
			}
			funs <- funsBack <- vector("list", length(samples))
			names(funs)<-samples
			names(funsBack)<-samples
			for(j in seq_along(funs)){
				funs[[samples[[j]]]] <- function(x) x - z
				e1 <- new.env(hash=TRUE)
				e1$z <- offsets[samples[[j]]]
				environment(funs[[samples[[j]]]]) <- e1
				funsBack[[samples[[j]]]] <- function(x) x + z
				e2 <- new.env(hash=TRUE)
				e2$z <- offsets[samples[[j]]]
				environment(funsBack[[samples[[j]]]]) <- e2
			}
		}else{ ## multiple peaks: warping
			if(is.null(target)){
				cat("Estimating warping functions\n")
				capture.output(regDens <- landmarkreg(fdobj, landmarks[[p]], WfdPar=WfdPar, monwrd=monwrd,...))  
			}else{
				cat("Estimating warping functions\n")
				capture.output(regDens <- landmarkreg(fdobj, landmarks[[p]],x0marks= apply(landmarks[[p]],2,jitter)[rownames(landmarks)%in%target,], WfdPar=WfdPar, monwrd=monwrd,...)) 
			}
			warpfdobj <- regDens$warpfd
			warpedX <- eval.fd(warpfdobj, argvals)
			warpedX[1,] <- head(argvals,1)
			warpedX[nrow(warpedX),] <- tail(argvals,1)
			## compute warping functions
			## funs <-  apply(warpedX, 2, function(y) approxfun(argvals, y))
			funs <-  apply(warpedX, 2, approxfun, argvals)
			funsBack <- apply(warpedX, 2, function(a, b) approxfun(b, a), argvals)
		}
		names(funs) <- names(funsBack) <- sampleNames(x)
		if(!warpFuns)
		{
			warpedLandmarks <- landmarks[[p]]
			leftBoard <- rightBoard <- vector("list",length(funs))
			#chunkleftBoard<-chunkrightBoard<-rep(list(length(funs)),max(1:length(funs)%/%chunksize)+1)
			newRange <- c(Inf, -Inf)
			## transform the raw data using the warping functions
			##TODO there may be some speed up to be gained here by using chunks.
			chunkgroups<-(1:length(funs))%/%chunksize
			chunkfuns<-split(funs,chunkgroups)
			chunksamples<-split(samples,chunkgroups)
			chunkranges<-split(ranges,chunkgroups)
			chunkleftBoard<-split(leftBoard,chunkgroups)
			chunkrightBoard<-split(rightBoard,chunkgroups)
			chunkindices<-split(1:length(funs),chunkgroups)
			for(k in seq_along(chunkfuns)){
				#read in the data first for the current chunk
				if(exists("thiDat")){
					rm(thisDat)
					gc(chunksize)
				}
				thisDat <- unlist(fsApply((expData[chunksamples[[k]]]),function(z)list(exprs(z)[,,drop=FALSE])),recursive=FALSE)
				#thisDat is a list
				newDat<-vector("list",chunksize)
				message("normalizing chunk ",k," of ",length(chunkfuns));
				for(i in seq_along(chunkfuns[[k]])){
					#thisDat <- try(exprs(expData[[samples[i]]][,p]))
					lb <- thisDat[[i]][,p] < chunkranges[[k]][[i]][1,p]+eps
					lb[is.na(lb)] <- TRUE
					chunkleftBoard[[k]][[i]] <- lb
					rb <- thisDat[[i]][,p] > chunkranges[[k]][[i]][2,p]-eps
					rb[is.na(rb)] <- TRUE   
					chunkrightBoard[[k]][[i]] <- rb
					#Don't exclude things beyond the range 
					#sel <- leftBoard[[i]] | rightBoard[[i]]
					sel<-rep(FALSE,length(thisDat[[i]][,p]))
					newDat[[i]] <- as.matrix(chunkfuns[[k]][[chunksamples[[k]][[i]]]](thisDat[[i]][!sel,p]))
					newDat[[i]][is.na(newDat[[i]])] <- thisDat[[i]][!sel,p][is.na(newDat[[i]])]
					#exprs(expData[[i]])[!sel,p] <- newDat
					#ind1<-ncdfFlow:::getIndices(expData,chunksamples[[k]][i])
					#ncdfFlow::updateIndices(expData,chunksamples[[k]][i],z=NA)
					#slice<-exprs(expData[[samples[i]]])
					thisDat[[i]][,p]<-newDat[[i]]
					#slice[ind1,][!sel,p]<-newDat
					
					warpedLandmarks[chunkindices[[k]][i], ] <- chunkfuns[[k]][[chunksamples[[k]][[i]]]](landmarks[[p]][chunkindices[[k]][i],])
					newRange[1] <- chunkranges[[k]][[i]][,p][1]
					newRange[2] <- chunkranges[[k]][[i]][,p][2]
					#this doesn't work for ncdf flow sets in this scenario since it negatively impacts plotting
					#newRange[1] <- min(newRange[1], min(exprs(expData[[samples[i]]])[,p], na.rm=TRUE))
					#newRange[2] <- max(newRange[2], max(exprs(expData[[samples[i]]])[,p], na.rm=TRUE))
				}
				#Writing after normalization
				cat("Writing to disk..\n");
				for(i in seq_along(chunkfuns[[k]])){
#					browser()
					curSampleName<-chunksamples[[k]][i]
					print(i)
					print(curSampleName)
					exprs(expData[[curSampleName]])<-thisDat[[i]]
					gc()
				}
				## make sure that edge envents are set to the extreme values
				## of the warped data range and update the parameters slot
				## accordingly
				for(i in seq_along(chunkfuns[[k]])){
					#samples<-names(chunkfuns[[k]])
					minSel <- chunkleftBoard[[k]][[i]]
					maxSel <- chunkrightBoard[[k]][[i]]
					
					ip <- match(p, pData(parameters(expData[[chunksamples[[k]][i]]]))$name)
					tmp <- parameters(expData[[chunksamples[[k]][i]]])
					oldRanges <- unlist(range(expData[[chunksamples[[k]][i]]],p))
					pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], newRange[1]),
							max(oldRanges[2], newRange[2]))
					expData[[chunksamples[[k]][i]]]@parameters <- tmp
					gc()
				}
				
			}
			lm[[p]] <- list(prior=landmarks[[p]], warped=warpedLandmarks,
					warpFun=funs, regions=regions, dists=dists,
					hasPeak=hasPeak, revWarpFuns=funsBack)
		}
	}
	if(warpFuns)
		return(funs)
	#regSet <- as(expData, "flowSet")
	#regSet<-expData
	phenoData(expData) <- phenoData(x)
	expDat <- expData[sampleNames(x)]
	attr(expData, "warping") <- lm
	rm(thisDat);rm(newDat);rm(x);
	gc()
	expData
}

# When isNew == FALSE, the original cdf is modified 
# when isNew == TRUE, a new cdf is created
warpSetNCDF <- function(x, stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
		peakNr=NULL, clipRange=0.01, nbreaks=11, fres, bwFac=2,
		warpFuns=FALSE,target=NULL,chunksize=10,isNew=FALSE,newNcFile=NULL,
		...)
{
     
	## Some type-checking first
	flowCore:::checkClass(x, "flowSet")
	flowCore:::checkClass(stains, "character")
	mt <- stains %in% colnames(x)
	if(!all(mt))
		stop("Invalid stain(s) not matching the flowSet:\n    ",
				paste(stains[!mt], collapse=", "))
    

	#expData should now be x...
	if(isNew){
		expData<-ncdfFlow:::clone.ncdfFlowSet(x,isNew=TRUE,isEmpty=FALSE,ncdfFile=newNcFile)
	}else{
		expData<-x;
	}
	samples<-sampleNames(expData)
	ranges<-vector("list",length(x));
	names(ranges)<-sampleNames(x)
	for(i in 1:length(x)){ranges[[i]]<-range(x[[i]])}
	#ranges <- fsApply(x, range) 
	if(!is.null(grouping))
		flowCore:::checkClass(grouping, "character", 1)
	if(!is.null(subsample))
	{
		flowCore:::checkClass(subsample, "numeric", 1)
		## subsample data set before all density estimation steps
		#Does Subset clobber anything in a permanent way for an ncdfFlowSet? Check this, otherwise use subsample.
		x <- Subset(x, sampleFilter(size=subsample))
	}
	flowCore:::checkClass(monwrd, "logical", 1)
	flowCore:::checkClass(bwFac, "numeric", 1)
	
	## find landmarks
	if(missing(fres))
	{
		fres <- list()
		for(p in stains){
			cat("\rEstimating landmarks for channel", p, "...")
			fres[[p]] <- filter(x[,p], curv1Filter(p, bwFac=bwFac))
		}
		cat("\n")    
	}
	else
	{
		if(!is.list(fres) || !all(stains %in% names(fres)))
			stop("The supplied list of filter results does not match ",
					"the channels to warp.")
		if(!all(sapply(fres, is, "filterResult")))
			stop("The argument 'fres' has to be a list of filterResultLists.")
		cat("Extracting landmarks from user supplied filterResults\n")
	}
	
	## define some variables
	nb <- 1001
	lm <- list()
	z <- NULL 
	
	## iterate over stains
	eps <- .Machine$double.eps
	for(p in stains)
	{
      thisX <- x[,p]
		## set up fda parameters
		extend <- 0.15
		from <- min(sapply(ranges, function(z) z[1,p]-diff(z[,p])*extend), na.rm=TRUE)
		to <- max(sapply(ranges, function(z) z[2,p]+diff(z[,p])*extend), na.rm=TRUE)
		wbasis <- create.bspline.basis(rangeval=c(from, to),
				norder=4, breaks=seq(from, to, len=nbreaks))
		WfdPar <- fdPar(wbasis, 1, 1e-4)
		densY <- t(fsApply(thisX, function(z){
							r <- range(z)[,p]
							z <- exprs(z)
							z <- z[z[,p]>r[1]+eps & z[,p]<r[2]-eps, p]
							density(z, from=from, to=to, n=nb, na.rm=TRUE)$y
						}))
		argvals <- seq(from, to, len=nb) 
		fdobj   <- Data2fd(argvals,densY, wbasis )
#                           argnames = c("x", "samples", "density"))
		
		## create matrix of landmarks from curv1Filter peaks
		cat("Registering curves for parameter", p, "...\n")
		landmarks <- landmarkMatrix(thisX, fres, p, border=clipRange, peakNr=peakNr,
				densities=densY, n=nb)
		if(inherits(landmarks,"logical")){
			if(landmarks==FALSE){
				#TODO return the unnormalized data
				message("Can't detect any significant landmarks. No normalization performed.")
				return(x);
			}
		}
#		rm(fres);gc(); #shouldn't remove it since next stain still needs it
		## check if we remove signal between groups
		sig <- 0.05
		if(!is.null(grouping)){
			if(!grouping %in% names(pData(thisX)))
				stop("'", grouping, "' is not a phenoData variable.")
			grps <- as.factor(pData(thisX)[,grouping])
			anv <- numeric(ncol(landmarks))
			for(i in seq_len(ncol(landmarks)))
				anv[i] <- anova(lm(landmarks[,i] ~ grps))$Pr[1]
			if(any(anv < sig))
				warning("Within-group variances are smaller than ",
						"across-group variances for stain ", p,
						".\nWarping might have removed signal!")
		}
		## fill NAs with column medians
		regions <- attr(landmarks, "regions")
		dists <- attr(landmarks, "cdists") 
		attr(landmarks, "regions") <- NULL
		attr(landmarks, "cdists") <- NULL
		hasPeak <- !is.na(landmarks)
		for(n in 1:ncol(landmarks))
		{
			nar <- is.na(landmarks[,n])
			landmarks[nar,n] <- mean(landmarks[,n], na.rm=TRUE)
			reg<-sapply(regions[!nar],function(r) if(is.null(dim(r))) r else r[n,])
			if(class(reg)=="list"){
				reg<-do.call(cbind,reg)
			}
			m <- matrix(apply(reg,1, mean, na.rm=TRUE), ncol=2)
			for(i in names(which(nar)))
				regions[[i]][n,] <- m
		}
		
		## register the densities
		if(ncol(landmarks)==1){ ## only one peak: offset
			if(is.null(target)){
				offsets <- landmarks-median(landmarks)
				names(offsets)<-sampleNames(expData)
			}else{				
				offsets<-landmarks-landmarks[sampleNames(expData)%in%target]
				names(offsets)<-sampleNames(expData)
			}
			funs <- funsBack <- vector("list", length(landmarks))
			names(funs)<-samples
			names(funsBack)<-samples
			for(j in seq_along(funs)){
				funs[[samples[[j]]]] <- function(x) x - z
				e1 <- new.env(hash=TRUE)
				e1$z <- offsets[samples[[j]]]
				environment(funs[[samples[[j]]]]) <- e1
				funsBack[[samples[[j]]]] <- function(x) x + z
				e2 <- new.env(hash=TRUE)
				e2$z <- offsets[samples[[j]]]
				environment(funsBack[[samples[[j]]]]) <- e2
			}
		}else{ ## multiple peaks: warping
			if(is.null(target)){
				capture.output(regDens <- landmarkreg(fdobj, landmarks, WfdPar=WfdPar, monwrd=monwrd,...))  
			}else{
				capture.output(regDens <- landmarkreg(fdobj, landmarks,x0marks= apply(landmarks,2,jitter)[rownames(landmarks)%in%target,], WfdPar=WfdPar, monwrd=monwrd,...)) 
			}
			warpfdobj <- regDens$warpfd
			warpedX <- eval.fd(warpfdobj, argvals)
			warpedX[1,] <- head(argvals,1)
			warpedX[nrow(warpedX),] <- tail(argvals,1)
			## compute warping functions
			## funs <-  apply(warpedX, 2, function(y) approxfun(argvals, y))
			funs <-  apply(warpedX, 2, approxfun, argvals)
			funsBack <- apply(warpedX, 2, function(a, b) approxfun(b, a), argvals)
		}
		names(funs) <- names(funsBack) <- sampleNames(thisX)
		if(!warpFuns)
		{
			warpedLandmarks <- landmarks
			leftBoard <- rightBoard <- vector("list",length(funs))
			#chunkleftBoard<-chunkrightBoard<-rep(list(length(funs)),max(1:length(funs)%/%chunksize)+1)
			newRange <- c(Inf, -Inf)
			## transform the raw data using the warping functions
			##TODO there may be some speed up to be gained here by using chunks.
			chunkgroups<-(1:length(funs))%/%chunksize
			chunkfuns<-split(funs,chunkgroups)
			chunksamples<-split(samples,chunkgroups)
			chunkranges<-split(ranges,chunkgroups)
			chunkleftBoard<-split(leftBoard,chunkgroups)
			chunkrightBoard<-split(rightBoard,chunkgroups)
			chunkindices<-split(1:length(funs),chunkgroups)
			for(k in seq_along(chunkfuns)){
				#read in the data first for the current chunk
#				if(exists("thiDat")){
#					rm(thisDat)
#					gc(chunksize)
#				}
                
                thisChunksample <- chunksamples[[k]]
                thisChunkFuns <- chunkfuns[[k]]
                thisChunkRanges <- chunkranges[[k]]
#                thisChunkleftBoard <- chunkleftBoard[[k]]
#                thisChunkrightBoard <- chunkrightBoard[[k]]
                thischunkindices <- chunkindices[[k]]
                
				
				for(i in seq_along(thisChunkFuns)){
					message("normalizing sample ",(k-1)*chunksize+i);
					
                    curChunkRange <- thisChunkRanges[[i]]
#                    curChunkleftBoard <- thisChunkleftBoard[[i]]
#                    curChunkrightBoard <- thisChunkrightBoard[[i]]
                    curChunksample <- thisChunksample[[i]]
                    curChunkFun <- thisChunkFuns[[curChunksample]]
                    curChunkindices <- thischunkindices[[i]]
                    
                    curfr <- expData[,p][[curChunksample]]
                    curData <- exprs(curfr)[,,drop = TRUE]
                    
					lb <- curData < curChunkRange[1,p]+eps
					lb[is.na(lb)] <- TRUE
#                    curChunkleftBoard <- lb
					rb <- curData > curChunkRange[2,p]-eps
					rb[is.na(rb)] <- TRUE   
#                    curChunkrightBoard <- rb
					#Don't exclude things beyond the range 
					#sel <- leftBoard[[i]] | rightBoard[[i]]
					sel<-rep(FALSE,length(curData))
					newDat <- curChunkFun(curData[!sel])
					newDat[is.na(newDat)] <- curData[!sel][is.na(newDat)]
					
					warpedLandmarks[curChunkindices, ] <- curChunkFun(landmarks[curChunkindices,])
					newRange[1] <- curChunkRange[,p][1]
					newRange[2] <- curChunkRange[,p][2]
                    
					
				
    				#Writing after normalization
                    exprs(curfr)[,1] <- newDat
                    expData[[curChunksample,only.exprs = TRUE]] <- curfr
    				## make sure that edge envents are set to the extreme values
    				## of the warped data range and update the parameters slot
    				## accordingly
				    srcFr <- expData@frames[[curChunksample]]
#					minSel <- curChunkleftBoard
#					maxSel <- curChunkrightBoard
					tmp <- parameters(srcFr)
					oldRanges <- unlist(range(curfr))
                    ip <- match(p, pData(tmp)[,"name"])
					pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], newRange[1])
                                                                    ,max(oldRanges[2], newRange[2])
                                                                    )
                    parameters(expData@frames[[curChunksample]]) <- tmp
				}
				
			}
			lm[[p]] <- list(prior=landmarks, warped=warpedLandmarks,
					warpFun=funs, regions=regions, dists=dists,
					hasPeak=hasPeak, revWarpFuns=funsBack)
		}
	}
	if(warpFuns)
		return(funs)
	expData
}
#
warpSet <- function(x, stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
		peakNr=NULL, clipRange=0.01, nbreaks=11, fres, bwFac=2,
		warpFuns=FALSE,target=NULL,
		...)
{
	## Some type-checking first
	flowCore:::checkClass(x, "flowSet")
	flowCore:::checkClass(stains, "character")
	mt <- stains %in% colnames(x)
	if(!all(mt))
		stop("Invalid stain(s) not matching the flowSet:\n    ",
				paste(stains[!mt], collapse=", "))
	expData <- as(x, "list") #Delete for ncdfFlowSet
	ranges <- fsApply(x, range) 
	if(!is.null(grouping))
		flowCore:::checkClass(grouping, "character", 1)
	if(!is.null(subsample))
	{
		flowCore:::checkClass(subsample, "numeric", 1)
		## subsample data set before all density estimation steps
		x <- Subset(x, sampleFilter(size=subsample))
	}
	flowCore:::checkClass(monwrd, "logical", 1)
	flowCore:::checkClass(bwFac, "numeric", 1)
	
	## find landmarks
	if(missing(fres))
	{
		fres <- list()
		for(p in stains){
			cat("\rEstimating landmarks for channel", p, "...")
			fres[[p]] <- filter(x, curv1Filter(p, bwFac=bwFac))
		}
		cat("\n")    
	}
	else
	{
		if(!is.list(fres) || !all(stains %in% names(fres)))
			stop("The supplied list of filter results does not match ",
					"the channels to warp.")
		if(!all(sapply(fres, is, "filterResult")))
			stop("The argument 'fres' has to be a list of filterResultLists.")
		cat("Extracting landmarks from user supplied filterResults\n")
	}
	
	## define some variables
	nb <- 1001
	lm <- list()
	z <- NULL 
	
	## iterate over stains
	eps <- .Machine$double.eps
	for(p in stains)
	{
		## set up fda parameters
		extend <- 0.15
		from <- min(sapply(ranges, function(z) z[1,p]-diff(z[,p])*extend), na.rm=TRUE)
		to <- max(sapply(ranges, function(z) z[2,p]+diff(z[,p])*extend), na.rm=TRUE)
#		browser()
		wbasis <- create.bspline.basis(rangeval=c(from, to),
				norder=4, breaks=seq(from, to, len=nbreaks))
		WfdPar <- fdPar(wbasis, 1, 1e-4)
		densY <- t(fsApply(x, function(z){
							r <- range(z)[,p]
							z <- exprs(z)
							z <- z[z[,p]>r[1]+eps & z[,p]<r[2]-eps, p]
							density(z, from=from, to=to, n=nb, na.rm=TRUE)$y
						}))
		argvals <- seq(from, to, len=nb) 
		fdobj   <- Data2fd( argvals,densY, wbasis) 
				#argnames = c("x", "samples", "density"))
		
		## create matrix of landmarks from curv1Filter peaks
		cat("Registering curves for parameter", p, "...\n")
		landmarks <- landmarkMatrix(x, fres, p, border=clipRange, peakNr=peakNr,
				densities=densY, n=nb)
		#check if the landmarks are valid (sometimes none are significant)
		if(inherits(landmarks,"logical")){
			if(landmarks==FALSE){
				#TODO return the unnormalized data
				message("Can't detect any significant landmarks. No normalization performed.")
				return(x);
			}
		}
		## check if we remove signal between groups
		sig <- 0.05
		if(!is.null(grouping)){
			if(!grouping %in% names(pData(x)))
				stop("'", grouping, "' is not a phenoData variable.")
			grps <- as.factor(pData(x)[,grouping])
			anv <- numeric(ncol(landmarks))
			for(i in seq_len(ncol(landmarks)))
				anv[i] <- anova(lm(landmarks[,i] ~ grps))$Pr[1]
			if(any(anv < sig))
				warning("Within-group variances are smaller than ",
						"across-group variances for stain ", p,
						".\nWarping might have removed signal!")
		}
		## fill NAs with column medians
		regions <- attr(landmarks, "regions")
		dists <- attr(landmarks, "cdists") 
		attr(landmarks, "regions") <- NULL
		attr(landmarks, "cdists") <- NULL
		hasPeak <- !is.na(landmarks)
		for(n in 1:ncol(landmarks))
		{
			nar <- is.na(landmarks[,n])
			landmarks[nar,n] <- mean(landmarks[,n], na.rm=TRUE)
			reg<-sapply(regions[!nar],function(r) if(is.null(dim(r))) r else r[n,])
			if(class(reg)=="list"){
				reg<-do.call(cbind,reg)
			}
			m <- matrix(apply(reg,1, mean, na.rm=TRUE), ncol=2)
			for(i in names(which(nar)))
				regions[[i]][n,] <- m
		}
		
		## register the densities
		if(ncol(landmarks)==1){ ## only one peak: offset
			if(is.null(target)){
				offsets <- landmarks-median(landmarks)
			}else{
				#TODO BUG names(expData doesn't match target name)
				offsets<-landmarks-landmarks[names(expData)%in%target]
			}
			funs <- funsBack <- vector("list", length(landmarks))
			for(j in seq_along(funs)){
				funs[[j]] <- function(x) x - z
				e1 <- new.env(hash=TRUE)
				e1$z <- offsets[j]
				environment(funs[[j]]) <- e1
				funsBack[[j]] <- function(x) x + z
				e2 <- new.env(hash=TRUE)
				e2$z <- offsets[j]
				environment(funsBack[[j]]) <- e2
			}
		}else{ ## multiple peaks: warping
			if(is.null(target)){
				capture.output(regDens <- landmarkreg(fdobj, landmarks, WfdPar=WfdPar, 
								monwrd=monwrd, ...))  
			}else{
				
				#add a small amount of noise 1% of sd (robust) to the target landmarks
				capture.output(regDens <- landmarkreg(fdobj, landmarks, x0marks=apply(landmarks,2,jitter)[rownames(landmarks)%in%target,],WfdPar=WfdPar, 
								monwrd=monwrd, ...))
			}
			warpfdobj <- regDens$warpfd
			warpedX <- eval.fd(warpfdobj, argvals)
			warpedX[1,] <- head(argvals,1)
			warpedX[nrow(warpedX),] <- tail(argvals,1)
			## compute warping functions
			## funs <-  apply(warpedX, 2, function(y) approxfun(argvals, y))
			funs <-  apply(warpedX, 2, approxfun, argvals)
			funsBack <- apply(warpedX, 2, function(a, b) approxfun(b, a), argvals)
		}
		names(funs) <- names(funsBack) <- sampleNames(x)
		if(!warpFuns)
		{
			warpedLandmarks <- landmarks
			leftBoard <- rightBoard <- list(length(funs))
			newRange <- c(Inf, -Inf)
			## transform the raw data using the warping functions
			for(i in seq_along(funs)){
				thisDat <- exprs(expData[[i]][,p])
				lb <- thisDat < ranges[[i]][1,p]+eps
				lb[is.na(lb)] <- TRUE
				leftBoard[[i]] <- lb
				rb <- thisDat > ranges[[i]][2,p]-eps
				rb[is.na(rb)] <- TRUE   
				rightBoard[[i]] <- rb
				#Include ALL data, none of this thresholding crap at borders.
				sel <- leftBoard[[i]] | rightBoard[[i]]
				#sel<-rep(FALSE,length(thisDat))
				#browser();
				newDat <- as.matrix(funs[[i]](thisDat[!sel,]))
				newDat[is.na(newDat)] <- thisDat[!sel,][is.na(newDat)]
				exprs(expData[[i]])[!sel,p] <- newDat
				warpedLandmarks[i, ] <- funs[[i]](landmarks[i,])
				newRange[1] <- min(newRange[1], min(exprs(expData[[i]])[,p], na.rm=TRUE))
				newRange[2] <- max(newRange[2], max(exprs(expData[[i]])[,p], na.rm=TRUE))
			}
			## make sure that edge envents are set to the extreme values
			## of the warped data range and update the parameters slot
			## accordingly
			for(i in seq_along(funs)){
				minSel <- leftBoard[[i]]
				maxSel <- rightBoard[[i]]
				exprs(expData[[i]])[minSel,p] <- as.matrix(rep(newRange[1],
								sum(minSel, na.rm=TRUE)),
						ncol=1)
				exprs(expData[[i]])[maxSel,p] <- as.matrix(rep(newRange[2],
								sum(maxSel, na.rm=TRUE)),
						ncol=1)
				ip <- match(p, pData(parameters(expData[[i]]))$name)
				tmp <- parameters(expData[[i]])
				oldRanges <- unlist(range(expData[[i]],p))
				pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], newRange[1]),
						max(oldRanges[2], newRange[2]))
				expData[[i]]@parameters <- tmp
			}
			lm[[p]] <- list(prior=landmarks, warped=warpedLandmarks,
					warpFun=funs, regions=regions, dists=dists,
					hasPeak=hasPeak, revWarpFuns=funsBack)
		}
	}
	if(warpFuns)
		return(funs)
	regSet <- as(expData, "flowSet")
	phenoData(regSet) <- phenoData(x)
	regSet <- regSet[sampleNames(x)]
	attr(regSet, "warping") <- lm
	regSet
}





## Some QA pots for a normalized data set
normQA <- function(data, morph=c("^fsc", "^ssc"), channels,
		odat=NULL, ask=names(dev.cur())!="pdf",
		grouping=NULL, tag.outliers=FALSE,
		peaksOnly=TRUE)
{
	if(! "warping" %in% names(attributes(data)))
		stop("This flowSet has not been normalized.")
	fsc <- grep(morph[1], colnames(data), ignore.case=TRUE)
	if(length(fsc)) fsc <- colnames(data)[min(fsc)] else
		stop("No parameter matching regular expression '", morph[1], "'.\n",
				"Please provide a valid backgating channel.")
	ssc <- grep(morph[2], colnames(data), ignore.case=TRUE)
	if(length(ssc)) ssc <- colnames(data)[min(ssc)] else
		stop("No parameter matching regular expression '", morph[2], "'\n",
				"Please provide a valid backgating channel.")
	ninfo <- attr(data, "warping")
	wchans <- if(missing(channels)) names(ninfo) else channels
	if(!all(wchans %in% names(ninfo)))
		stop("No normalization information available for one or more channels.")
	if(is.null(grouping))
	{
		grouping <- factor(rep(1, length(data)))
	}
	else
	{
		if(is.character(grouping)){
			if(! grouping %in% names(pData(data)))
				stop("'grouping' must be a covariate in the phenoData slot of ",
						"'data' or a factor of the same length as 'data'.")
			grouping <- factor(pData(data)[,grouping])
		}
		else
		{
			grouping <- factor(grouping)
		}
		if(length(grouping) != length(data))
			stop("'grouping' must be a covariate in the phenoData slot of ",
					"'data' or a factor of the same length as 'data'.")
	}
	
	## Plot the amount of warping for each landmark.
	if(all(c("prior", "warped", "hasPeak") %in% names(ninfo[[1]])))
	{
		## We first create a dataFrame ameanable for trellis graphics.
		prior <- data.frame()
		m <- lapply(ninfo, function(x) apply(x$warped, 2, mean, na.rm=TRUE))
		for(p in wchans)
		{
			pr <- ninfo[[p]]$prior
			hp <- ninfo[[p]]$hasPeak
			nc <- ncol(pr)
			nr <- nrow(pr)
			peak <- factor(rep(seq_len(nc), each=nr))
			present <- sapply(split(hp, peak), function(x)
						sprintf("(%d/%d)", sum(x), nr))
			peakLong <- factor(paste(peak, present[peak]))
			prior <- rbind(prior,
					data.frame(sample=factor(sampleNames(data),
									levels=sampleNames(data)),
							group=factor(as.integer(grouping)),
							value=as.vector(pr),
							peak=peak, peakLong=peakLong,
							channel=factor(p),
							hasPeak=as.vector(hp),
							dists=if(length(ninfo[[p]]$dists))
										1-as.vector(ninfo[[p]]$dists) else NA))
			sel <- is.na(prior$value)
			prior[sel, "value"] <- m[[p]][prior[sel, "peak"]]
		}
		## A regular stripplot with dotted lines extending from the
		## adjusted landmark positions for each peak.
		par(ask=ask)
		on.exit(par(ask=FALSE))
		myPanelStrip <- function(x,y, groupMeans, datSet, groups, ng, ...)
		{
			chan <- levels(datSet$channel)[panel.number()]
			thisP <- datSet[datSet$channel==chan,] 
			nc <- length(groupMeans[[chan]])
			nr <- nlevels(y)
			col <- rep(trellis.par.get("superpose.symbol")$col[seq_len(ng)], ng)
			pch <- if(ng==1) rep(c(4,20), 2) else rep(c(4, 20), each=ng)
			panel.abline(v=groupMeans[[chan]], lty="dotted", col="gray")
			for(i in seq_len(nc)){
				peak <- thisP[thisP$peak==i,"value"]
				panel.segments(y0=seq_len(nr), y1=seq_len(nr), x0=peak,
						x1=groupMeans[[chan]][i], col="gray", lty="dotted")
			}
			panel.stripplot(x,y, pch=pch, col=col, groups=groups, ...)
		}
		col <- trellis.par.get("superpose.symbol")$col[seq_len(nlevels(prior$group))]
		ng <- nlevels(prior$group)
		key <- list()
		if(ng>1)
			key <- append(key, list(points=list(col=col, pch=20),
							text=list(paste("Group", levels(grouping)))))
		if(!all(prior$hasPeak))
			key <- append(key, list(points=list(pch=4, col=1),
							text=list("no peak detected"),
							rep=FALSE))
		print(stripplot(sample~value|channel, prior, ng=ng,
						groups=interaction(group, hasPeak),
						panel=myPanelStrip, groupMeans=m, datSet=prior, xlab=NULL,
						key=if(length(key)) c(key, cex=0.8, between=1) else NULL,
						main="Amount of Landmark Adjustment\n"))
	}
	
	## Plot the confidence of the landmark registration step
	if("dists" %in% names(ninfo[[1]]))
	{
		present <- lapply(split(prior, prior$channel), function(x)
					sapply(split(x$hasPeak, x$peak), function(x)
								sprintf("(%d/%d)", sum(x), length(data))))
		col <- trellis.par.get("superpose.symbol")$col[seq_len(nlevels(prior$group))]
		key <- list()
		ng <- nlevels(prior$group)
		if(ng>1)
			key <- append(key, list(text=list(paste("Group", levels(grouping))),
							points=list(col=col, pch=1)))
		print(stripplot(dists~peak |channel, prior, ylim=c(0,1.1),
						ylab="Clustering Confidence", xlab="Peak",
						groups=group, myLabel=present,
						panel=function(myLabel, datSet, ...){
							panel.stripplot(...)
							chan <- levels(datSet$channel)[panel.number()]
							panel.text(seq_along(myLabel[[chan]]), 1.05, myLabel[[chan]],
									cex=0.7)
							panel.abline(h=1, col="gray")
						}, datSet=prior,
						key=if(length(key)) c(key, cex=0.8) else NULL,    
						main="Landmark Registration\n"))
	}
	
	## Plot the warping functions. 
	if("warpFun" %in% names(ninfo[[1]]))
	{
		## We first create a dataFrame for the trellis plotting
		wfRes <- data.frame()
		for(p in wchans)
		{
			np <- 30
			xvals <- seq(range(data[[1]], p)[1,], range(data[[1]], p)[2,], len=np)
			yvals <- sapply(ninfo[[p]]$warpFun, function(fun) fun(xvals))
			wfRes <- rbind(wfRes,
					data.frame(sample=rep(factor(sampleNames(data),
											levels=sampleNames(data)),
									each=np),
							group=rep(factor(as.integer(grouping)),each=np),
							original=xvals, normalized=as.vector(yvals),
							channel=p))
		}
		col <- trellis.par.get("superpose.symbol")$col[seq_len(nlevels(wfRes$group))]
		key <- list()
		ng <- nlevels(wfRes$group)
		if(ng>1)
			key <- append(key, list(text=list(paste("Group", levels(grouping))),
							lines=list(col=col, pch=1)))
		print(xyplot(normalized~original|channel, wfRes, type="l", alpha=0.4,
						main="Warping Functions\n", groups=group,
						key=if(length(key)) c(key, cex=0.8) else NULL))
	}
	
	
	## Plot the projection in the FCS vs. SSC space.
	if("regions" %in% names(ninfo[[1]]))
	{
		## We first create rectangleGates from the peak regions and
		## compute norm2Filters in the backgating dimensions
		cat("Computing the backgating information. Please wait...")
		backGates <- cents <- vector(mode="list", length=length(wchans))
		names(backGates) <- names(cents) <-  wchans
		stand <- apply(range(data[[1]])[,c(fsc, ssc)], 2, diff)
		vars <- data.frame()
		noOdat <- !is.null(ninfo[[1]]$warpFun) && is.null(odat)
		if(!noOdat && is.null(odat))
			stop("The un-nomalized data set needs to be supplied as argument 'odat'.")
		if(noOdat)
			odat <- flowCore:::copyFlowSet(data)
		for(p in wchans)
		{
			regions <- ninfo[[p]]$regions
			if(is.null(names(regions)))
				names(regions) <- sampleNames(data)
			wfuns <- if(noOdat) ninfo[[p]]$warpFun else {
						wf <- list()
						for(n in names(regions))
							wf[[n]] <- function(x) x
						wf}
			np <- ncol(ninfo[[p]]$warped)
			bg <- cs <- vector(mode="list", length=np)
			for(i in seq_len(np))
			{
				g <- list()
				for(j in names(regions))
				{
					g[[j]] <- if(ninfo[[p]]$hasPeak[which(sampleNames(data)==j),i])
								rectangleGate(matrix(wfuns[[j]](regions[[j]][i,]), ncol=1,
												dimnames=list(NULL, p)))
							else rectangleGate(matrix(c(-Inf, Inf),  ncol=1,
												dimnames=list(NULL, p)))
				}
				subDat <- Subset(Subset(odat, filterList(g)), boundaryFilter(c(fsc, ssc)))
				hp <- ninfo[[p]]$hasPeak[,i]
				bg[[i]] <- filter(subDat, norm2Filter(c(fsc, ssc)))
				bg[[i]][fsApply(subDat, nrow)<30] <- NA
				cs[[i]] <-  fsApply(subDat, function(x){
							if(nrow(x)>30)
								apply(x[,c(fsc, ssc)], 2, mean, na.rm=TRUE)
							else as.numeric(rep(NA,2))}, use.exprs=TRUE)
				if(peaksOnly)
				{
					bg[[i]][!hp] <- NA
					cs[[i]][!hp] <- rep(NA,2)
				}
				vars <- rbind(vars,
						data.frame(value=sapply(bg[[i]], function(x)
											if(is(x, "filterResult") && !is.na(filterDetails(x)[[1]]$cov) &&
													length(x@subSet) > 30 )
												prod(sqrt(eigen(filterDetails(x)[[1]]$cov)$values)/stand)*pi
											else NA),
								channel=p, peak=factor(paste("Peak", i)),
								sample=factor(sampleNames(data),
										levels=sampleNames(data))))
			}
			backGates[[p]] <- bg
			cents[[p]] <- cs
		}
		dummy <- data.frame(channel=factor(rep(names(backGates),
								times=length(data)*listLen(backGates))),
				peak=factor(unlist(sapply(listLen(backGates),
										function(x)
											rep(seq_len(x),
													each=length(data))))),
				sample=factor(sampleNames(data), levels=sampleNames(data)),
				group=factor(as.integer(grouping)),
				x=0, y=0)
		cat("\n")
		myXYPanel <- function(datSet, backGates, groups, ...)
		{
			wp <- which.packet()
			p <- levels(datSet$channel)[wp[1]]
			ng <- nlevels(groups)
			bg <- backGates[[p]]
			if(length(wp)==2)
				bg <- bg[wp[2]]
			for(i in seq_along(bg))
			{
				crgb <- col2rgb(trellis.par.get("superpose.symbol")$col[as.integer(groups)])
				col <- rgb(t(crgb), alpha=255*0.12, maxColorValue=255)
				colc <- rgb(t(crgb), alpha=255*0.2, maxColorValue=255)
				for(f in seq_along(bg[[i]]))
					if(is(bg[[i]][[f]], "filterResult") &&
							!is.na(filterDetails(bg[[i]][[f]])[[1]]$cov) &&
							length(bg[[i]][[f]]@subSet) > 30 )
						glpolygon(bg[[i]][[f]], gpar=list(gate=list(fill=col[f],
												col=colc[f])))
			}
		}
		flowViz:::plotType("gsmooth", c(fsc, ssc))
		key <- list()
		ng <- nlevels(dummy$group)
		col <- trellis.par.get("superpose.symbol")$col[seq_len(ng)]
		if(ng>1)
			key <- append(key, list(text=list(paste("Group", levels(grouping))),
							rectangles=list(col=col, pch=1)))
		print(xyplot(x~y|channel+factor(paste("Peak", peak)), dummy, xlab=fsc, ylab=ssc,
						xlim=range(data[[1]])[,fsc], groups=group,
						ylim=range(data[[1]])[, ssc],
						panel=myXYPanel, datSet=dummy, backGates=backGates,
						key=if(length(key)) c(key, cex=0.8) else NULL,
						main="Backgating Shape\n"))
		
		
		##  tmp <- sapply(regions, rowMeans)
		##         expval <- signif(if(!is.null(dim(tmp))) rowMeans(tmp) else mean(tmp),3)
		##         legend("topleft", pch=20, col=seq_len(nc)+1,
		##                legend=sprintf("Peak %d (mean %s=%g)", seq_len(nc), p, expval),
		##                cex=0.8, bty="n")
		
		## Now only plot the centroids of the ellipses
		alldat <- Subset(as(odat, "flowFrame")[,c(fsc, ssc)],
				boundaryFilter(c(fsc, ssc)) %subset%
						sampleFilter(10000))   
		
		myXYPanelCent <- function(datSet, centroids, alldat, labels, groups, ...)
		{
			wp <- which.packet()
			p <- levels(datSet$channel)[wp[1]]
			ct <- centroids[[p]]
			if(length(wp)==2)
				ct <- ct[wp[2]]
			exp <- exprs(alldat)[,1:2]
			xr <- range(exp[,1], na.rm=TRUE)
			yr <- range(exp[,2], na.rm=TRUE)
			bw <- diff(apply(exp, 2, quantile, probs=c(0.05, 0.95),
							na.rm=TRUE)) / 25
			range <- list(xr+c(-1,1)*bw[1]*2.5, yr+c(-1,1)*bw[2]*2.5)
			z <- bkde2D(exp, bw, c(65,65), range.x=range)
			ll <- contourLines(z$x1, z$x2, z$fhat, nlevels=10)
			for(pg in ll)
				panel.polygon(pg$x, pg$y, border="lightgray")
			col <- trellis.par.get("superpose.symbol")$col[as.integer(groups)]
			for(i in seq_along(ct))
			{
				panel.points(ct[[i]], col=col, pch=i)
				if(tag.outliers)
				{
					tg <- subset(datSet, channel==p & peak==i)$group
					x <- split(ct[[i]][,1], tg)
					y <- split(ct[[i]][,2], tg)
					outliers <- unlist(mapply(function(xx, yy){
										tmp <- cbind(xx, yy)
										sel <- apply(tmp, 1, function(z) any(is.na(z)))
										res <- rep(FALSE, nrow(tmp))
										ol <- !pcout(tmp[!sel,], outbound=0.1)$wfinal01
										res[!sel] <- ol
										return(res)}, x, y, SIMPLIFY=FALSE))
					panel.text(ct[[i]][,1][outliers], ct[[i]][,2][outliers],
							labels[outliers], adj=c(1.5,1.5), cex=0.6)
				}
			}
		}
		key <- list()
		np <- nlevels(dummy$peak)
		if(ng>1 || np>1)
		{
			labs <- as.character(levels(interaction(paste("Group", levels(grouping)),
									paste("Peak", levels(dummy$peak)),
									sep=" ")))
			key <- append(key, list(text=list(labs),
							points=list(col=rep(col, nlevels(grouping)),
									pch=rep(seq_len(nlevels(dummy$peak)),
											each=nlevels(grouping))),
							columns=ng))
		}
		print(xyplot(x~y|channel, dummy, xlab=fsc, ylab=ssc,
						xlim=range(data[[1]])[,fsc], groups=group,
						ylim=range(data[[1]])[, ssc],
						panel=myXYPanelCent, datSet=dummy, centroids=cents,
						main="Backgating Location\n", alldat=alldat, labels=sampleNames(data),
						key=if(length(key)) c(key, cex=0.8) else NULL))
		
		
		cvSum <- data.frame()
		for(i in seq_along(cents))
			for(j in seq_along(cents[[i]]))
			{
				mc <- apply(cents[[i]][[j]], 2, median, na.rm=TRUE)
				md <- t(cents[[i]][[j]]) - mc
				mds <- colMeans(md/stand)
				cvSum <- rbind(cvSum, data.frame(sample=factor(names(mds),
										levels=unique(names(mds))),
								group=factor(as.integer(grouping)),
								channel=wchans[[i]],
								peak=factor(paste("Peak",j)),
								location=mds))
			}
		cvSum <- cbind(cvSum, variation=vars[, "value"])
		myXYPanelSum <- function(x, y, labels, datSet, groups, ...)
		{
			wp <- which.packet()
			tg <- subset(datSet, channel==levels(factor(datSet$channel))[wp[1]] &
							peak==paste("Peak", wp[2]))$group
			if(length(x))
			{
				col <- trellis.par.get("superpose.symbol")$col[as.integer(tg)]
				panel.points(x,y, col=col, pch=1)
				if(tag.outliers)
				{
					xs <- split(x, tg)
					ys <- split(y, tg)
					outliers <- unlist(mapply(function(xx, yy){
										tmp <- cbind(xx, yy)
										sel <- apply(tmp, 1, function(z) any(is.na(z)))
										res <- rep(FALSE, nrow(tmp))
										ol <- !pcout(tmp[!sel,], outbound=0.1)$wfinal01
										res[!sel] <- ol
										return(res)}, xs, ys, SIMPLIFY=FALSE))
					col <- trellis.par.get("superpose.symbol")$col[as.integer(tg)]
					panel.text(x[outliers], y[outliers], labels[outliers], adj=c(1.5,1.5),
							cex=0.6)
				}
			}
		}
		key <- list()
		col <- trellis.par.get("superpose.symbol")$col[seq_len(ng)]
		if(ng>1)
			key <- append(key, list(text=list(paste("Group", levels(grouping))),
							rectangles=list(col=col, pch=1)))
		print(xyplot(location ~ variation|channel+peak, cvSum, panel=myXYPanelSum,
						labels=sampleNames(data), main="Backgating Summary\n",
						groups=group, datSet=cvSum,
						key=if(length(key)) c(key, cex=0.8) else NULL))
		
	}
	## Plot the relative ellipse volumes for the backgated channels
	##  layout(matrix(1:length(wchans), ncol=1))
	##     for(p in wchans)
	##     {
	##         nc <- ncol(ninfo[[p]]$warped)
	##         nr <- nrow(ninfo[[p]]$warped)
	##         plot(seq_len(nc), rep(0, nc),
	##              ylim=c(0,max(5, max(unlist(allV[[p]])))), type="n",
	##              xlim=c(0, nc+1), xlab="peaks", ylab="variation", xaxt="n",
	##              pch=pch, col=rep(seq_len(nc), each=nr),
	##              main=paste("Backgating Variation", p))
	##         axis(1, at=seq_len(nc), labels=TRUE)
	##         box()
	##         for(i in seq_len(nc))
	##             points(rep(i, nr), allV[[p]][[i]]) 
	##     }
	invisible()
}



























## An optimized object size function that is able to deal with embedded
## environments.
## setGeneric("objectSize", function(x, ...) standardGeneric("objectSize"))
## setMethod("objectSize", signature(x="ANY"), definition=function(x, ...)
##       {
##           realObjSize <- function(object, size=0, done=list())
##           {
##               ## recursively go though environments and record the total size
##               ## We also pass on pointers to all the environments that have
##               ## already been counted in order to avoid infinite loops or
##               ## counting things multiple times.
##               objSize <- size
##               if(is.environment(object) && !any(sapply(done, identical, object)))
##               {
##                   for(i in ls(object))
##                       objSize <- realObjSize(get(i,object),objSize,
##                                              done=append(done, object))
##               }
##               else
##               {
##                   objSize <- object.size(object, ...)+size
##               }
##               ## check if any of the slots are evironments and pass through those
##               slots <- slotNames(object)
##               envs <- sapply(slots, function(s) is.environment(slot(object, s)))
##               for(i in slots[envs])
##                   objSize <- realObjSize(slot(object,i), objSize, done)
##               return(objSize)
##           }
##           realObjSize(x)
##       })
