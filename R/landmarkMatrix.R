## create matrix of landmarks by aligning modes of the data density estiamtes
## Number of modes per sample doesn't need to be the same, the algorithm will
## try to guess wich of the modes are to be aligned.
## Currently this uses k-means clustering
landmarkMatrix <- function(data, fres, parm, border=0.05, peakNr=NULL, densities=NULL, n=201)
{
    ## Some type-checking first
    flowCore:::checkClass(data, "flowSet")
    flowCore:::checkClass(fres, "list")
    flowCore:::checkClass(parm, "character")
    flowCore:::checkClass(border, "numeric", 1)
    ## get peaks of significant regions from curv1Filter results
    ranges <- fsApply(data, range)
    from <- min(sapply(ranges, function(z) z[1,parm]-diff(z[,parm])*0.15), na.rm=TRUE)
    to <- max(sapply(ranges, function(z) z[2,parm]+diff(z[,parm])*0.15), na.rm=TRUE)
    peaks <- list()
    eps <- .Machine$double.eps
    for(i in sampleNames(data)){
        if(is.null(densities))
        {
            dens <- densities
            x <- exprs(data[[i]][,parm])
            r <- ranges[[i]][,parm]
            x <- x[x>r[1]+eps & x<r[2]-eps]
        }
        else
        {
            dens <- densities[,i]
            x <- NULL
        }
        peaks[[i]] <- curvPeaks(fres[[parm]][[i]], x, from=from, to=to, n=n,
                                borderQuant=border, densities=dens)[["peaks"]][,"x"]
    }
    ## are multiple peaks reasonable?
    nrPeaks <- table(listLen(peaks))
    fnrPeaks <- as.numeric(max(names(which(nrPeaks/sum(nrPeaks) > 0.1))))
    if(!is.null(peakNr))
	peakNr <- min(as.numeric(names(which.max(nrPeaks)), peakNr))
    clustCenters <- if(fnrPeaks>1 && !is.null(peakNr)){
        fnrPeaks <- peakNr
        apply(matrix(unlist(peaks[(listLen(peaks)==peakNr)]), ncol=peakNr, byrow=T) ,2, mean)
    }else fnrPeaks
    if(fnrPeaks==1){
        single <- which(listLen(peaks)==1)
        med <- median(unlist(peaks[single]), na.rm=TRUE)
        resPeaks <- numeric(length(peaks))
        resPeaks[single] <- unlist(peaks[single])
        resPeaks[-single] <- unlist(sapply(peaks[-single],
                                 function(x) x[which.min(abs(x-med))]))
        resPeaks[is.na(resPeaks)] <- med
        return(matrix(resPeaks, ncol=1))
    }  
    ## cluster peaks in k cluster where k is max number of peaks for a sample
    mat <- matrix(nrow=length(peaks), ncol=fnrPeaks)
    rownames(mat) <- sampleNames(data)
    pvect <- unlist(peaks)
    names(pvect) <- rep(names(peaks), listLen(peaks))
    sel <- !is.na(pvect)
    km <- kmeans(pvect[sel], clustCenters)
    km$cluster <- match(km$cluster, order(km$centers))
    clusterDist <- numeric(length(km$cluster))
    for(i in seq_along(km$cluster))
        clusterDist[i] <- abs(pvect[i] - km$centers[km$cluster[i]])
    cList <- split(data.frame(cluster=km$cluster, dist=clusterDist,
                       landmark=pvect[sel]), names(pvect[sel]))
    cList <- lapply(cList, function(x) x[order(x$dist),][1:min(nrow(x),
                                                               fnrPeaks),])
    ## put peaks in matrix according to clustering where row=sample and
    ## col=cluster
    for(i in names(cList)){
        cl <- cList[[i]]
        mat[i,cl$cluster] <- cl$landmark
    }
    return(mat)
}

