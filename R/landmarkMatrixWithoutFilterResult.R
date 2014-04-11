#TODO fix formal arguments
landmarkMatrixWithoutFilterResult <- function(data, peaksandregions, parm, border=0.05, peakNr=NULL,
                           indices=FALSE)
{
        ##TODO maybe code to organzie the peaksandregions variable should go here..        
        peaks<-do.call(c,lapply(peaksandregions,function(x)x[[parm]][[1]]))
        regions<-do.call(c,lapply(peaksandregions,function(x)x[[parm]][[2]]))
        
        ranges <- fsApply(data, range)
        nrPeaks <- table(listLen(peaks))
        fnrPeaks <- as.numeric(max(names(which(nrPeaks/sum(nrPeaks) > 0.1))))
        if(!is.null(peakNr))
    	peakNr <- min(as.numeric(names(which.max(nrPeaks)), peakNr))
        clustCenters <- if(fnrPeaks>1 && !is.null(peakNr)){
            fnrPeaks <- peakNr
            apply(matrix(unlist(peaks[(listLen(peaks)==peakNr)]), ncol=peakNr, byrow=T) ,2, mean)
        }else fnrPeaks
         resRegions <- list()
		 if(is.na(fnrPeaks))
			 return(FALSE)
        if(fnrPeaks==1){
            single <- which(listLen(peaks)==1)
            med <- median(unlist(peaks[single]), na.rm=TRUE)
            rmed <- apply(t(sapply(regions[single], c)), 2, median, na.rm=TRUE)
            resPeaks <- numeric(length(peaks))
            resPeaks[single] <- unlist(peaks[single])
            resRegions[single] <- regions[single]
            resPeaks[-single] <- unlist(sapply(peaks[-single],
                                               function(x) x[which.min(abs(x-med))]))
            resRegions[-single] <- mapply(function(x,y) y[which.min(abs(x-med)),,drop=FALSE],
                                          x=peaks[-single],
                                          y=regions[-single], SIMPLIFY=FALSE)
            resPeaks[is.na(resPeaks)] <- med
            resRegions[is.na(resPeaks)] <- matrix(rmed, ncol=2)
            names(resRegions) <- sampleNames(data)
            if(!indices){
                m <- matrix(resPeaks, ncol=1)
                attr(m, "regions") <- resRegions
                attr(m, "cdists") <- matrix(0, nrow=length(data), ncol=1,
                                           dimnames=list(sampleNames(data),NULL))
                return(m)
            }
            else
            {
                return(matrix(1, ncol=1, nrow=length(resPeaks),
                              dimnames=list(sampleNames(data), NULL)))
            }
        }  
        ## cluster peaks in k cluster where k is max number of peaks for a sample
        mat <- matD <- matrix(nrow=length(peaks), ncol=fnrPeaks)
        rownames(mat) <-  rownames(matD) <- sampleNames(data)
        pvect <- unlist(peaks)
        names(pvect) <- rep(names(peaks), listLen(peaks))
        sel <- !is.na(pvect)
        km <- kmeans(pvect[sel], clustCenters)
        km$cluster <- match(km$cluster, order(km$centers))
        clusterDist <- numeric(length(km$cluster))
        for(i in seq_along(km$cluster))
            clusterDist[i] <- abs(pvect[i] - km$centers[km$cluster[i]])
        cList <- split(data.frame(cluster=km$cluster, dist=clusterDist,
                           landmark=pvect[sel], index=unlist(sapply(peaks, seq_along))[sel]),
                       names(pvect[sel]))
        cList <- lapply(cList, function(x) x[order(x$dist),][1:min(nrow(x),
                                                                   fnrPeaks),])
        ## put peaks in matrix according to clustering where row=sample and
        ## col=cluster
        for(i in rownames(mat)){
            cl <- cList[[i]]
            tmp <- matrix(NA, ncol=2, nrow=fnrPeaks)
            tmp[cl$cluster,] <- regions[[i]][cl$index,,drop=FALSE]
            resRegions[[i]] <- matrix(tmp, ncol=2)
            mat[i,cl$cluster] <- if(!indices) cl$landmark else cl$index
            matD[i,cl$cluster] <- cl$dist/diff(ranges[[1]][,parm])
        }
        resRegions <- resRegions[rownames(mat)]
        attr(mat, "regions") <- resRegions
        attr(mat, "cdists") <- matD
        return(mat)
    }

