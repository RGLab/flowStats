## use backgeting to identify the features and define a set of reference
## features which are needed by function gpaSet

idFeatures <- function(bg, merge.peak.cluster=TRUE, thres.merge=1/5, plot=FALSE)
{
    ## 'bg' is the return value of function flowStats:::backGating

    regFeatures <- list()
    refFeatures <- list()
    ch = split(bg, factor(bg$channel))
    ## identify features on the projection of each bg channels
    for (i in names(ch))
    { 
        perChannel <- ch[[i]]
       	centInd <- which(perChannel$type=="centroid")
    	center  <- perChannel[centInd, ]
        ## cluster the centriod based on peak+population. For example,
        ## CD8 has two peaks, one peak has 3 population and one has 1 population.
        ## Then, CD8 has four cluster center
        center  <- labelCentroid(center) # add a slot: center$cluster
        
        ## let reference features = median of cluster members
        refFeatures[[i]] <- t(sapply(split(center, center$cluster),
                                function(a) apply(a[, 1:2], 2, median)))
                  
        perSample = split(center, factor(center$sample))

        ## classifyFeatures based upon the cluster labelling done by labelCluster
        regFeatures[[i]] <- 
             lapply(perSample, classifyFeatures, refFeatures[[i]],
                    by="distance")

        #if (plot) .plotResults(perChannel, center)
    }

    ## combine the members of regFeatures for each samples
   
    combf <- list()
    for (m in levels(factor(bg$sample))) {
         combf[[m]] <- do.call(rbind, 
                     lapply(regFeatures, function(a, b) {a[[b]]}, m))
    }
    
    reference <- do.call(rbind, refFeatures)
    rownames(reference) <- rownames(combf[[1]])
    reference <- as.data.frame(reference)
    reference$bogus <- rep(FALSE, nrow(reference))
    reference$cluster <- rownames(combf[[1]])
    combf <- lapply(combf, function(x)
                    {x$cluster=rownames(x); x})
    combf[["reference"]] <- reference

    
    if (merge.peak.cluster) {
        ## merge clusters and re-calculate the reference points
        combf <- mergeCluster(combf, thres=thres.merge)
        ## re-calculate the reference?
    }
    
    return(combf)
}

###################
## classifyFeatures
###################
classifyFeatures <- function(eachSample, refF, by="distance")
{
  ## make the size of refF and column names consistant with refF 
  ## initialize parameter
  reg <- data.frame(apply(refF, 2, function(x) rep(NA, length(x))),
                    bogus=rep(FALSE, nrow(refF)))

  perClust <- split(eachSample, factor(eachSample$cluster))
       for (j in 1:length(perClust))
       {
           lab <- as.numeric(perClust[[j]]$cluster)[1]
	   if (nrow(perClust[[j]]) == 1) {
       	       reg[lab, c("x", "y")] <- as.matrix(perClust[[j]][, c("x", "y") ])
	   }
	   else {
	       tmp <- rbind(refF[lab, c("x", "y")],
                     perClust[[j]][, c("x", "y")])
               mindist <-
                   which.min(as.numeric(dist(tmp))[1:nrow(perClust[[j]])])
	       reg[lab, c("x","y")] <- as.matrix(
	          perClust[[j]][mindist, c("x", "y") ])
           }
       }

  ## 2. deal with NA peaks, give the refFeatures and cluster status as bogus
    naIdx <- which(is.na(reg[, 1]))
  
    if (!identical(naIdx, integer(0))) {
       reg[naIdx, c("x", "y")] <- refF[naIdx, c("x", "y")]
       reg[naIdx, "bogus"] <- TRUE
    }
  
  return(reg)
}

#################
###############
## plot results
###############
.plotResult <- function(perChannel, center)
{
    channel = perChannel$channel[1]
    print(xyplot(x ~ y | factor(perChannel$sample), 
                 data=perChannel,
                 pch=20,
                 main=paste("backgating on", channel)))

   print(xyplot(x~y, data=center,
                main=paste("all centroids", channel)))
    
   print(xyplot(x~y|factor(center$peak),
                data=center, 
                main=paste("centriod per popupations", channel)))
            
   print(xyplot(x~y|factor(paste("peak", center$peak)) +
                  factor(paste("population", center$population)), 
                data=center, 
                main=paste("centriod per peak", channel)))
    	                
   print(xyplot(x ~ y|which,
                do.call(make.groups, lapply(regFeatures[[i]], as.data.frame)),
                as.table=TRUE, 
                main=paste("registered features with backgating on", channel)))
            
   #Features=regFeatures[[i]]
   #Features$reference=refFeatures[[i]][, 1:2]
   #print(xyplot(x ~ y,
   #             data=do.call(make.groups, lapply(Features, as.data.frame)),
   #             pch=0:(length(Features)-1), 
   #   	        group=which, col="black",
   #   	        key=list(space="right", title="samples",
   #                text=list(names(Features)), 
   #                points=list(pch=1:length(Features)-1 )),
   #   	        main=paste("registered features with backgating on", channel)))
}

useCluster <- function(center) {
  ## decide k by BIC score
  BIC <- mclustBIC(center[, 1:2])
  plot(BIC)
  
  #a <- kmeans(center[, 1:2], centers=2)
  #a <- Mclust(center[, 1:2]) #a$classification
  
  center$label <- factor(a$cluster)
  design <- model.matrix(~0+center$label)
  colnames(design) <- paste("clust", seq(ncol(design)), sep="")
  aa <- apply(design, 2, function(x, y) y[as.logical(x),], center[, 1:2])
  x11();xyplot(x~y, data=do.call(make.groups, lapply(aa, as.data.frame)),
         group=which, auto.key=list(space="right"))

}

###############################################
## labelling the centroid by kmeans clustering
###############################################
labelCentroid <- function(center) {

 
    perPeak = split(center, factor(center$peak))
    nclust <- 0
    km <- NULL
    for (i in 1:length(perPeak)) {
        ## if there are more than one population -> no clustering, get median
        if (length(unique(perPeak[[i]]$population))==1) {
            medPop <- apply(perPeak[[i]][, 1:2], 2, median)
            perPeak[[i]]$cluster <- 1 + nclust 
        }
        else {
            perPop <- split(perPeak[[i]], factor(perPeak[[i]]$population))
            medPop <- lapply(perPop, function(a) apply(a[, 1:2], 2, median))
            km <- kmeans(perPeak[[i]][, 1:2], do.call(rbind, medPop))
            perPeak[[i]]$cluster <- km$cluster+nclust
        }
        nclust <- nclust + length(medPop) 
    }

    
    mPeak <- do.call(rbind, perPeak)
    
    return(do.call(rbind, perPeak))
}

## merge cluster
#mergeCluster <- mergeCluster(mPeak)
#{
#    perClust <- split(mPeak, factor(mPeak$cluster))
#    md <- t(sapply(perClust, function(x) apply(x[, 1:2], 2, median)))
#    dst <- dist(md, diag=TRUE, upper=TRUE)
#    ## use Hierarchical cluster analysis to merge the equivalent clusters
#    ## and give label accordingly
#    hdst <- hclust(dst)
#
#}

mergeCluster <- function(feats, thres=1/5)
{
    sampleF <- feats[-which(names(feats)=="reference")]
    ref <- feats$reference
    dst <- dist(ref[, 1:2])
    hdst <- hclust(dst)
    JMAX = 7
    nTimes = 1
    #plot(hdst)
    ## while loop: continue if the ratio of the shortest and longest distance
    ## is less than 1/15
    while(nTimes <= JMAX | min(dst)/max(dst) < thres) {
        mergeIdx <- which(rowMeans(sign(hdst$merge)) < 0)
        dst <- as.matrix(dst)
        
        ## fix reference (ref) if needed
        ## fix index (mergeIdx) at which sampleF needed to be altered
        fixIdx <- NULL
        for (i in 1:length(mergeIdx)) {
            nclust <- -hdst$merge[mergeIdx[i], ]
            dratio <- dst[nclust[1], nclust[2]] / max(dst)   
            if (dratio < thres) ## replace by the mean values
                ref[nclust[1], 1:2] <-
                    (ref[nclust[1], 1:2]+ref[nclust[2], 1:2])/2
            else
                fixIdx=cbind(fixIdx, i)
        }
        if (!is.null(fixIdx)) mergeIdx <- mergeIdx[-fixIdx]
        
        ## if mergeIdx is empty, then the distances between the reference are
        ## sufficient large. -> no need to continue the loop
        if (identical(mergeIdx, integer(0))) break
            
        ## fix the sampleF: get the one that is closer to the reference features
        sampleF <- lapply(sampleF, function(sF, ref, cidx) {
                          ## make sure cidx is an matrix
                          if (!is.matrix(cidx)) cidx <- matrix(cidx, nrow=1)
                          
                          for (j in 1:nrow(cidx)) {
                            x <- cidx[j,]
                            dst <- dist(rbind(ref[x[1],], sF[x, ])[, 1:2])
                            sF[x[1], 1:2] <-
                              sF[x[which.min(as.numeric(dst)[1:2])], 1:2]
                          }
                      return(sF)
                      }, ref, -hdst$merge[mergeIdx, ])
        ## delete hdst$merge[mergeIdx, 2] from ref and sampleF
        del <- -hdst$merge[mergeIdx,2]
        ref <- ref[-del, ]
        sampleF <- lapply(sampleF, function(x, del) x[-del, ], del)

        ## get the hierarchical cluster analysis and distance
        dst <- dist(ref[, 1:2])
        hdst <- hclust(dst)
        nTimes <- nTimes + 1
    }
    sampleF[["reference"]] <- ref
    return(sampleF)
}

is.allNegative <- function(x) {
    mean(sign(x))<0

}
.plotFeatures <- function(combf, reference) {

   ## per sample
   print(xyplot(x ~ y,
                data=do.call(make.groups, lapply(combf, as.data.frame)),
                group=which,
                auto.key=list(space="right")))

   a = do.call(rbind, combf)
   xyplot(x~y|factor(a$channel), a)

   #Features=regFeatures[[i]]
   #Features$reference=refFeatures[[i]][, 1:2]
   #print(xyplot(x ~ y,
   #             data=do.call(make.groups, lapply(Features, as.data.frame)),
   #             pch=0:(length(Features)-1), 
   #   	        group=which, col="black",
   #   	        key=list(space="right", title="samples",
   #                text=list(names(Features)), 
   #                points=list(pch=1:length(Features)-1 )),
   #   	        main=paste("registered features with backgating on", channel)))
}
