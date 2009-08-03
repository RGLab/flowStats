## use backgeting to find the reference features and register the features
## used by function gpaSet

idFeatures <- function(bg, plot=FALSE)
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
        ## cluster the centriod based on peak+population
        center  <- labelCentroid(center)
        ## let reference features = median of cluster members
        refFeatures[[i]] <- t(sapply(split(center, center$cluster),
                                function(a) apply(a[, 1:2], 2, median)))
                  
        perSample = split(center, factor(center$sample))
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
         ## rownames(combf[[m]]) <- as.character(1:nrow(combf[[m]]))
    }
    
    reference <- do.call(rbind, refFeatures)
    rownames(reference) <- rownames(combf[[1]])
    reference <- as.data.frame(reference)
    reference$bogus <- rep(FALSE, nrow(reference))
    reference$cluster <- rownames(combf[[1]])
 
    
    combf <- lapply(combf, function(x)
                    {x$channel=rownames(x); x})
    combf[["reference"]] <- reference
    #pts <- do.call(rbind, combf)
    #reference <- t(sapply(split(do.call(rbind, combf), factor(pts$channel)),
    #                      function(x) sapply(x[, 1:2], median)))
    #reference <- reference[rownames(combf[[1]]$channel,]
    ## this reference = do.cal(rbind, refFeatures)
    #if (plot) .plotFeatures(combf, reference)
    
    return(combf)
}

###################
## classifyFeatures
###################
classifyFeatures <- function(eachSample, refF, by="distance")
{
  ## make the size of refF and column names consistant with refF (later)
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
        
        ##plot the result for peak 1
        ##peaks <- perPeak[[i]]
        ##peaks$label <- km$cluster
        ##x11();print(xyplot(x~y, data=peaks, group=peaks$label))
    }

    return(do.call(rbind, perPeak))
    
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
