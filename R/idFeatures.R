## identify and labeling significant features using backgating return values and
## divisive clustering method -- diana.
idFeatures <- function(bg, nDim, thres.merge=0.2, weight.mfeature=FALSE,
                           plot=FALSE)
{
    ## 'bg' is the return value of function flowStats:::backGating

    if (thres.merge > 1) {
        thres.merge <- 0.2
        message("idFeatures: thres.merge must between 0 to 1. Defaulting
        to 0.2.\n")
    }
    
    regFeatures <- list()
    refFeatures <- list()
    ch <- split(bg, factor(bg$channel))
    lambda <- 0.1
    
    ## identify features on the projection of each bg channels
    for (i in names(ch))
    { 
        perChannel <- ch[[i]]
        center <- perChannel[perChannel$type=="centroid", ]
        ## use diversive hierarchical clustering (diff from agglomerative)
        center<-.useDivCluster(center, 2, thres.merge, h.cluster="diana") 

        ## flag the outliers
        center <- .filterOutlier(center, lambda, channel=i, keep.outlier=TRUE)
        
        refFeatures[[i]] <- .calcReference(center)
                                                       
        perSample = split(center, factor(center$sample))
        
        ## if cluster labeled as an outlier...
        idx <- grep("outlier", refFeatures[[i]]$cluster)
        if (length(idx) > 0)
            refF <- refFeatures[[i]][-idx, ]
        else
            refF <- refFeatures[[i]]

        ## classify features for each sample
        regFeatures[[i]] <- 
            lapply(perSample, .classifyFeatures,
                   refF=refF, nDim=nDim)

    }

    ## combine all registered features and group features for each sample
    combf <- list()
    for (m in levels(factor(bg$sample))) {
         combf[[m]] <- do.call(rbind, 
                     lapply(regFeatures, function(a, b) {a[[b]]}, m))
    }
    ## plot the result
    if (plot)
        print(lattice::xyplot(FSC~SSC|sample, group=cluster,
                              data=do.call(rbind,combf),
                              main="feature labeling",
                              auto.key=list(space="right")
                              ))
    ## combind all the 
    combf <- .useDivCluster(do.call(rbind, combf),
                            nDim=2, h.cluster="diana",
                            thres.merge=thres.merge,
                            naming.channel=FALSE)
    ## plot the result
    if (plot)
        print(lattice::xyplot(FSC~SSC|sample, group=cluster,
                              data=combf,
                              main="feature labeling: using diana",
                              auto.key=list(space="right")
                              ))           
    
    ## reduce outlier
    combf <- .filterOutlier(combf, lambda=lambda, keep.outlier=FALSE)
    reference <- .calcReference(combf)
    reference$bogus <- rep(FALSE, nrow(reference))

    
    ## need to classify the cluster center again, pick the feature that is
    ## closer to the reference
    combf <-lapply(split(combf, factor(combf$sample)),
              function(a, refF, nDim) {
                 perCluster <- split(a, factor(a$cluster))
                 reg <- lapply(perCluster, function(b) {
                     ## take false if there is any
                     if (sum(!b$bogus)>0)
                         b <- b[!b$bogus, ]
                     
                     tmp <- rbind(refF[b$cluster[1], 1:nDim], b[, 1:nDim])
                     minD <- b[which.min(as.numeric(dist(tmp))[1:nrow(b)]), ]
                 })
                 do.call(rbind, reg)
              }, reference, nDim)
    ## check each sample has the same amont of feature
    ## combf <- .checkFeatureNum(combf, reference, nDim)
    if (plot)
        print(lattice::xyplot(FSC ~ SSC | sample, group=cluster,
                        data=do.call(rbind, combf),
                        auto.key=list(title="cluster",space="right"),
                        main="reduced features"))
        #print(lattice::xyplot(FSC ~ SSC | sample, group=bogus,
        #                data=do.call(rbind, combf),
        #                auto.key=list(title="Bogus",space="right"),
        #                main="reduced features"))
    combf[["reference"]] <- reference

     return(combf)
}

#################################
## filter out outlier. 
#################################
.filterOutlier <-function(center, lambda, channel=NULL, keep.outlier=FALSE) {
    nclust <- sapply(split(center, factor(center$cluster)),
                         function(x) nrow(x))
    outlier_thres <- round(length(unique(center$sample)) * lambda)
    outlier <- names(nclust)[nclust <= outlier_thres]

    if (length(outlier) > 0) { ## if there are outliters
        if (keep.outlier) {
            center$cluster[center$cluster %in% outlier] <-
                 paste(channel, "outlier", sep=".")
        }
        else
            center <- center[!center$cluster %in% outlier, ]
     }
    return(center)
}

###########################
## calculate the refenrence
###########################
.calcReference <- function(center, method="median")  {
      ref <- do.call(rbind, lapply(split(center, center$cluster),
                 function(a) {
                    ref <- as.data.frame(t(apply(a[, 1:2], 2, median)))
                    ref$cluster <- a$cluster[1]
                    return(ref)
                  }))
}


#################################################
## use divisive clustering method: diana or diasy
#################################################
.useDivCluster <- function(center, nDim, naming.channel=TRUE,
                           thres.merge=0.2, h.cluster="diana") {
    ## sign (or re-sign the cluster label)
    clustD  <- cluster::diana(center[, 1:nDim])
    heightThres <- max(clustD$height)*thres.merge

    if (naming.channel)
        center$cluster <- paste(center$channel, 
                                stats::cutree(as.hclust(clustD),
                                              h=heightThres),
                            sep=".")
    else
        center$cluster <- stats::cutree(as.hclust(clustD), h=heightThres)
    
    return(center)
}

###################
## classifyFeatures
###################
.classifyFeatures <- function(eachSample, refF, nDim)
{
  ## make the size of refF and column names consistant with refF 
  ## initialize parameter

  reg <- data.frame(apply(refF, 2, function(x) rep(NA, length(x))),
                    sample=rep(eachSample$sample[1], nrow(refF)),
                    bogus=rep(TRUE, nrow(refF)))
  reg$cluster <- refF$cluster
  
  perClust <- split(eachSample, factor(eachSample$cluster))
  if (!identical(grep("outlier", names(perClust)), integer(0)))
      perClust <- perClust[-grep("outlier", names(perClust))]

  dim <- seq_len(nDim)
  
  for (j in 1:length(perClust)) {
       lab <- perClust[[j]]$cluster[1]
       if (nrow(perClust[[j]]) == 1) 
            reg[reg$cluster==lab, dim] <- as.matrix(perClust[[j]][, dim])
       else {
           tmp <- rbind(refF[lab, dim], perClust[[j]][,dim])
           mindist <-
                   which.min(as.numeric(dist(tmp))[1:nrow(perClust[[j]])])
           reg[reg$cluster==lab, dim] <- as.matrix(perClust[[j]][mindist, dim])
       }
       reg$bogus[reg$cluster==lab] <- FALSE
   }

  ## 2. deal with NA peaks, give the refFeatures and cluster status as bogus
    naIdx <- which(is.na(reg[, 1]))
  
    if (!identical(naIdx, integer(0))) {
       reg[naIdx, dim] <- refF[naIdx, dim]
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
    	                
   #print(xyplot(x ~ y|which,
   #             do.call(make.groups, lapply(regFeatures[[i]], as.data.frame)),
   #             as.table=TRUE, 
   #             main=paste("registered features with backgating on", channel)))
            
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


###################################################
## check if the number of features for
## each sample is equal to that of reference
## if the feature is missing, replace with a bogus
###################################################
.checkFeatureNum <- function(combf, reference, nDim) {
    lapply(combf, function(x) {
        if (nrow(x) != nrow(reference)) {
           i <- nrow(x)
           x[i+1,] <- x[i,]
           x[i+1, 1:(nDim+1)] <-
               reference[!(reference$cluster %in% x$cluster), 1:(nDim+1)]
           x$bogus[i+1] = TRUE
           x <- x[sort(x$cluster), ]
        }
        x
    })
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
