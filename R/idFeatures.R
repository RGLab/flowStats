## Internal use only.
## Identify and labeling significant features using backgating return values and
## divisive clustering method (diana from the cluster package).
##
## input argument:
##     bg    A data.frame, must be a return value of the backgating function
##     nDim  The number of dimension of the data subject to normalization
##     thres.merge Threshold (defult to 0.2). A proportion to the high
##                 (distance) between potential features. Used to divide
##                 clusters of features.
##     lambda      Threshold (defult to 0.05). A proportion to the numbe of
##                 samples. Used to decide outliers of features.
##     weight.mfeature Weight bogus features as zero
##     plot.workflow   Plot results of workflow of 'idFeatures'
##
## return value:
##     register        A list containg features for each samples.
##

idFeaturesByBackgating <- function(bg, nDim, thres.sigma=2.5, lambda=0.1,
                                   reference.method="median",
                                   plot.workflow=FALSE,
                                   ask=names(dev.cur())!="pdf")
{   
    ## 'bg' is the return value of function flowStats:::backGating
    if (lambda > 1) {
        message("idFeatures: lambda must be between 0 to 1. Defaulting to
                 0.1. \n")
        lambda <- 0.1
    }
    
    if (! reference.method %in% c("median", "mean"))
        stop(paste(reference.method, "is not a valid method."))
    
    ## 1. take all the centoids from all the backgating and apply 'diana'
    ##    to get the clusters of features
    cent <- bg[bg$type=="centroid", ]
    # cent.translate <- .centerFigure(cent, nDim)
    
    ## cluster centered centroids using 'diana' 
    cent.clust <- .useDivCluster(cent, nDim, thres.sigma, cut.by.height=TRUE)

    ## 2. label the outliers
    filt.clust <- .filterOutlier(cent.clust, lambda=lambda, keep.outlier=TRUE,
                                 label.outlier=".outlier")

    ## 3. get reference features, do not include outliers
    ## reference <- .calcReference(filt.clust, nDim) # if keep.outlier=FALSE
    reference <- .calcReference(
                    filt.clust[!filt.clust$cluster %in% ".outlier", ], nDim,
                    method=reference.method)

    ## 4. register features for each samples (with bogus features)
    perSample <- split(filt.clust, factor(filt.clust$sample))
    register  <- lapply(perSample, .classifyFeatures,
                        refF=reference, nDim=nDim)
    ## 5. prepare return value register (containing the reference features)
   
    register$reference <- reference
    
    if (plot.workflow)
        .plotWorkFlow(cent=cent, cent.clust=cent.clust,
                      filt.clust=filt.clust, register=register, ask=ask)
    
     return(register)
}

###############################
## get the centered figure
###############################
.centerFigure <- function(cent, nDim)
{
  ## cent: centroids of the return values of backgating
  ## return translated centroids (called centered figures)
  TransMatrix <- lapply(split(cent, cent$sample),
                        function(x) mean(x[, 1:nDim], na.rm=TRUE))
  cF <- lapply(split(cent, cent$sample),
               function(x, TransMatrix) {
                   I <- matrix(1, nrow=nrow(x), ncol=1)
                   x[, 1:nDim] <- x[, 1:nDim] - I %*% TransMatrix[[x$sample[1]]]
                   x
               }, TransMatrix)
  
  return(list(centered=do.call(rbind, cF), TransMatrix=TransMatrix))
}

#################################################
## filter out outlier cluster and cluster members
## threshold = lambda * number of total features
#################################################
.filterOutlier <-function(center, lambda=0.1,
                          keep.outlier=FALSE, label.outlier=NULL)
{
    if (keep.outlier & is.null(label.outlier))
        label.outlier <- ".outlier"
    
    nclust <- table(center$cluster)
    outlier_thres <- lambda * sum(nclust)
    #outlier_thres <- median(table(center$cluster)) * lambda
    outlier <- names(nclust)[nclust < outlier_thres]
    if (length(outlier) == length(nclust))
        stop("idFeature: chose smaller lambda.")

    ifelse(keep.outlier,
       center$cluster[center$cluster %in% outlier] <- label.outlier,
       center <- center[!center$cluster %in% outlier, ])

    attr(center, "outlier_thres") <- outlier_thres
    return(center)
}

##############################
## calculate the refenrence ##
##############################
.calcReference <- function(center, nDim=NULL, method="mean")  {
  
    ## return value: ref is a data frame with columns of channels of interest,
    ## cluster, sample, and bogus.
  
    if (is.null(nDim)) stop("Dimension of the variables must be provided")
    
    ref <- do.call(rbind, lapply(split(center, center$cluster),
                   function(a) {
                      ref <- as.data.frame(t(apply(a[, 1:nDim], 2, method)))
                      ref$cluster <- a$cluster[1]
                      return(ref)
                  }))
    ref$sample <- "reference"
    ref$bogus <- rep(FALSE, nrow(ref))
    return(ref)
}


###########################################################
## use divisive hierarchical clustering, diana-algorithm ##
###########################################################
.useDivCluster <- function(center, nDim, thres.sigma=2.5,
                           cut.by.height=TRUE, k=2) {
    ## sign (or re-sign the cluster label)
    clustD  <- cluster::diana(center[, 1:nDim])
    
    ## cut into n groups by 'heightThres' or k
    if (cut.by.height) {
        ## heightThres indicates where the tree should be cut
        heightThres <- sd(clustD$height) * thres.sigma
        center$cluster <- stats::cutree(as.hclust(clustD), h=heightThres)
        attr(center, "height") <- heightThres
        
    }
    else {
        center$cluster <- stats::cutree(as.hclust(clustD), k=k)
        attr(center, "k") <- k
    }

    attr(center, "cluster") <- clustD
    return(center)
}

######################
## classifyFeatures ##
######################
.classifyFeatures <- function(eachSample, refF, nDim)
{
  ## initialize parameter: dim and dim names of reg are equeal xzto that of refF
  reg <- refF
  reg[, 1:nDim] <- NA
  reg$sample <- rep(eachSample$sample[1], nrow(refF))
  
  #reg <- data.frame(apply(refF, 2, function(x) rep(NA, length(x))),
  #                  sample=rep(eachSample$sample[1], nrow(refF)),
  #                  bogus=rep(TRUE, nrow(refF)))
  #reg$cluster <- refF$cluster
  
  perClust <- split(eachSample, factor(eachSample$cluster))

  ## ignore outliers (the label not in refF$clust
  label <- unique(eachSample$cluster)
  anyOutlier <- label[!label %in% refF$cluster]
  if (length(anyOutlier))
      perClust <- perClust[!names(perClust) %in% anyOutlier]


  dim <- seq_len(nDim)
  
  for (j in 1:length(perClust)) {
       lab <- perClust[[j]]$cluster[1]
       if (nrow(perClust[[j]]) == 1) 
            reg[reg$cluster==lab, dim] <- as.matrix(perClust[[j]][, dim])
       else {
           tmp <- rbind(refF[lab, dim], perClust[[j]][,dim])
           
           #we don't need to call compositions::dist here since tmp is data.frame
           # and it always falls back to the stats::dist without any
          # compositions::cdt transformation
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


################################
## plot the results of workflow
################################

.plotWorkFlow <- function(cent, cent.clust, filt.clust, register, ask)
{
  par(ask=ask)
  on.exit(par(ask=FALSE))
  ## 1. cent
  fo <- as.formula(paste(names(cent)[1], "~", names(cent)[2]))
  mainTitle <- do.call(paste, as.list(unique(cent$channel)))
  mainTitle <- paste("Centroids of high-density areas yielded by
                     backgating on channel", mainTitle)
  print(lattice::xyplot(fo, group=channel,
                        data=cent,
                        auto.key=list(cex=0.8, title="channel", space="right"),
                        main=mainTitle,
                        sub="potential features"))
  
  col <- rainbow(length(unique(cent$sample)))
  key <- list(text=list(levels(factor(cent$sample))),
              points=list(col=col, pch=20, alpha=0.7),
              cex=0.7, title="sample")
  
  print(lattice::xyplot(fo, group=sample,
                        data=cent,
                        col=col, pch=20, alpha=0.7,
                        key=c(key, space="right"),
                        main=mainTitle,
                        sub="potential features"))
  
  ## 2. cent.translate$centered
  #print(lattice::xyplot(fo, group=sample,
  #                     data=cent.translate$centered,
  #                     auto.key=list(cex=0.8, title="sample", space="right"),
  #                     main="Centered figure: translate the center of
  #                           centroids to the origin"))
  
  ## 2. cent.clust
  subTitle <- paste("Clusters containing <",
                    format(attr(cent.clust, "height"), digit=3),
                    "members are considered outliers.")
  
  
  print(plot(attr(cent.clust, "cluster"), which.plot=2,
             xlab="Potential Features",
             main=paste("Cut the tree at height",
                         format(attr(cent.clust, "height"), digit=3))
             ))
  
  d <- density(attr(cent.clust, "cluster")$height)
  print(plot(d, main="Distribution of the cluster (diana) height"))
  abline(v=attr(cent.clust, "height"), col="grey")
  text(x= attr(cent.clust, "height"), y=mean(d$y)+2*sd(d$y),
       labels="cutting the\n tree here", pos=4)
  
  print(lattice::xyplot(fo, group=cluster,
                   data=cent.clust,
                   auto.key=list(cex=0.8,title="cluster", space="right"),
                   main="Clusters of potential features",
                   sub=subTitle))
  ## plot dendrogram of attr(cent.clust, "cluster")
  
  
  ## 3. filt.clust
  print(lattice::xyplot(fo, group=cluster,
                        data=filt.clust,
                        auto.key=list(cex=0.8, title="channel", space="right"),
                        sub=subTitle,
                        main="Clusters of features and outliers"))
                        
  ## 4. register -- final features and reference features
  fo <- as.formula(paste(names(cent)[1], "~", names(cent)[2],
                         " | sample"))
  print(lattice::xyplot(fo , group=cluster,
                        data=do.call(rbind, register),
                        auto.key=list(cex=0.8, title="label",space="right"),
                        main="Labeled features"))
  print(lattice::xyplot(fo, group=bogus,
                        data=do.call(rbind, register),
                        auto.key=list(cex=0.8, title="Bogus",space="right"),
                        main="Labeled features",
                        sub="Missing features are those with bogus=TURE."))
        
}
