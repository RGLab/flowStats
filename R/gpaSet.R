## This function impliments generalized Procrestes analysis (GPA) to align
## high-density regions in a flowSet. The GPA method is concerned with
## multi-dimensional normalization wheares the warping method with
## one-dimentional. 
gpaSet <- function(x, params, 
       	  	   register="backgating", bgChannels=NULL,
		   ref.method="mean",
                   rotation.only=TRUE, iter.max=10, thres=1e-3, fig=FALSE)
{
    ## check valid arguments
    flowCore:::checkClass(x, "flowSet")

    flowCore:::checkClass(params, "character")
    .checkChannel(params, colnames(x))
    if (length(params) < 2) {
        stop("At least two params are required to apply multi-dimensional
              normalization \n") }

    flowCore:::checkClass(register, "character")
    ## what if register=="curv1Filter"
    if (is.null(bgChannels)) {
       cat("Argument 'bgChannels' is missing. Default with channels other than",
           params, ".\n")
       bgChannels <- setdiff(colnames(x), c(params, "time", "Time"))
    }
    else {
        flowCore:::checkClass(bgChannels, "character")
        .checkChannel(bgChannels, colnames(x))
    }
 
    flowCore:::checkClass(rotation.only, "logical")
    flowCore:::checkClass(ref.method, "character")
    flowCore:::checkClass(iter.max, "numeric")
    flowCore:::checkClass(thres, "numeric")
    flowCore:::checkClass(fig, "logical")  ## for debugging purposes

    ## 1. identifying and registering (labelling) features
    if (register=="backgating") {
        cat("Backgating ... \n")
        bg <- backGating(x, xy=params, channels=bgChannels)
        regFeatures <- useBackGating(bg, xy = params, fig=fig) 
    }
    else { ## use Curve1Filter and landmarkMatrix to find features for each
           ## channels for each flowFrames 
    }
    
    cat("Procrustes analysis ... \n")
    ## 2. translating the centroid of the registered features of each flowFrames
    ## to the origin
    translate <- lapply(regFeatures, colMeans, na.rm=TRUE)
    tfeatures <- list() ## translated registered features, debugging
    I <- matrix(1, nrow=nrow(regFeatures[[1]]), ncol=1)
    for (i in names(regFeatures))
      tfeatures[[i]] <- regFeatures[[i]] - I %*% translate[[i]]

    ## find the reference feature to which the translated features to be aligned
    fbar <- .getRefFeatures(tfeatures, ref.method=ref.method)
 
    ## 3. applying SVD to find rotation matrix and scalling factor
    SVD <- lapply(tfeatures, iProcrustes, fbar, rotation.only=rotation.only)
    
    ## eliminate boundary points
    for (i in params) {
        bd <- flowCore::boundaryFilter(i)
        x <- flowCore::Subset(x, bd)
    }

    ## define and assign some variables
    expData <- as(x, "list")
    newRange <- matrix(c(Inf, -Inf), nrow=2, ncol=length(params), 
                       dimnames=list(c("min", "max"), params))
  
    ## 4. Alignment: translation, rotation, and re-scaling the raw data using
    ##    GPA 
    for (i in names(expData))
    {
        newSet <- exprs(expData[[i]])[, params]
        newSet <- t(t(newSet) - translate[[i]]) 
        newSet <- SVD[[i]]$scal * (newSet %*% SVD[[i]]$Q)
        exprs(expData[[i]])[, params] <- newSet
        ## ranges for the translated data across samples
        rng <- apply(newSet, 2, range)
        dimnames(rng) <- dimnames(newRange)
        newRange <- rbind(min=pmin(newRange["min", ], rng["min", ]),
	                  max=pmax(newRange["max", ], rng["max", ]),
			  deparse.level=2)
    }

    ## update parameters slot: minRange and maxRange
    for (i in sampleNames(x))
    {
        ip <- match(params, pData(parameters(expData[[i]]))$name)
        tmp <- parameters(expData[[i]])
	oldrng <- as.matrix(range(expData[[i]], params))
	pData(tmp)[ip, c("minRange", "maxRange")] <- t(newRange)
        expData[[i]]@parameters <- tmp
    }
    
    ## set up flowSets GPA attribute
    chstr <- lapply(params, as.name)
    chstr$sep <- "/"
    chstr <- do.call(paste, chstr, quote=TRUE)
    gpa <- SVD
    attr(gpa, "trans.feature") <- tfeatures
    attr(gpa, "prior") <- regFeatures
    attr(gpa, "channels") <- chstr
    attr(gpa, "refFeature") <- fbar
    class(gpa) <- "GPA"
    ## add backgating channels to attribute

    ## wrap up the result
    regSet <- as(expData, "flowSet")
    phenoData(regSet) <- phenoData(x)
    regSet <- regSet[sampleNames(x)]
    attr(regSet, "GPA") <- gpa

    return(regSet)
}

## get reference features (multi-dimensional)
.getRefFeatures <- function(features, ref.method="mean")
{
    
    fbar <- array(0, dim=dim(features[[1]]), 
                     dimnames=dimnames(features[[1]]))
    if (ref.method != "mean")
       cat("Only mean-value method is avaliable.\n")   

    for (i in rownames(features[[1]]))
         fbar[i, ] <- 
            rowMeans(sapply(features, function(x) x[i, ]), na.rm=TRUE)

    return(fbar)
}

.checkChannel<- function(ch, allch) {
       mc <- ch %in% allch
       if(!all(mc))
           stop("Invalid parameters not mathcing the flowSet:\n   ",
                 paste(ch[!mc], collapes=", "))
}


