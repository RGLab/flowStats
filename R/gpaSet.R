## This function impliments generalized Procrestes analysis (GPA) to align
## high-density regions in a flowSet. The GPA method is concerned with
## multi-dimensional normalization wheares the warping method with
## one-dimentional. 
gpaSet <- function(x, params, register="backgating", bgChannels=NULL,
		   ref.method="median", rotation.only=TRUE,
                   merge.peak.cluster=TRUE,
                   iter.max=10,
                   thres=1e-3, plot=FALSE, ...) 
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
    flowCore:::checkClass(plot, "logical")  ## for debugging purposes
  
    ## 1. identifying and registering (labelling) features
    if (register=="backgating") {
        cat("Backgating ... \n")
        bg <- flowStats:::backGating(x, xy=params, channels=bgChannels)
        features <- flowStats:::idFeatures(bg, merge.peak.clust=TRUE) ## $sample** and $reference
    }
    else { ## use Curve1Filter and landmarkMatrix to find features for each
           ## channels for each flowFrames
        stop("gpaSet: Only Backgating method is available")
    }

    cat("Procrustes analysis ... \n")
    ## 2. translation: translate the centroid of the labelled features of
    ##    each flowFrames to the origin
    translate <- lapply(features, function(x) colMeans(x[, 1:2], na.rm=TRUE))
    
    tfeatures <- list() ## translated registered features
    I <- matrix(1, nrow=nrow(features[[1]]), ncol=1)
    for (i in names(features))
      tfeatures[[i]] <- features[[i]][, 1:2] - I %*% translate[[i]]
 
    ## 3. transformation: applying SVD to find rotation matrix and scalling factor
    
    SVD <- lapply(tfeatures[-which(names(tfeatures)=="reference")],
                  iProcrustes,
                  tfeatures$reference,
                  rotation.only=rotation.only)
    
    ## eliminate boundary points
    for (i in params) {
        bd <- flowCore::boundaryFilter(i)
        x <- flowCore::Subset(x, bd)
    }

    ## define and assign some variables
    expData <- as(x, "list")
    newRange <- matrix(c(Inf, -Inf), nrow=2, ncol=length(params), 
                       dimnames=list(c("min", "max"), params))
  
    ## 4. alignment: translation, rotation, and re-scaling the raw data using
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
    cstr <- function(strs) {
        chstr <- lapply(strs, as.name)
        chstr$sep <- "/"
        chstr <- do.call(paste, chstr, quote=TRUE)
        return(chstr)
    }
    
    chstr <- lapply(params, as.name)
    chstr$sep <- "/"
    chstr <- do.call(paste, chstr, quote=TRUE)
    gpa <- SVD
    attr(gpa, "prior.features") <- features
    attr(gpa, "trans.features") <-
      tfeatures[-which(names(tfeatures)=="reference")]
    attr(gpa, "trans.refFeatures") <- tfeatures$reference
    attr(gpa, "norm.channels") <- cstr(params)
    attr(gpa, "backgating.channels") <- cstr(bgChannels)
    
    class(gpa) <- "GPA"
    ## add backgating channels to attribute

    ## wrap up the result
    regSet <- as(expData, "flowSet")
    phenoData(regSet) <- phenoData(x)
    regSet <- regSet[sampleNames(x)]
    
    if (!is.null(attr(x, "warping")))
        attr(regSet, "warping") <- attr(x, "warping")

    attr(regSet, "GPA") <- gpa

    return(regSet)
}


.checkChannel<- function(ch, allch) {
       mc <- ch %in% allch
       if(!all(mc))
           stop("Invalid parameters not mathcing the flowSet:\n   ",
                 paste(ch[!mc], collapes=", "))
}


