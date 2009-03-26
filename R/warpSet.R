## Align data in a flowSet by estimating high density regions and using this
## information as landmarks. This works separately on each parameter.
warpSet <- function(x, stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
                    peakNr=NULL, clipRange=0.01, nbreaks=11, fres, ...)
{
    ## Some type-checking first
    flowCore:::checkClass(x, "flowSet")
    flowCore:::checkClass(stains, "character")
    mt <- stains %in% colnames(x)
    if(!all(mt))
        stop("Invalid stain(s) not matching the flowSet:\n    ",
             paste(stains[!mt], collapse=", "))
    expData <- as(x, "list")
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
     
    ## find landmarks
    if(missing(fres))
    {
        fres <- list()
        for(p in stains){
            cat("\rEstimating landmarks for channel", p, "...")
            fres[[p]] <- filter(x, curv1Filter(p, bwFac=1.3))
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
        fdobj   <- data2fd(densY, argvals, wbasis, 
                           argnames = c("x", "samples", "density"))

        ## create matrix of landmarks from curv1Filter peaks
        cat("Registering curves for parameter", p, "...\n")
        landmarks <- landmarkMatrix(x, fres, p, border=clipRange, peakNr=peakNr,
                                    densities=densY, n=nb)
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
        for(n in 1:ncol(landmarks))
            landmarks[is.na(landmarks[,n]),n] <-
                mean(landmarks[,n], na.rm=TRUE)
        ## register the densities
        if(ncol(landmarks)==1){ ## only one peak: offset
            offsets <- landmarks-median(landmarks)
            funs <- vector("list", length(landmarks))
            for(j in seq_along(funs)){
                funs[[j]] <- function(x) x - z
                e1 <- new.env()
                e1$z <- offsets[j]
                environment(funs[[j]]) <- e1
            }
        }else{ ## multiple peaks: warping
            capture.output(regDens <- landmarkreg(fdobj, landmarks, WfdPar=WfdPar, 
                                   monwrd=monwrd, ...))  
            warpfdobj <- regDens$warpfd
            warpedX <- eval.fd(warpfdobj, argvals)
            ## compute warping functions
            ## funs <-  apply(warpedX, 2, function(y) approxfun(argvals, y))
            funs <-  apply(warpedX, 2, approxfun, argvals)
        }
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
            sel <- leftBoard[[i]] | rightBoard[[i]]
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
                        warpFun=funs)
    }
    regSet <- as(expData, "flowSet")
    phenoData(regSet) <- phenoData(x)
    regSet <- regSet[sampleNames(x)]
    attr(regSet, "warping") <- lm
    regSet
}

