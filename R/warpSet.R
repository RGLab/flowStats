## Align data in a flowSet by estimating high density regions and using this
## information as landmarks. This works separately on each parameter.
warpSet <- function(x, stains, grouping=NULL, monwrd=TRUE, subsample=NULL,
                    peakNr=NULL, clipRange=0.01, nbreaks=11, fres, warpFuns=FALSE,
                    ...)
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
    wfuns <- list()
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
        regions <- attr(landmarks, "regions")
        dists <- attr(landmarks, "cdists") 
        attr(landmarks, c("regions", "cdists")) <- NULL
        hasPeak <- !is.na(landmarks)
        for(n in 1:ncol(landmarks))
        {
            nar <- is.na(landmarks[,n])
            landmarks[nar,n] <- mean(landmarks[,n], na.rm=TRUE)
            m <- matrix(apply(sapply(regions[!nar],
                                     function(r) if(is.null(dim(r))) r else r[n,]),
                              1, mean, na.rm=TRUE), ncol=2)
            for(i in names(which(nar)))
                regions[[i]][n,] <- m
        }
            
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
        names(funs) <- sampleNames(x)
        wfuns[[p]] <- funs
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
                            warpFun=funs, regions=regions, dists=dists,
                            hasPeak=hasPeak)
        }
    }
    if(warpFuns)
        return(wfuns)
    regSet <- as(expData, "flowSet")
    phenoData(regSet) <- phenoData(x)
    regSet <- regSet[sampleNames(x)]
    attr(regSet, "warping") <- lm
    regSet
}





## Some QA pots for a normalized data set
normQA <- function(data, morph=c("^fsc", "^ssc"), channels, ask=TRUE)
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

    ## Plot the amount of warping for each landmark
    layout(matrix(1:length(wchans), ncol=1))
    for(p in wchans)
    {
        nc <- ncol(ninfo[[p]]$warped)
        nr <- nrow(ninfo[[p]]$warped)
        matplot(ninfo[[p]]$prior, type="n",
                ylab=p, xlab="sample", ylim=range(data[[1]])[,p],
                main=paste("Landmark Warping", p))
        pch <- ifelse(ninfo[[p]]$hasPeak, 20, 4)
        for(i in seq_len(nc))
            points(ninfo[[p]]$prior[,i], pch=pch[,i], col=i+1)
        m <- apply(ninfo[[p]]$warped, 2, mean, na.rm=TRUE)
        abline(h=m, col=seq_len(nc)+1, lty="dotted")
        diffs <- abs(ninfo[[p]]$prior - ninfo[[p]]$warped)
        for(i in seq_len(nc))
          segments(x0=seq_len(nr), x1=seq_len(nr), y0=ninfo[[p]]$prior[,i],
                   y1=m[i], col=i+1, lty="dotted")
    }
    par(ask=ask)
    on.exit(par(ask=FALSE))

    ## Plot the confidence of the landmark registration step
    layout(matrix(1:length(wchans), ncol=1))
    for(p in wchans)
    {
        nc <- ncol(ninfo[[p]]$warped)
        nr <- nrow(ninfo[[p]]$warped)
        dists <- 1-ninfo[[p]]$dists
        #dists[is.na(dists)] <- 1
        pch <- ifelse(ninfo[[p]]$hasPeak, 1, 4)
        plot(rep(seq_len(nc), each=nr), as.vector(dists),
             ylim=c(min(0.5, min(dists, na.rm=TRUE)), 1.05),
             xlim=c(0, nc+1), xlab="peaks", ylab="confidence", xaxt="n",
             pch=pch,
             main=paste("Landmark Registration", p))
        axis(1, at=seq_len(nc), labels=TRUE)
        box()
        present <- apply(ninfo[[p]]$dists, 2, function(x) sum(!is.na(x)))
        text(seq_len(nc), rep(1.05, nc), paste(present, nr, sep="/"),
             cex=0.8)
        
    }

    ## Plot the warping functions
    layout(matrix(1:length(wchans), ncol=1))
    for(p in wchans)
    {
       xvals <-  seq(range(data[[1]], p)[1,], range(data[[1]], p)[2,], len=50)
       yvals <- sapply(ninfo[[p]]$warpFun, function(fun) fun(xvals))
       matplot(yvals, xvals, type="l", xlab="warped", ylab="original",
               main=paste("Warping Functions", p))
    }

    ## Plot the projection in the FCS vs. SSC space
    alldat <- Subset(as(data, "flowFrame")[,c(fsc, ssc)],
                     boundaryFilter(c(fsc, ssc)))
    allM <- allV <- list()
    layout(matrix(1:length(wchans), ncol=1))
    for(p in wchans)
    {
        nc <- ncol(ninfo[[p]]$warped)
        nr <- nrow(ninfo[[p]]$warped)
        ## First create rectangleGates from the peak regions
        regions <- ninfo[[p]]$regions
        wfuns <- ninfo[[p]]$warpFun
        datc <-  Subset(data, boundaryFilter(c(fsc, ssc)))
        nfc <- filter(datc, norm2Filter(c(fsc, ssc)))
        vcc <- sapply(nfc, function(x)
                         prod(sqrt(eigen(filterDetails(x)[[1]]$cov)$values))*pi)
        contour(alldat, xlab=fsc, ylab=ssc,
                main=paste("Backgating Shape", p), col="transparent")
        tmp <- sapply(regions, rowMeans)
        expval <- signif(if(!is.null(dim(tmp))) rowMeans(tmp) else mean(tmp),3)
        legend("topleft", pch=20, col=seq_len(nc)+1,
               legend=sprintf("Peak %d (mean %s=%g)", seq_len(nc), p, expval),
               cex=0.8, bty="n")
        for(i in seq_len(nc))
        {
            g <- list()
            for(j in names(regions))
            {
                g[[j]] <- rectangleGate(matrix(wfuns[[j]](regions[[j]][i,]), ncol=1,
                                               dimnames=list(NULL, p)))
            }
            dat <- Subset(Subset(data, filterList(g)), boundaryFilter(c(fsc, ssc)))
      
            crgb <- col2rgb(i+1)
            col <- rgb(crgb[1,], crgb[2,], crgb[3,], 255*0.12, maxColorValue=255)
            colc <- rgb(crgb[1,], crgb[2,], crgb[3,], 255*0.2, maxColorValue=255)
            nf <- filter(dat, norm2Filter(c(fsc, ssc)))
            for(f in nf)
                gpolygon(f, col=col, border=colc)
            m <- fsApply(dat, function(x) apply(x[,c(fsc, ssc)], 2, mean, na.rm=TRUE),
                         use.exprs=TRUE)
            v <- fsApply(dat, function(x) apply(x[,c(fsc, ssc)], 2, sd, na.rm=TRUE),
                         use.exprs=TRUE)
            vc <- sapply(nf, function(x)
                         prod(sqrt(eigen(filterDetails(x)[[1]]$cov)$values))*pi)
            #segments(x0=m[,1]-v[,1], x1=m[,1]+v[,1], y0=m[,2], y1=m[,2], col=i,
            #         lwd=1)
            #segments(x0=m[,1], x1=m[,1], y0=m[,2]-v[,2], y1=m[,2]+v[,2], col=i,
            #         lwd=1)
            #points(m, col=i)
            relV <- vc/vcc
            allV[[p]][[i]] <- relV
            allM[[p]][[i]] <- m
        }
    }

    for(p in wchans)
    {
        nc <- ncol(ninfo[[p]]$warped)
        nr <- nrow(ninfo[[p]]$warped)
        contour(alldat, col="lightgray", xlab=fsc, ylab=ssc,
                main=paste("Backgating Centroids", p))
        regions <- ninfo[[p]]$regions
        tmp <- sapply(regions, rowMeans)
        expval <- signif(if(!is.null(dim(tmp))) rowMeans(tmp) else mean(tmp),3)
        legend("topleft", pch=20, col=seq_len(nc)+1,
               legend=sprintf("Peak %d (mean %s=%g)", seq_len(nc), p, expval),
               cex=0.8, bty="n")
        for(i in seq_len(nc))
            points(allM[[p]][[i]], col=i+1)
    }

    ## Plot the relative ellipse volumes for the backgated channels
    layout(matrix(1:length(wchans), ncol=1))
    for(p in wchans)
    {
        nc <- ncol(ninfo[[p]]$warped)
        nr <- nrow(ninfo[[p]]$warped)
        plot(seq_len(nc), rep(0, nc),
             ylim=c(0,max(5, max(unlist(allV[[p]])))), type="n",
             xlim=c(0, nc+1), xlab="peaks", ylab="variation", xaxt="n",
             pch=pch, col=rep(seq_len(nc), each=nr),
             main=paste("Backgating Variation", p))
        axis(1, at=seq_len(nc), labels=TRUE)
        box()
        for(i in seq_len(nc))
            points(rep(i, nr), allV[[p]][[i]]) 
    }
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
