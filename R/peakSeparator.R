## This is the main function to perform backgating. 'data' is supposed to be a flowSet,
## 'xy' is a character vector of length 2 giving the dimensions of interest and 'channels'
## are the additional channels used for backgating. 'fres' is an optional list of filterResults
## of the initial curv1Filter foer each parameter in channels and can be supplied in order to
## speed things up. the output of the function depends on the value of 'incidence'. If this
## is TRUE, it will be a list of incidence matrices for each sample in the set, where each matrix
## indicates whether a particular event was assign to a given peak. If 'incidence' is FALSE, a
## hook function can be supplied which directly uses the subset of events under a peak
## (as a flowFrame in the first argument). The value of xy is passed in as a second argument.
## The output is supposed to be a data frame, potentially of landmark coordinates and further
## annotation. The hook function is called for each frame in the set, for each peak in each
## channel, and the respective results are concatenated using rbind. The peak, channel and
## sample information is added automatically. The final output of the backGating function is
## the concatenated data frame in this case.
backGating <- function(data, xy, channels=setdiff(colnames(data), c(xy, "time", "Time")),
                                                  fres, hookFunction="curvRegs", incidence=FALSE)
{
    if(missing(fres))
    {
        fres <- lapply(channels, function(x) filter(data, curv1Filter(x)))
        names(fres) <- channels
    }
    allres <- NULL
    subsets <- vector(mode="list", length=length(data))
    names(subsets) <- sampleNames(data)
    dataColl <- as(data, "flowFrame")
    for(i in seq_along(channels))
    {
        fresColl <- filter(dataColl, curv1Filter(channels[i]) %subset% sampleFilter(50000))
        np <- length(fresColl)-1
        lm <- landmarkMatrix(data, fres, channels[i], indices=TRUE,
                             peakNr=np)
        if(!all(lm==0))
        {
	    for(j in seq_len(ncol(lm)))
            {
                for(k in sampleNames(data))
                {
                    cat("\r                                                ")
                    cat("\rProcessing channel", channels[i], "peak", j, "sample", k)
                    subset <- split(data[[k]], fres[[channels[i]]][[k]])[[lm[k,j]+1]]
                    if(!is.null(subset))
                    {
                        tfr <-  fres[[channels[i]]][[k]][[lm[k,j]+1]]
                        subsets[[k]] <- cbind(subsets[[k]], tfr@subSet)
                        colnames(subsets[[k]])[ncol(subsets[[k]])] <-
                        paste(channels[[i]], j, sep="_")
                        if(!incidence)
                        {
                            res <- do.call(hookFunction, args=list(subset, xy))
                            if(!is.null(res))
                                allres <- rbind(allres, cbind(res, sample=k, peak=j,
                                                              channel=channels[i]))
                        }
                    }
                }
            }
        }
    }
    cat("\n")
    return(if(incidence) subsets else allres)  
}







## Compute high-density regions using a curv2Filter, try to combine
## neighboring areas if they are within a certain threshold and return
## the polygon vertices as well as theirt centroids.
curvRegs <- function(dat, p)
{
    r <- range(dat)[, p]
    c2 <- curv2Filter(p, bwFac=1.8)
    c2Res <- filter(dat, c2)
    pols <- attr(c2Res@subSet, "polygons")
    polM <- lapply(pols, function(x) sapply(x[2:3], c))
    if(length(polM)>1)
        polM <- pruneRegions(polM, max(apply(r, 2, diff))/5,
                             prod(apply(r, 2, diff))/250)
    res <- NULL
    for(i in seq_along(polM))
    {
        tmp <- polM[[i]]
        colnames(tmp) <- p
        cent <- apply(tmp, 2, mean)
        ntmp <-   data.frame(x=cent[1], y=cent[2], population=i,
                                      type="centroid")
        colnames(ntmp)[seq_along(p)] <- p
        res <- rbind(res,
                     rbind(data.frame(tmp, population=i,
                                      type="polygon"),
                           ntmp))
    }
    return(res)
}






## Some helpers below here...


## Find and score separators between peaks of a flowSet. The scoring
## is based on the overlap of normal approximations for the underlying
## data of neighbouring peaks.

peakSeparator <- function(data, parm, fres)
{
    if(missing(fres))
    {
        fres <- lapply(parm, function(x) filter(data, curv1Filter(x)))
        names(fres) <- parm
    }
    lm <- landmarkMatrix(data, fres, parm, indices=TRUE)
    ## Initially we only take peaks that are present in all samples
    ## lm <- lm[,!apply(lm, 2, function(x) any(is.na(x))), drop=FALSE]
    ##regions <- getRegions(data, parm, fres, lm)
    subs <- getSubsets(data, parm, fres, lm)
    props <- getProportions(data, parm, fres, lm)
    norm <- lapply(subs, normApprox)
    ol <- computeOverlap(norm)
    plotting(data, norm, parm, props, ol)
    
}


getSubsets <- function(data, parm, fres, mat)
{
    res <- list()
    for(i in sampleNames(data))
        res[[i]] <- split(data[[i]][,parm], fres[[parm]][[i]])[mat[i,]+1]
    return(res)
}



getRegions <- function(data, parm, fres, mat)
{
    res <- list()
    for(i in sampleNames(data))
        res[[i]] <- curvPeaks(fres[[parm]][[i]], exprs(data[[i]])[,parm])$regions[mat[i,],]
    return(res)
}


getProportions <- function(data, parm, fres, mat)
{
    res <- list()
    for(i in sampleNames(data))
        res[[i]] <- summary(fres[[parm]][[i]], toTable=TRUE)$p[mat[i,]+1]
    return(res)
}

                      
normApprox <- function(x)
{
    res <- matrix(ncol=2, nrow=length(x))
    colnames(res) <- c("mu", "sd")
    for(i in seq_along(x))
    {
        if(!is.null(x[[i]]))
           res[i,] <- c(median(exprs(x[[i]]), na.rm=TRUE), mad(exprs(x[[i]]), na.rm=TRUE))
        
    }
    return(res)
}



plotting <- function(data, norm, parm, props, ol)
{
    opar <- par(ask=TRUE)
    on.exit(par(opar))
    for(i in sampleNames(data))
    {
        plot(density(exprs(data[[i]][,parm])), main=parm)
        r <- range(exprs(data[[i]]))[,parm]
        for(j in 1:nrow(norm[[i]]))
        {
            sq <- seq(r[1], r[2], len=500)
            tt <- dnorm(sq, mean=norm[[i]][j,1], sd=norm[[i]][j,2])*props[[i]][j]
            lines(sq, tt, col=j+1)
            if(!is.na(ol[[i]][j]))
            {
                xl <- if(nrow(norm[[i]]) >= j+1) mean(c(norm[[i]][j,1], norm[[i]][j+1,1]))
                else norm[[i]][j,1]
                text(xl, par("usr")[4]-(par("usr")[4]/5), signif(ol[[i]][j],2))
            }
        }
    }
    
}



computeOverlap <- function(norm)
{
    res <- list()
    for(n in names(norm))
    {
        nn <- norm[[n]]
        lr <- min(nn[,1] - 10*nn[,2], na.rm=TRUE)
        rr <- max(nn[,1] + 10*nn[,2], na.rm=TRUE)
        tmp <- NULL
        if(nrow(nn)>1)
        {
            for(j in 2:nrow(nn))
            {
                sq <- seq(lr, rr, len=500)
                d1 <-  dnorm(sq, mean=nn[j-1,1], sd=nn[j-1,2])
                d2 <-  dnorm(sq, mean=nn[j,1], sd=nn[j,2])
                tmp <- c(tmp) 
            }
        }
        res[[n]] <- tmp
    }
    return(res)
}





## Compute the area within a polygon by triangular approximation
area.in.polygon <- function(x)
{
    n<-nrow(x)
    D <- diff(x[c(1:n, 1),  ])
    A <- x[, 1] * D[, 2] - x[, 2] * D[, 1]
    abs(sum(A)/2)
}


## A simple Hausdorff distance algorithm for two sets of polygon
## vertices.  The proper solution is to compute distances between
## edges, however this is only defined for convex polygons (which
## could be forced by first calling chull...)
## Note that this makes sense only if the two dimensions are on a
## similar scale, e.g., FSC vs SSC
simpleHd <- function (A, B) 
{
    hd <- 0
    n1 <- nrow(A)
    n2 <- nrow(B)
    d1 <- c()
    d2 <- c()
    if (n1 * n2 == 0) 
        hd <- -1
    else {
        stopifnot(is.numeric(A))
        stopifnot(is.numeric(B))
        for (i in 1:n1) d1[i] <- min(sqrt((B[,1]-A[i,1])^2 + (B[,2]-A[i,2])^2))
        for (i in 1:n2) d2[i] <- min(sqrt((A[,1]-B[i,1])^2 + (A[,2]-B[i,2])^2))
        hd <- max(max(d1), max(d2))
    }
    hd
}


## The pariwise distance matrix for a list of polygons, based on the
## Hausdorff distance metrix.
hdMat <- function(pols)
{
    mat <- matrix(ncol=length(pols), nrow=length(pols))
    for(i in seq_along(pols))
        for(j in seq_along(pols))
            mat[i,j] <- simpleHd(pols[[i]], pols[[j]])
    return(mat)
}


## Collapse two polygons to their common convex hull.
combinePols <- function(p1, p2)
{
    pc <- rbind(p1, p2)
    pc[chull(pc),]
}


## Take a list of polygons and collapse all areas that are closer than
## 'dt' in the Hausdorff distance metric. Also drop all polygons with
## areas smaller than 'st'.
pruneRegions <- function(pols, dt, st)
{
    ssel <- sapply(pols, area.in.polygon) > st
    if(any(ssel))
       pols <- pols[ssel]
    else
       return(pols)
    d <- hdMat(pols)
    sel <- which(hdMat(pols) < dt, arr.ind=TRUE)
    sel <- sel[apply(sel, 1, diff)!=0, ]
    while(length(sel))
    {
        pols[[sel[1,1]]] <- combinePols(pols[[sel[1,1]]],
                                        pols[[sel[1,2]]])
             
        pols <- pols[-sel[1,2]]
        d <- hdMat(pols)
        sel <- which(hdMat(pols) < dt, arr.ind=TRUE)
        sel <- sel[apply(sel, 1, diff)!=0, ]
    }
    return(pols)
}
