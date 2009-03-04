## parse output of curv1Filter and find modes and midpoints of the density
## regions
curvPeaks <- function(x, dat, borderQuant=0.01, n=201, from, to)
{
    ## Some type-checking first
    flowCore:::checkClass(x, "multipleFilterResult")
    flowCore:::checkClass(dat, "numeric")
    flowCore:::checkClass(borderQuant, "numeric", 1)
    flowCore:::checkClass(n, "numeric", 1)
    if(missing(from))
        from <- min(dat)
    flowCore:::checkClass(from, "numeric", 1)
    if(missing(to))
        to <- max(dat)
    flowCore:::checkClass(to, "numeric", 1)
    ## extract boundaries
    bound <- attr(x@subSet, "boundaries")
    dens <- density(dat, n=n, from=from, to=to)
    ## iterate over regions
    regPoints <- list()
    peaks <- midpoints <- regions <- densFuns <- NULL
    i <- 1
    if(!all(is.na(bound[[1]]))){
        oo <- options(warn=-1)
        on.exit(options(oo))
        for(b in bound){
            ## discard regions on the margins
            if(b[2] > quantile(c(from,to), borderQuant) &&
               b[1] < quantile(c(from,to), 1-borderQuant)){
                ## approximate density by function
                afun <- approxfun(dens)
                sel <- seq(b[1], b[2], len=50)
                regPoints[[i]] <- cbind(x=sel, y=afun(sel))
                ## compute maximum of function
                m <- optimize(afun, b, maximum=TRUE)
                peaks <- rbind(peaks, cbind(x=m$maximum, y= m$objective))
                regions <- rbind(regions, cbind(left=b[1], right=b[2]))
                midpoints <- c(midpoints, mean(b))
                densFuns <- c(densFuns, afun)
                i <- i+1
            }
        }
    }
    if(i==1)
        return(list(peaks=cbind(x=NA, y=NA), regions=cbind(left=NA, right=NA),
                    midpoints=NA, regPoints=list(cbind(x=NA, y=NA)),
                    densFuns=NA))
    return(list(peaks=peaks, regions=regions, midpoints=midpoints,
                regPoints=regPoints, densFuns=densFuns))
}

