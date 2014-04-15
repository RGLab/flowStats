#should be okay as is
#Pulls out the peaks of significant regions from curv1Filter results. Can be used in code that chunks the samples.
#It will be up to the calling function to organize the results for the landmarkMatrix call
getPeakRegions <- function(data, fres, parm, border=0.05, peakNr=NULL, densities=NULL, n=201,
                           indices=FALSE)
{
    ## Some type-checking first
    flowCore:::checkClass(data, "flowSet")
    flowCore:::checkClass(fres, "list")
    flowCore:::checkClass(parm, "character")
    flowCore:::checkClass(border, "numeric", 1)
    ## get peaks of significant regions from curv1Filter results
    ranges <- fsApply(data, range)
    from <- min(sapply(ranges, function(z) z[1,parm]-diff(z[,parm])*0.15), na.rm=TRUE)
    to <- max(sapply(ranges, function(z) z[2,parm]+diff(z[,parm])*0.15), na.rm=TRUE)
    peaks <- list()
    regions <- list()
    eps <- .Machine$double.eps
    for(i in sampleNames(data)){
        if(is.null(densities))
        {
            dens <- densities
            x <- exprs(data[[i]][,parm])
            r <- ranges[[i]][,parm]
            x <- x[x>r[1]+eps & x<r[2]-eps]
        }
        else
        {
            dens <- densities[,i]
            x <- NULL
        }
        tmp <- curvPeaks(fres[[parm]][[i]], x, from=from, to=to, n=n,
                                borderQuant=border, densities=dens)
        peaks[[i]] <- tmp[["peaks"]][,"x"]
        regions[[i]] <- tmp[["regions"]]
    }
    return(list(peaks,regions))
}