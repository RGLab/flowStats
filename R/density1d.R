## Find most likely separator between peaks in 1D
density1d <- function(x, stain, alpha="min", sd=2, plot=FALSE, borderQuant=0.1,
                      absolute=TRUE, inBetween=FALSE, refLine=NULL, ...)
{
    ## some type checking first
    flowCore:::checkClass(x, c("flowFrame", "flowSet"))
    flowCore:::checkClass(stain, "character", 1)
    if(!stain %in% colnames(x))
        stop("'", stain,"' is not a valid parameter in this flowFrame")
    flowCore:::checkClass(alpha, c("character", "numeric"), 1)
    flowCore:::checkClass(sd, "numeric", 1)
    flowCore:::checkClass(plot, "logical", 1)
    flowCore:::checkClass(borderQuant, "numeric", 1)
    flowCore:::checkClass(absolute, "logical", 1)
    if (!is.null(refLine))
        flowCore:::checkClass(refLine, "numeric", 1)
    flowCore:::checkClass(inBetween, "logical", 1)
    
    ## collapse to flowFrame
    if(is(x, "flowSet"))
       x <- as(x, "flowFrame")
    ## clip to data range (if needed) and compute density peaks
    if(absolute){
        vrange <- range(x, na.rm=TRUE)[, stain]
        vrange[is.nan(vrange)] <- 0
        vrange[1] <- min(vrange[1], min(exprs(x)[,stain], na.rm=TRUE))
    }
    else  
        vrange <- range(exprs(x)[,stain], na.rm=TRUE)

    vrange[1] <- ifelse(is.null(refLine), vrange[1], max(vrange[1], refLine))
    inc <- diff(vrange)/1e5
    exprf <- char2ExpressionFilter(sprintf("`%s` > %s & `%s` < %s", stain,
                                           vrange[1]+inc, stain, vrange[2]-inc))
    fres <- filter(tmp <- Subset(x, exprf), curv1Filter(stain))
    bnds <- curvPeaks(fres, exprs(tmp)[, stain], borderQuant=borderQuant)
    ## define 'positive' and 'negative' peaks by proximity to anchors
    anchors <- quantile(as.numeric(vrange), c(0.25, 0.5))
    anc <- abs(sapply(anchors, "-", bnds$peaks[, "x", drop=FALSE]))
    dens <- density(exprs(tmp)[, stain])
    ## only one peak: use robust estimation of mode and variance for boundary
    est <- hubers(exprs(tmp[, stain]))
    if(is.null(nrow(anc)))
    {
        class <- which.min(abs(est$mu - anchors))
    }
    else if(nrow(anc)==2 && inBetween)
    {
        left <- bnds$regions[1,"right"]
        right <-  bnds$regions[2,"left"]
        class <- 1:2
    }
    else
    {
        ## multiple peaks: classify as pos or neg and set boundary between
        class <- apply(anc,1, which.min)
        left <- ifelse(sum(class==1)>0, max(bnds$regions[class==1, "right"]),
                       min(bnds$regions[class==2, "left"]))
        right <- ifelse(sum(class==2)>0, min(bnds$regions[class==2, "left"]),
                        max(bnds$regions[class==1, "right"]))
    }
    if(all(class==1))
        bound <- est$mu + sd * est$s
    else if(all(class==2))
        bound <- est$mu - sd * est$s
    else{ 
        if(is.character(alpha)){
            if(alpha=="min"){
                sel <- (dens$x > left & dens$x <= right)
                bound <- dens$x[sel][which.min(dens$y[sel])]
            }else
            stop("unknown value for alpha")
        }else
        bound <- (1-alpha)*left + alpha*right    
    }
    ## create output if needed
    if(plot){
        plot(dens, main=paste("breakpoint for parameter", stain), cex.main=1, ...)
        regs <- bnds$regions
        for(i in 1:nrow(regs)){
            sel <- dens$x >= regs[i,1] & dens$x <= regs[i,2]
            polygon(c(dens$x[min(which(sel))], dens$x[sel],
                      dens$x[max(which(sel))]), c(0, dens$y[sel], 0),
                    col="lightgray", border=NA)
        }
        lines(dens, ...)
        abline(v = bound, col = 2, lwd = 2)
        legend("topright", c("breakpoint", "dens region"),
               fill=c("red", "lightgray"), bty="n")
    }
    return(bound)
}


## A wrapper around density1D directly creating a range gate.
rangeGate <- function(x, stain, alpha="min", sd=2, plot=FALSE, borderQuant=0.1,
                     absolute=TRUE, filterId="defaultRectangleGate",
                     positive=TRUE, refLine=NULL, ...)
{
    loc <- density1d(x=x, stain=stain, alpha=alpha, sd=sd, plot=plot,
                     borderQuant=borderQuant, absolute=absolute, refLine=refLine, ...) 
    bounds <- if(positive) list(c(loc, Inf)) else list(c(-Inf, loc))
    names(bounds) <- stain
    rectangleGate(bounds, filterId=filterId)
}

## An old alias, this will go away soon
oneDGate <- function(...){
    .Deprecated(new="rangeGate")
    rangeGate(...)
}

## =================================================================================
## rangeFilter
## ---------------------------------------------------------------------------------
## This is basically an abstraction of the lymphGate function. It allows us
## to use it as a regular gate object.
## ---------------------------------------------------------------------------------
setClass("rangeFilter",
         representation=representation(alpha="character",
                                       sd="numeric",
                                       borderQuant="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultRangeFilter", 
         alpha="min",
         sd= 2,
         borderQuant=0.1))
         
## Constructor. We allow for the following inputs:
##  alpha is always a character and sd and borderQuant both are always
##     numerics of length 1.
##  stain is a list of characters and/or transformations, y is missing         
rangeFilter <- function(stain, alpha="min", sd=2, borderQuant=0.1, 
                       filterId="defaultRangeFilter")
{
    flowCore:::checkClass(filterId, "character", 1)
    flowCore:::checkClass(alpha, "character", 1)
    flowCore:::checkClass(sd, "numeric", 1)
    flowCore:::checkClass(borderQuant, "numeric", 1)        
    new("rangeFilter", parameters=stain, alpha=alpha,
        sd=sd, borderQuant=borderQuant, filterId=as.character(filterId))    
}

setMethod("%in%",
          signature=signature("flowFrame",
                              table="rangeFilter"),
          definition=function(x, table)
      {
          
          if(length(parameters(table)) != 1)
              stop("range filters require exactly one parameters.")
          tmp <- rangeGate(x, stain=parameters(table),
                           alpha=table@alpha,
                           sd=table@sd,
                           borderQuant=table@borderQuant,
                           filterId=table@filterId )
          filter(x, tmp)
        })


