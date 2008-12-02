## Do quad gating of a flowSet or flowFrame
quadrantGate <- function(x, stains, alpha = c("min", "min"), sd=c(2,2),
                         plot=FALSE, filterId="defaultQuadGate", ...)
{
    ## some type checking first
    flowCore:::checkClass(x, c("flowFrame", "flowSet"))
    flowCore:::checkClass(stains, "character", 2)
    mt <- stains %in% colnames(x)
    if(!all(mt))
        stop("Invalid stain(s) not matching the flowSet:\n    ",
             paste(stains[!mt], collapse=", "))
    alpha <- rep(alpha, length = 2)[1:2]
    flowCore:::checkClass(alpha, c("character", "numeric"), 2)
    sd <- rep(sd, length=2)[1:2]
    flowCore:::checkClass(sd, "numeric", 2)
    flowCore:::checkClass(plot, "logical", 1)
    flowCore:::checkClass(filterId, "character", 1)
    ## set up plotting device (if needed)
    if(plot){
        opar <- par(mfrow=c(2,2))
        on.exit(par(opar))
    }
    ## collapse to flowFrame and compute gates for each of the two channels
    if(is(x, "flowSet"))
        x <- as(x, "flowFrame")
    boundX <- density1d(x, stains[1], alpha = alpha[1], plot, sd=sd[1], ...)
    boundY <- density1d(x, stains[2], alpha = alpha[2], plot, sd=sd[2], ...)
    ## add final plot if needed
    if(plot){
        plot(x[,stains], main=paste("Quad-gate for parameters\n",
                         paste(stains, collapse=" and ")), cex.main=1)
        abline(h = boundY, v = boundX)
    }
    ## create the quadGate object
    gate <- c(boundX, boundY)
    names(gate) <- stains
    return(quadGate(.gate=gate, filterId=filterId))    
}

