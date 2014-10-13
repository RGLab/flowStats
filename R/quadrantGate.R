## Do quad gating of a flowSet or flowFrame
quadrantGate <- function(x, stains, alpha = c("min", "min"), sd=c(2,2),
                         plot=FALSE, filterId="defaultQuadGate",
                         refLine.1=NULL, refLine.2=NULL
				 		,rare=c(FALSE,FALSE)
				 		,sig=c(NULL,NULL)
						,...)
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
	rare <- rep(rare, length=2)[1:2]
	sig <- rep(sig, length=2)[1:2]
	
    flowCore:::checkClass(sd, "numeric", 2)
    
    if (!is.null(refLine.1))
        flowCore:::checkClass(refLine.1, "numeric",1)
    if (!is.null(refLine.2))
        flowCore:::checkClass(refLine.2, "numeric",1)
    
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
    boundX <- density1d(x, stains[1], alpha = alpha[1], plot, sd=sd[1],
                        refLine = refLine.1,rare=rare[1],sig=sig[1], ...)
    boundY <- density1d(x, stains[2], alpha = alpha[2], plot, sd=sd[2],
                        refLine = refLine.2,rare=rare[2], sig=sig[2],...)
				

    ## add final plot if needed
    if(plot){
        ## FIXME: Ugly hack foe name space issue
        plot <- selectMethod("plot", c("flowFrame", "character"))
        plot(x, stains, main=paste(basename(description(x)$FIL),"\n Quad-gate for parameters\n",
                                   paste(stains, collapse=" and ")), cex.main=1)
        abline(h = boundY, v = boundX)
    }
    ## create the quadGate object
    gate <- c(boundX, boundY)
    names(gate) <- stains
    return(quadGate(.gate=gate, filterId=filterId))    
}

