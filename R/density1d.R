hubers1<-function (y, k = 1.5, mu, s, initmu = median(y), tol = 1e-06) 
{
	mmu <- missing(mu)
	ms <- missing(s)
	y <- y[!is.na(y)]
	n <- length(y)
	if (mmu) {
		mu0 <- initmu
		n1 <- n - 1
	}
	else {
		mu0 <- mu
		mu1 <- mu
		n1 <- n
	}
	if (ms) {
		s0 <- mad(y)
		if (s0 == 0) 
			return(list(mu = mu0, s = 0))
	}
	else {
		s0 <- s
		s1 <- s
	}
	th <- 2 * pnorm(k) - 1
	beta <- th + k^2 * (1 - th) - 2 * k * dnorm(k)
	for (i in 1:30) {
		yy <- pmin(pmax(mu0 - k * s0, y), mu0 + k * s0)
		if (mmu) 
			mu1 <- sum(yy)/n
		if (ms) {
			ss <- sum((yy - mu1)^2)/n1
			s1 <- sqrt(ss/beta)
		}
		if ((abs(mu0 - mu1) < tol * s0) && abs(s0 - s1) < tol * 
				s0) 
			break
		mu0 <- mu1
		s0 <- s1
	}
	list(mu = mu0, s = s0)
}
density1d_simple <- function(x, stain, alpha="min", sd=2, plot=FALSE, borderQuant=0.1,
		absolute=TRUE, inBetween=FALSE, refLine=NULL,rare=FALSE,bwFac = 1.2
		,sig=NULL
		,peakNr=NULL
		,cutoff=0.05
		, ...)
{
	
	## some type checking first
	flowCore:::checkClass(x, c("flowFrame", "flowSet"))
	flowCore:::checkClass(stain, "character", 1)
	if(!stain %in% colnames(x))
		stop("'", stain,"' is not a valid parameter in this flowFrame")
	
	flowCore:::checkClass(alpha, c("character", "numeric"), 1)
#	browser()
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
#   browser()
#   den<-density(exprs(x)[,stain])
	#get density curve peaks	
	fres<-filter(x,curv1Filter(stain,bwFac=bwFac))
	fres<-filterDetails(fres)
	fhat<-fres[[1]]$fsObj$fhat
	
#	bnds <- flowStats:::curvPeaks(fres, exprs(tmp)[, stain], borderQuant=borderQuant)
	#fitted curve points
	curve<-data.frame(x=fhat[[1]][[1]],y=fhat$est)
#	browser()
	curve_filtered<-subset(curve,y>=cutoff)
	#peak boudaries
	boundaries<-fres[[1]]$boundaries
	#filter peaks as well
	isNoise<-unlist(lapply(boundaries,function(curbn){
#				browser()
				subCurve<-subset(curve,x>=curbn[1]&x<=curbn[2])
				peak_mode<-subCurve[which.max(subCurve$y),]	
				peak_mode$y<cutoff
			}))
	boundaries<-boundaries[!isNoise]
	nPeaks<-length(boundaries)
#	plot(curve)
#	plot(curve_filtered)
#	browser()
	if(rare||nPeaks==1)
	{
		if(nPeaks==1)
		{
			peak_mode<-curve_filtered[which.max(curve_filtered$y),]
			#duplicate the lhs of the peak and try to get correct estimation of mean and sd
			#without the distortion of another peak
			if(!is.null(sig))
				if(sig=="L")
				{
					
					halfS<-curve[curve$x<=peak_mode[,"x"],]
					
				}else if(sig=="R")
				{
					halfS<-curve[curve$x>=peak_mode[,"x"],]
				}else
					stop("invalid value for 'sig' !")
#			plot(halfS)
		}else  
		{
			##pick LHS of far left peak or RHS of far right peak 
			#when there are more than 1 peak
			if(!is.null(sig))
				if(sig=="L")
				{
					lpeak_bnd<-boundaries[[1]]
					subCurve<-subset(curve,x>=lpeak_bnd[1]&x<=lpeak_bnd[2])
					peak_mode<-subCurve[which.max(subCurve$y),]	
					halfS<-curve[curve$x<=peak_mode$x,]
				}else if(sig=="R")
				{
					rpeak_bnd<-boundaries[[nPeaks]]
					subCurve<-subset(curve,x>=rpeak_bnd[1]&x<=rpeak_bnd[2])
					peak_mode<-subCurve[which.max(subCurve$y),]	
					halfS<-curve[curve$x>=peak_mode$x,]
				}else
					stop("invalid value for 'sig' !")
#			plot(halfS)
			
		}
		if(is.null(sig))
		{
			#use entire signal to estimate mode and sd
			est <- hubers(curve$x)
		}else
		{
#			browser()
			sd1<-sqrt(mean((halfS$x-peak_mode[,"x"])^2*halfS$y))*2
			
			est<-list(mu=peak_mode[,"x"],s=sd1)	
		}
		
		loc <- est$mu + sd * est$s
			
	}else
	{
		#pick lowest valley
		valleys<-do.call(rbind
				,lapply(1:(nPeaks-1),function(i){
							#		browser()
							#get valley boundaries
							valleyLbound<-max(boundaries[[i]])
							valleyRbound<-min(boundaries[[i+1]])
							
							#
							subCurve<-subset(curve,x>=valleyLbound&x<=valleyRbound)
							#		plot(subCurve)
							subCurve[which.min(subCurve$y),]
							
						})
		)   	
		
		loc<-valleys[which.min(valleys$y),"x"]
			
			
	}
	## create output if needed
#	browser()
	if(plot){
		plot(curve,type="l", main=paste("breakpoint for parameter", stain), cex.main=1, ...)
		abline(v = loc, col = 2, lwd = 2)
		
	}
#	browser()
	return(loc)
}

## Find most likely separator between peaks in 1D
density1d <- function(x, stain, alpha="min", sd=2, plot=FALSE, borderQuant=0.1,
                      absolute=TRUE, inBetween=FALSE, refLine=NULL,rare=FALSE,bwFac = 1.2
			  		,sig=NULL
					,peakNr=NULL
			  		, ...)
{
	
    ## some type checking first
    flowCore:::checkClass(x, c("flowFrame", "flowSet"))
    flowCore:::checkClass(stain, "character", 1)
    if(!stain %in% colnames(x))
        stop("'", stain,"' is not a valid parameter in this flowFrame")
	
    flowCore:::checkClass(alpha, c("character", "numeric"), 1)
#	browser()
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
#   browser()
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
						   
    fres <- filter(tmp <- Subset(x, exprf), curv1Filter(stain,bwFac=bwFac))
    bnds <- flowStats:::curvPeaks(fres, exprs(tmp)[, stain], borderQuant=borderQuant)
	##when peakNr is present,drop the less significant peaks by their heights
	
	if(!is.null(peakNr))
	{
		if(peakNr<nrow(bnds$peaks))
		{
			bigPks<-order(bnds$peaks[,"y"],decreasing=T)[1:peakNr]
			bnds$peaks<-bnds$peaks[bigPks,]
			bnds$regions<-bnds$regions[bigPks,]
			bnds$densFuns<-bnds$densFuns[bigPks]
			
		}
	}
	
    ## define 'positive' and 'negative' peaks by proximity to anchors
    anchors <- quantile(as.numeric(vrange), c(0.25, 0.5))
    anc <- abs(sapply(anchors, "-", bnds$peaks[, "x", drop=FALSE]))
    dens <- density(exprs(tmp)[, stain])
    ## only one peak: use robust estimation of mode and variance for boundary
	signals<-exprs(tmp[, stain])
		
#		browser()
    if(is.null(nrow(anc)))
    {
		if(!is.null(sig))
		{
			if(sig=="L")
			{
				halfS<-signals[signals<=bnds$peak[,"x"]]
				
			}
			if(sig=="R")
			{
				halfS<-signals[signals>=bnds$peak[,"x"]]
			}
			
			varS<-sqrt(mean((halfS-bnds$peak[,"x"])^2))
			
			est<-list(mu=bnds$peak[,"x"],s=varS)
		}else
		{
			est <- hubers(signals)
			
		}
		
        class <- which.min(abs(est$mu - anchors))
    }else 
	{
		est <- hubers(signals)
		if(nrow(anc)==2 && inBetween)
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
	}
	
	
		
	##if rare==true then assume there is only one major peak and the rare population is at the right tail
#    browser()
	if(rare||all(class==1))
	{
		bound <- est$mu + sd * est$s
		
	}
    else if(all(class==2))
        bound <- est$mu - sd * est$s
    else{ 
        if(is.character(alpha)){
            if(alpha=="min"){
                sel <- (dens$x > left & dens$x <= right)
                bound <- dens$x[sel][which.min(dens$y[sel])]
            }else
			{
				alpha<-as.numeric(alpha)
				if(is.na(alpha))
					stop("unknown value for alpha")
				else
					bound <- (1-alpha)*left + alpha*right
			}
            
        }else
        bound <- (1-alpha)*left + alpha*right    
    }
    ## create output if needed
#	browser()
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
                     positive=TRUE, refLine=NULL,simple=FALSE, ...)
{
	if(simple)
		loc <- density1d_simple(x=x, stain=stain, alpha=alpha, sd=sd, plot=plot,
				borderQuant=borderQuant, absolute=absolute, refLine=refLine, ...)
	else	
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


