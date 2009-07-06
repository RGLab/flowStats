## use backgeting to find the reference features and register the features
## used by function gpaSet

useBackGating <- function(bg, xy, plot=FALSE)
{
    ## 'bg' is the return value of function flowStats:::backGating

    regFeatures <- list()
    ch = split(bg, factor(bg$channel))
 
    for (i in 1:length(ch))
    { 
        perChannel <- ch[[i]]
       	centInd <- which(perChannel$type=="centroid")
    	center  <- perChannel[centInd, ]
        perPeak <- split(center, factor(center$peak)) 
        ## identify potential features from each peak and its population
	ref <- list()
        for (j in 1:length(perPeak)) 
	{
            if (length(unique(perPeak[[j]]$population))==1)
	        ref[[j]] <- data.frame(
                              matrix(apply(perPeak[[j]][, c("x", "y")], 2,
                                     median), nrow=1,
                                     dimnames=list(NULL, c("x", "y"))),
                              population=perPeak[[j]]$population[[1]],
 			      peak=perPeak[[j]]$peak[1],
                              channel=perPeak[[j]]$channel[1])
            else {
    	        eachPop <- split(perPeak[[j]], factor(perPeak[[j]]$population))
		## how many featurs in this peak population??
		## use cluster to find out: consider distance? and BIC?
                keepF <- sapply(eachPop, function(a, b) {
	                    (nrow(a) > 0.15 * b)}, 
                            nrow(perPeak[[j]]))
                ## ok for now, but ... 
                ## still have to do clustering to diffirenciate the clusters
                eachPop <- eachPop[keepF]

                ref[[j]] <- data.frame(t(sapply(eachPop, function(x) 
                             apply(x[, c("x", "y")], 2, median))),
                             population=sapply(eachPop,
                               function(x) x$population[1]),
                             peak=rep(perPeak[[j]]$peak[1], 
                                length(eachPop)),
                             channel=rep(perChannel$channel[1],
                                length(eachPop)))
	    }
        }
        
        refFeatures[[i]] <- do.call(rbind, ref)

    	#refFeatures <- sapply(eachPop, function(x) {
        #            km=kmeans(x[, c("x", "y")], 2)
	#	    f = km$center[which.max(km$size), ]})
        #refFeatures <- t(refFeatures)

        ## register features for each sample. Each sample's
    	## feature matrix must has the same dimension.
    	perSample = split(center, factor(center$sample))
    	nPeaks <- nrow(refFeatures[[i]])

    	## regFeatures <- lapply(perSample, registerFeatures, reference,
    	## by="population") 

    	regFeatures[[i]] <- 
             lapply(perSample, registerFeatures, refFeatures[[i]],
                    by="distance")

	     ## I want different color for each group
	if (plot) ## for debugging purposes
	{   
    	    channel = perChannel$channel[1]
    	    print(xyplot(x ~ y | factor(perChannel$sample), 
                data=perChannel, 
                main=paste("backgating on", channel)))
            
    	    print(xyplot(x~y|factor(center$population), 
                  data=center, 
                  main=paste("centriod per popupations", channel)))
            
    	    print(xyplot(x ~ y|which,
                do.call(make.groups, lapply(regFeatures[[i]], as.data.frame)),
                as.table=TRUE, 
                main=paste("registered features with backgating on", channel)))
            
            Features=regFeatures[[i]]
            Features$reference=refFeatures[[i]][, 1:2]
            print(xyplot(x ~ y,
                data=do.call(make.groups, lapply(Features, as.data.frame)),
                pch=0:(length(Features)-1), 
      	        group=which, col="black",
      	        key=list(space="right", title="samples",
                   text=list(names(Features)), 
                   points=list(pch=1:length(Features)-1 )),
      	        main=paste("registered features with backgating on", channel)))
         }

    }


    ## combine the members of regFeatures!
    combf <- list()
    for (m in names(regFeatures[[1]])) {
         combf[[m]] <- do.call(rbind, 
                         lapply(regFeatures, function(a, b) {a[[b]]}, m))
         rownames(combf[[m]]) <- as.character(1:nrow(combf[[m]]))
    }
    reference <- do.call(rbind, refFeatures)
    rownames(reference) <- as.character(1:nrow(reference))

    return(list(register=combf, reference=reference))
}
