## use backgeting to register and label features

useBackGating <- function(bg, xy, fig=FALSE)
{

    regFeatures <- list()
    ch = split(bg, factor(bg$channel))

    for (k in 1:length(ch))
    { 
        perChannel = ch[[k]]
    	head(perChannel)
    	centroid = which(perChannel$type=="centroid")
    	center=perChannel[centroid, ]

	## keep the population that has more than 15% of the centroid point
    	eachPop <- split(center, factor(center$population))
        keepF <- sapply(eachPop, function(a, b) {
	                   (nrow(a) > 0.10 * b) & (nrow(a) >1) }, 
                        nrow(center))
        eachPop <- eachPop[keepF]

    	refFeatures <- sapply(eachPop, function(x) {
                    k=kmeans(x[, c("x", "y")], 2)
		    f = k$center[which.max(k$size), ]})
    		    refFeatures <- t(refFeatures)

        ## register (labelling) features for each samples. each sample's
    	## feature matrix must has the same dimension.
    	perSample = split(center, factor(center$sample))
    	nPeaks <- nrow(refFeatures)

    	## regFeatures <- lapply(perSample, registerFeatures, refFeatures,
    	## by="population") 

    	regFeatures[[k]] <- 
             lapply(perSample, registerFeatures, refFeatures, by="distance")

	     ## I want different color for each group
	if (fig) ## for debugging only
	{   
    	    channel = perChannel$channel[1]
    	    print(xyplot(x ~ y | factor(perChannel$sample), 
                data=perChannel, 
                main=paste("backgating", channel)))
    	    print(xyplot(x~y|factor(center$population), 
                  data=center, 
                  main=paste("centriod per popupations", channel)))
    	    print(xyplot(x ~ y|which,
                  do.call(make.groups, lapply(regFeatures[[k]], as.data.frame)),
                  as.table=TRUE, 
                  main=paste("registered features", channel)))

            print(xyplot(x ~ y,
                do.call(make.groups, lapply(regFeatures[[k]], as.data.frame)),
      	        group=which, 
      	        auto.key=list(space="right"), 
      	        main=paste("registered features", channel)))
         }

    }
    ## combine the members of regFeatures!
    combf <- list()
    for (m in names(regFeatures[[1]])) {
         combf[[m]] <- do.call(rbind, 
                         lapply(regFeatures, function(a, b) {a[[b]]}, m))
         rownames(combf[[m]]) <- as.character(1:nrow(combf[[m]]))
    }
    
    return(combf)
}
