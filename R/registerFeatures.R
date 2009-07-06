registerFeatures <- function(eachSample, refFeatures, by="distance")
{  ## eachSample is output of backgating, see useBackGating.R for details

   ## debugging, need to change the size of reg
   reg <- matrix(NA, ncol=2, nrow=nrow(refFeatures), 
                 dimnames=list(rownames(refFeatures), c("x", "y")))
   nf = nrow(refFeatures)
   if (by == "population") 
   {
       npop <- split(eachSample, eachSample$population)
       ## 1. deal with single or multiple centroids for each popoulation 
       for (j in 1:length(npop))
       {  
          pop <- as.numeric(npop[[j]]$population)[1]
          if (nrow(npop[[j]]) == 1){
              reg[pop, ] = as.matrix(npop[[j]][1, c("x", "y")])
          }		   
          else {
             tmp = rbind(refFeatures[pop, ], npop[[j]][, c("x", "y")])
             ## get the one with smallest distance
 	     reg[pop, ] = 
               as.matrix(npop[[j]][which.min(dist(tmp)[1:nf]), c("x", "y")])
          }   
       }
   }
   else 
   {   ## by=="distance"
       ## between the centroid of each population and refFeatures. 
       eachSample$label <- apply(eachSample[, c("x", "y")], 1, 
           function(m, refF) {
                              tmp = rbind(m, refF[, c("x", "y")])
                              return(which.min(dist(tmp)[1:nf]))}, 
           refFeatures)
       nlab <- split(eachSample, factor(eachSample$label))
       for (j in 1:length(nlab))
       {
           lab <- as.numeric(nlab[[j]]$label)[1]
	   if (nrow(nlab[[j]]) == 1) {
       	       reg[lab, c("x", "y")] <- as.matrix(nlab[[j]][, c("x", "y") ])
	   }
	   else {
	       tmp <- rbind(refFeatures[lab, c("x", "y")],
                     nlab[[j]][, c("x", "y")])
	       reg[lab, c("x","y")] <- as.matrix(
	          nlab[[j]][which.min(dist(tmp)), c("x", "y") ])
           }
       }
       
    }
    ## 2. deal with NA peaks, give the refFeatures
    naIdx <- which(is.na(reg[, 1]))
    if (!identical(naIdx, integer(0)))
           reg[naIdx, ] <- as.matrix(refFeatures[naIdx, c("x", "y")])

    return(reg)
}
