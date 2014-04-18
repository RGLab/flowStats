## normalizes data per-channel basis by locating high density regions
## on each channel and use these as landmarks to adjust the distance
## between populations.  Inputs: -flowset: the input flowset
## containing the flowcytomery expriments.  -channel.names: the list
## of channels to be normalized.  -max.lms: a vector in which each
## element is the maximum number of landmarks used in the
## normalization of the corresponding channel. The default values is 2
## for all channels.  -debug: if set TRUE draws before and after
## normalization plots for each sample. The plot of the i-th sample is
## stored in "paste(fname, i)" file.  -base.lms: list of base
## landmarks for each channel. i.e. base.lms[[p]] is a vector of base
## landmars for channel p.  If not specified the base.lms is computed.
## The peaks with density value less than peak.density.thr*maximum.peak.density are discarded.
## and of the peaks with distance less than peak.distance.thr*range.data only one is considered.

##Normalize a GatingSet, allowing for normalization on a subset of gates
gaussNormGS <- function(gatingset, channel.names, max.lms=2, base.lms=NULL,peak.density.thr=0.05,peak.distance.thr=0.05,debug=FALSE,fname='',gate=NULL){
	if(!inherits(gatingset,"GatingSet"))
		stop("gatingset must be of class GatingSet")
	if(!gatingset[[1]]@isNcdf)
		stop("gatingset not gated using netcdf. Just use the regular gaussNorm function, please");
	if(is.null(gate)){
		ncflowset<-graph::nodeData(gatingset[[1]]@tree,gatingset[[1]]@nodes[1],"data")[[1]][["data"]]$ncfs
	}else{
		ncflowset<-Subset(gatingset,gate)
	}
	#Update the column names on the ncflowset so that they are consistent with the column names on the GatingSet
	colnames(ncflowset)<-colnames(getData(gatingset[[1]]))
	normed<-gaussNormNC(ncflowset=ncflowset,channel.names=channel.names,max.lms=max.lms,base.lms=base.lms,peak.density.thr=peak.density.thr,peak.distance.thr=peak.distance.thr,debug=debug,fname=fname);
	return(normed);
}
## Output: -returns the normalized flowset.
gaussNormNC <- function(ncflowset, channel.names, max.lms=2, base.lms=NULL,  peak.density.thr=0.05, peak.distance.thr=0.05, debug=FALSE, fname=''){
  remb.flowset=remBoundaryNC(ncflowset, channel.names) #This method needs to be rewritten to deal with large data sets.
  expr.list=c(1:length(	ncflowset))  
  normalized<-ncdfFlow::clone.ncdfFlowSet(ncflowset)
  if(length(max.lms)==1){
    max.lms=rep(max.lms, times=length(channel.names))
  }
  if(length(max.lms)!=length(channel.names)){
    cat('Error: length of max.lms and channel.names doesn\'t match\n')
    return(NULL)
  }
  names(max.lms)=channel.names
  lms=extract.landmarksNC(ncflowset,remb.flowset, expr.list, channel.names, max.lms,  peak.density.thr, peak.distance.thr)
  if(is.null(base.lms))
    base.lms=extract.base.landmarks(lms$filter, channel.names, max.lms)
  
  ## finds the best matching between landmarks and base landmarks.
  ## the score of a matching is defined as a function of the distance between landmark and the base landmark and the score of the landmark itself.
  matched.lms=match.all.lms(lms, base.lms, channel.names, max.lms)
  confidence=compute.confidence(matched.lms, base.lms)
  cat('\nAdjusting the distance between landmarks\n')  
  newRange=matrix(ncol=length(channel.names), nrow=2)
  colnames(newRange)=channel.names
  newRange[,channel.names] <- c(Inf, -Inf)  
  for(i in expr.list){
    cat('.')
    if(fname != '')
      file.name=paste(fname, as.character(i), sep='.')
    else
      file.name=''
    ## normalizing sample i. 	
    tmp<-flowFrame(normalize.one.exprNC(ncflowset,remb.flowset$index[[i]],i, base.lms, lms$filter[,i], lms$original[,i], matched.lms[,i], channel.names, max.lms, file.name, debug))
	parameters(tmp)<-parameters(ncflowset[[i]])
	description(tmp)<-description(ncflowset[[i]])
	normalized[[sampleNames(ncflowset)[i]]] <- tmp
    for(p in channel.names){
      newRange[1,p] <- min(newRange[1,p], min(exprs(normalized[[i]])[,p], na.rm=TRUE))
      newRange[2,p] <- max(newRange[2,p], max(exprs(normalized[[i]])[,p], na.rm=TRUE))
    } 
  }
  cat('\n')
  restoreBoundaryNC(ncflowset,normalized, remb.flowset, channel.names)
  #updating the new ranges.
  for(i in expr.list){
    for(p in channel.names){
      ip <- match(p, pData(parameters(normalized[[i]]))$name)
      tmp <- parameters(normalized[[i]])
      oldRanges <- unlist(range(ncflowset[[i]],p))
      pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], newRange[1,p]),
                                                     max(oldRanges[2], newRange[2,p]))
      normalized[[i]]@parameters <- tmp
      
    }
  }
  list(flowset=normalized, confidence=confidence)
}

gaussNorm <- function(flowset, channel.names, max.lms=2, base.lms=NULL,  peak.density.thr=0.05, peak.distance.thr=0.05, debug=FALSE, fname=''){
  remb.flowset=remBoundary(flowset, channel.names)
  expr.list=c(1:length(flowset))  
  if(length(max.lms)==1){
    max.lms=rep(max.lms, times=length(channel.names))
  }
  if(length(max.lms)!=length(channel.names)){
    cat('Error: length of max.lms and channel.names doesn\'t match\n')
    return(NULL)
  }
  names(max.lms)=channel.names
  lms=extract.landmarks(remb.flowset$res, expr.list, channel.names, max.lms,  peak.density.thr, peak.distance.thr)
  if(is.null(base.lms))
    base.lms=extract.base.landmarks(lms$filter, channel.names, max.lms)
  
  ## finds the best matching between landmarks and base landmarks.
  ## the score of a matching is defined as a function of the distance between landmark and the base landmark and the score of the landmark itself.
  matched.lms=match.all.lms(lms, base.lms, channel.names, max.lms)
  confidence=compute.confidence(matched.lms, base.lms)
  cat('\nAdjusting the distance between landmarks\n')  
  newRange=matrix(ncol=length(channel.names), nrow=2)
  colnames(newRange)=channel.names
  newRange[,channel.names] <- c(Inf, -Inf)  
  for(i in expr.list){
    cat('.')
    if(fname != '')
      file.name=paste(fname, as.character(i), sep='.')
    else
      file.name=''
    ## normalizing sample i.  
    exprs(remb.flowset$res[[i]])=normalize.one.expr(exprs(remb.flowset$res[[i]]), base.lms, lms$filter[,i], lms$original[,i], matched.lms[,i], channel.names, max.lms, file.name, debug)
    for(p in channel.names){
      newRange[1,p] <- min(newRange[1,p], min(exprs(remb.flowset$res[[i]])[,p], na.rm=TRUE))
      newRange[2,p] <- max(newRange[2,p], max(exprs(remb.flowset$res[[i]])[,p], na.rm=TRUE))
    } 
  }
  cat('\n')
  restoreBoundary(flowset, remb.flowset, channel.names)
  #updating the new ranges.
  for(i in expr.list){
    for(p in channel.names){
      ip <- match(p, pData(parameters(remb.flowset$res[[i]]))$name)
      tmp <- parameters(remb.flowset$res[[i]])
      oldRanges <- unlist(range(flowset[[i]],p))
      pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], newRange[1,p]),
                                                     max(oldRanges[2], newRange[2,p]))
      remb.flowset$res[[i]]@parameters <- tmp
      
    }
  }
  list(flowset=remb.flowset$res, confidence=confidence)
}

## shifts the data in such a way that the peak at matched.lms[i] is moved to matched.lms[i+1] for each i.
## base.lms, lms.lis and lms.original are passed for the debug purposes only. 
##Works on ncdfFlowSet.
normalize.one.exprNC <- function(data,indices,i,base.lms,lms.list,lms.original,matched.lms,channel.names,max.lms,fname='',debug=FALSE){
	norm.data=exprs(data[[i]]);
	norm.data[indices,]<-NA;
	if(debug==TRUE){
		if(fname=='')
		x11()
		else
		png(filename=fname, bg="white",width=1800,height=1000)      
		par(mfrow=c(1, length(which(max.lms!=0))))
	}
	cn=1
	i=0
	## normalizing each channel.
	for (p in channel.names){
		i=i+1
		if(max.lms[i]==0){
			norm.data[,p]=data[[i]][,p]
		}
		else{
			## registering the data of channel p given the matching landmarks.
			norm.data[,p]=register.channel(norm.data[,p], matched.lms[[p]])
			## drawing the pre- and post- normalization density curves.
			if(debug==TRUE){
				if(length(lms.list[[p]])==0){
					A1=density(na.omit(data[,p]))
					xlim=c(min(A1$x, na.rm=T), max(A1$x, na.rm=T))
					ylim=c(min(A1$y, na.rm=T), max(A1$y, na.rm=T))
					par(mfg=c(1, cn))
					plot(A1, type="l", col= "black", xlim=xlim, ylim=ylim, xlab=p, main='No lms is found')
				}
				else{
					A1=density(na.omit(data[,p]))
					A2=density(na.omit(norm.data[,p]))      
					xlim=c(min(min(A1$x, na.rm=T), min(A2$x, na.rm=T)), max(max(A1$x, na.rm=T), max(A2$x, na.rm=T)))
					ylim=c(min(min(A1$y, na.rm=T), min(A2$y, na.rm=T)), max(max(A1$y, na.rm=T), max(A2$y, na.rm=T)))
					par(mfg=c(1, cn))
					plot(A1, type="l", col= "blue", xlim=xlim, ylim=ylim, xlab=p, main='')
					Lms=matched.lms[[p]][seq(1, length(matched.lms[[p]]), by=2)]
					M.Lms=matched.lms[[p]][seq(2, length(matched.lms[[p]]), by=2)]
					M.Lms.index=1
					for(k in 1:(length(M.Lms)))
					M.Lms.index[k]=which(base.lms[[p]]==M.Lms[k])
					points(Lms, returny(A1, Lms), pch=19, col="red")          
					text(Lms, returny(A1, Lms), as.character(M.Lms.index), pos=3, col='red')          
					points(lms.original[[p]], returny(A1, lms.original[[p]]), pch=21, col="black") ##all the peaks
					par(mfg=c(1, cn))
					plot(A2, type="l", col="red", xlim=xlim, ylim=ylim, xlab=p)   
					points(M.Lms, returny(A2, M.Lms), pch=15, col="blue")          
					text(M.Lms, returny(A2, M.Lms), as.character(M.Lms.index), pos=3, col='blue')          
				}
				cn=cn+1        
			}      
		}
	}
	if(debug){
		if(fname=='')
		readline()
		dev.off()
	}
	return (norm.data)
}


## shifts the data in such a way that the peak at matched.lms[i] is moved to matched.lms[i+1] for each i.
## base.lms, lms.lis and lms.original are passed for the debug purposes only. 
normalize.one.expr <- function(data, base.lms, lms.list, lms.original, matched.lms, channel.names, max.lms, fname='', debug=FALSE){
  norm.data=data
  if(debug==TRUE){
    if(fname=='')
      x11()
    else
      png(filename=fname, bg="white",width=1800,height=1000)      
    par(mfrow=c(1, length(which(max.lms!=0))))
  }
  cn=1
  i=0
  ## normalizing each channel.
  for (p in channel.names){
    i=i+1
    if(max.lms[i]==0){
      norm.data[,p]=data[,p]
    }
    else{
      ## registering the data of channel p given the matching landmarks.
      norm.data[,p]=register.channel(data[,p], matched.lms[[p]])
      ## drawing the pre- and post- normalization density curves.
      if(debug==TRUE){
        if(length(lms.list[[p]])==0){
          A1=density(na.omit(data[,p]))
          xlim=c(min(A1$x, na.rm=T), max(A1$x, na.rm=T))
          ylim=c(min(A1$y, na.rm=T), max(A1$y, na.rm=T))
          par(mfg=c(1, cn))
          plot(A1, type="l", col= "black", xlim=xlim, ylim=ylim, xlab=p, main='No lms is found')
        }
        else{
          A1=density(na.omit(data[,p]))
          A2=density(na.omit(norm.data[,p]))      
          xlim=c(min(min(A1$x, na.rm=T), min(A2$x, na.rm=T)), max(max(A1$x, na.rm=T), max(A2$x, na.rm=T)))
          ylim=c(min(min(A1$y, na.rm=T), min(A2$y, na.rm=T)), max(max(A1$y, na.rm=T), max(A2$y, na.rm=T)))
          par(mfg=c(1, cn))
          plot(A1, type="l", col= "blue", xlim=xlim, ylim=ylim, xlab=p, main='')
          Lms=matched.lms[[p]][seq(1, length(matched.lms[[p]]), by=2)]
          M.Lms=matched.lms[[p]][seq(2, length(matched.lms[[p]]), by=2)]
          M.Lms.index=1
          for(k in 1:(length(M.Lms)))
            M.Lms.index[k]=which(base.lms[[p]]==M.Lms[k])
          points(Lms, returny(A1, Lms), pch=19, col="red")          
          text(Lms, returny(A1, Lms), as.character(M.Lms.index), pos=3, col='red')          
          points(lms.original[[p]], returny(A1, lms.original[[p]]), pch=21, col="black") ##all the peaks
          par(mfg=c(1, cn))
          plot(A2, type="l", col="red", xlim=xlim, ylim=ylim, xlab=p)   
          points(M.Lms, returny(A2, M.Lms), pch=15, col="blue")          
          text(M.Lms, returny(A2, M.Lms), as.character(M.Lms.index), pos=3, col='blue')          
        }
        cn=cn+1        
      }      
    }
  }
  if(debug){
    if(fname=='')
      readline()
    dev.off()
  }
  return (norm.data)
  
}

#Works on ncdfFlowSet rather than flowSet
## computes the channel landmarks for each flowset experiment in expr.list.
## max.lms is the maximum number of landmarks we want to shift when normalizing the data.
## for each landmark a score is assigned which is a function of its sharpness and its density value.
## output:
##       lms.list$original: list of all landmarks
##       lms.list$score:    the score of lms.list$original landmarks
##       lms.list$filtered: max.lms top score landmarks.
extract.landmarksNC<-function(ncflowset,indices,expr.list, channel.names, max.lms,  peak.density.thr, peak.distance.thr){
	## defining output variables.
	lms.list=list()    
	lms.list$original=matrix(vector("list"), length(channel.names), length(expr.list))
	lms.list$score=matrix(vector("list"), length(channel.names), length(expr.list))  
	lms.list$filter=matrix(vector("list"), length(channel.names), length(expr.list))
	row.names(lms.list$original)=channel.names
	row.names(lms.list$score)=channel.names  
	row.names(lms.list$filter)=channel.names  
	max.lms.num=0
	c=1
	## iterating over channels.
	for( p in channel.names){
		if(max.lms[c]!=0){
			j=1
			for(i in expr.list){        
				data=exprs(ncflowset[[i]])
				## finding the landmarks of channel p of sample i.
				data[indices$index[[i]],]<-NA
				lms=landmarker(data, p, max.lms[c])
				lms.list$original[p, j]=list(lms)
				## returns the max.lms[c] top score landmarks.
				filtered=filter.lms(lms, data[,p], max.lms[c],  peak.density.thr, peak.distance.thr)
				lms.list$filter[p, j]=list(filtered$lms)
				lms.list$score[p, j]=list(filtered$score)
				j=j+1
			}
		}
		c=c+1
	}
	return(lms.list)
}

## computes the channel landmarks for each flowset experiment in expr.list.
## max.lms is the maximum number of landmarks we want to shift when normalizing the data.
## for each landmark a score is assigned which is a function of its sharpness and its density value.
## output:
##       lms.list$original: list of all landmarks
##       lms.list$score:    the score of lms.list$original landmarks
##       lms.list$filtered: max.lms top score landmarks.

extract.landmarks <- function(flowset, expr.list, channel.names, max.lms,  peak.density.thr, peak.distance.thr){
  ## defining output variables.
  lms.list=list()    
  lms.list$original=matrix(vector("list"), length(channel.names), length(expr.list))
  lms.list$score=matrix(vector("list"), length(channel.names), length(expr.list))  
  lms.list$filter=matrix(vector("list"), length(channel.names), length(expr.list))
  row.names(lms.list$original)=channel.names
  row.names(lms.list$score)=channel.names  
  row.names(lms.list$filter)=channel.names  
  max.lms.num=0
  c=1
  ## iterating over channels.
  for( p in channel.names){
    if(max.lms[c]!=0){
      j=1
      for(i in expr.list){        
        data=exprs(flowset[[i]])
        ## finding the landmarks of channel p of sample i.
        lms=landmarker(data, p, max.lms[c])
        lms.list$original[p, j]=list(lms)
        ## returns the max.lms[c] top score landmarks.
        filtered=filter.lms(lms, data[,p], max.lms[c],  peak.density.thr, peak.distance.thr)
        lms.list$filter[p, j]=list(filtered$lms)
        lms.list$score[p, j]=list(filtered$score)
        j=j+1
      }
    }
    c=c+1
  }
  return(lms.list)
}

## computes the base landmarks.
## output: a vector of base landmarks for each channel p.
extract.base.landmarks <- function(lms, channel.names, max.lms){
  lms.list=list()
  max.lms.num=0
  for( p in channel.names){
    if(max.lms[p]!=0){
      ## if the total number of landmarks for channel p is less than max.lms[p] return them as base landmarks.
      if(length(unlist(lms[p,])) <= max.lms[p]){
        lms.list[[p]]=sort(unlist(lms[p,]))
      }
      else{
        
        if (max.lms[p]==1){
          lms.list[[p]]=median(unlist(lms[p,]))
        }
        else {
          ## first identify samples that have exactly max.lms[p] landmarks on their channel p. These landmarks are labeled from 1 to max.lms
          ## the landmarks samples with less than max.lms[p] landmarks are stored in short.lms vector.
          short.lms=vector()
          lms.class=list()
          for(jj in 1:max.lms[p])
            lms.class[[jj]]=vector()
          for (ii in 1:length(lms[p,])){
            lms.p.ii=lms[p,ii][[1]]
            if(length(lms.p.ii)==max.lms[p]){
              for (jj in 1:max.lms[p])
                lms.class[[jj]]=append(lms.class[[jj]], lms.p.ii[jj])
            }
            else{
              short.lms=append(short.lms, lms.p.ii)
            }
          }

          if(length(lms.class[[1]])==0){
            cat('No curve with ', max.lms[p], ' landmarks was found for channel ', p, '. Decrease max.lms[',p,'] and try again.\n')
            stop()
          }
          mean.lms.class=0
          for(jj in 1:max.lms[p]){
            mean.lms.class[jj]=mean(lms.class[[jj]])
          }

          ## assigning short.lms landmarks to the class of labeled landmarks wich are closer to. 
          if(length(short.lms)!=0){
            for(jj in 1:length(short.lms)){
              kk=which(abs(mean.lms.class-short.lms[jj])==min(abs(mean.lms.class-short.lms[jj])))
              kk=kk[1]
              lms.class[[kk]]=append(lms.class[[kk]], short.lms[jj])
            }
          }
          lms.class.len=0
          lms.class.med=0
          for(jj in 1:max.lms[p]){
            lms.class.len[jj]=length(lms.class[[jj]])
            lms.class.med[jj]=median(lms.class[[jj]], na.rm=T)            
          }
          s=sd(lms.class.len, na.rm=T)
          m=mean(lms.class.len, na.rm=T)
          ll=which(m-lms.class.len > 2*s)
          ## if a class of landmarks has too few landmarks in it just ignore this class.
          if(length(ll)!=0){
            cat('warning: fewer landmark classes found in channel ', p, '\n')
            tmp.lms=list()
            kk=1
            for(jj in (1:max.lms[p])[-ll]){
              tmp.lms[[kk]]=lms.class[[jj]]
              kk=kk+1
            }
            lms.class=tmp.lms
            
          }

          ## returning the median of each class as the base landmark.
          for(jj in 1:length(lms.class)){
            lms.class.med[jj]=median(lms.class[[jj]], na.rm=T)            
          }

          lms.list[[p]]=lms.class.med
          
        }
      }
    }
  }
  return(lms.list)
}

## finds the best matching between landmarks (lms) and the base landmarks (base.lms).
## the score of a matching is defined as a function of the distance between the base landmark and the landmark and the score of the landmark itself.
## returns: matched.lms-for each channel p of a sample i matched.lms[p, i] is a list of pairs of the form (lms, base.lms)
match.all.lms <- function(lms, base.lms, channel.names, max.lms){  
  n=length(lms$filter[channel.names[1],]) ##number of samples.
  matched.lms=matrix(vector("list"), length(channel.names), n)
  row.names(matched.lms)=channel.names
  lms.class=list()
  lms.class.median=list()
  for(p in channel.names){
    lms.class[[p]]=list()
    lms.class.median[[p]]=list()
  }
  for (p in channel.names){
    if(max.lms[[p]]==0)
      next
    for(i in 1:n){
      matched.lms[p, i][[1]]=best.match(lms$original[p, i][[1]], lms$score[p, i][[1]], base.lms[[p]], max.lms[[p]])
    }
  }
  return(matched.lms)  
}

best.match <- function(lms, lms.score, base.lms, max.lms){
  comb=choose.nk(max(length(lms), length(base.lms)), min(length(lms), length(base.lms)))
  sc=list()
  k=1
  max.s=-1
  if(length(lms)<length(base.lms)){
    for(i in 1:(dim(comb)[1])){
      if(length(which(comb[i,]==0))!=0)
        c=comb[i,][-which(comb[i,]==0)]
      else
        c=comb[i,]      
      d=combinations.itr(length(comb[i,]), length(c))
      for(j in 1:(dim(d)[1])){
        s=match.score(lms[d[j,]], base.lms[c], lms.score[d[j,]])
        if(max.s<s){
          lms.index=d[j,]
          base.lms.index=c
          max.s=s
        }
        sc[[k]]=list(score=s, lms.index=d[j,], base.lms.index=c)
        k=k+1
      }
    }
    
  }
  else{
    for(i in 1:(dim(comb)[1])){
      if(length(which(comb[i,]==0))!=0)      
        c=comb[i,][-which(comb[i,]==0)]
      else
        c=comb[i,]
      d=combinations.itr(length(comb[i,]), length(c))
      for(j in 1:(dim(d)[1])){
        s=match.score(lms[c], base.lms[d[j,]], lms.score[c])
        if(max.s<s){
          lms.index=c
          base.lms.index=d[j,]
          max.s=s
        }
        
        k=k+1        
      }      
    }    
  }
  res=1:(2*length(lms.index))
  res[seq(1,2*length(lms.index), by=2)]=lms[lms.index]
  res[seq(2,2*length(lms.index), by=2)]=base.lms[base.lms.index]
  return(res)
}

## defines the score of a matching between landmarks (lms) and base landmarks (base.lms)
match.score <- function(lms, base.lms, lms.score){
  d=abs(lms-base.lms)
  if(length(which(d==0))!=0)
    d[which(d==0)]=0.01
  return(sum(1/d*lms.score*lms.score))
}

choose.nk <- function(n, k){
  v=matrix(0, ncol=k, nrow=0)
  for(j in 1:k){
    c=combinations.itr(n, j)
    if(j<k)
      c=cbind(c, matrix(0, nrow=dim(c)[1], ncol=k-j))
    v=rbind(v, c)
  }
  return(v)    
}

combinations.itr <- function(n, k){
  v=1
  res=NULL
  while(length(v)!=0){
    if(v[length(v)]>n){
      v=v[-length(v)]
      v[length(v)]=v[length(v)]+1
    }
    else if(length(v)==k){

      res=rbind(res, v)
      v[length(v)]=v[length(v)]+1
    }
    else{
      v=append(v, v[length(v)]+1)
    }
  }
  return(res)
}

## manipulates the data in such a way that the landmark at matched.lms[i] is moved to matched.lms[i+1] for each i.
register.channel <- function(data, matched.lms){
  if(length(matched.lms)==0){
    cat('*')
    return (data)
  }
  s=m=shift=vector()
  lms=vector()
  for(i in seq(1,length(matched.lms), by=2)){
    shift=append(shift, matched.lms[i+1]-matched.lms[i])
    lms=append(lms, matched.lms[i])
    s=append(s, sd(na.omit(data)))
  }
  r.data=register.function(data, s, lms, shift)
  return(r.data)
}

returny <- function(A, X){
  y=vector()
  i=1;
  for( x in X){
    y[i]=A$y[which(abs(A$x-x)==min(abs(A$x-x)))[1]]
    i=i+1
  }
  return (y)
}

register.function <- function(data, s, m, shift){    
  sum=0
  if(length(m)==1){
    return(data+shift)
  }
  if(length(m)==2){
    sh=(shift[1]-shift[2])
    data=data+gau(data, s[1], m[1])*(sh/2)
    data=data-gau(data, s[2], m[2])*(sh/2)
    return(data+shift[1]-sh/2)
  }
  max.shift=which(abs(shift)==max(abs(shift)))[1]
  if(shift[max.shift]>0){
    sh=(shift[max.shift]-(shift[max.shift]-min(shift[-max.shift]))/2)
  }
  else{
    sh=(shift[max.shift]-(shift[max.shift]-max(shift[-max.shift]))/2)
  }
  data=data+sh
  shift=shift-sh

  for(i in 1:length(m))
    data=data+gau(data, s[i], m[i])*shift[i]    
  return (data)
}


## gaussian function used in shifting the data.
gau <- function(d, s, m){
  return(2.7182^(-((d-m)^2)/(2*s*s)))
}

## computes the confidence value based on the matched landmarks.
compute.confidence <- function(matched.lms, base.lms){
  confidence=rep(1, times=length(matched.lms[1,]))
  # for each channel c.
  for(c in 1:length(base.lms)){
    for(l in 1:length(base.lms[[c]])){
      lms.list=0
      for ( i in 1:length(matched.lms[1,])){
        ind=which(matched.lms[c,i][[1]]==base.lms[[c]][l])
        if(length(ind)==0){
          lms.list[i]=NA
        }
        else{
          if(ind[1] %% 2 != 0)
            ind[1]=ind[1]+1
          lms.list[i]=matched.lms[c,i][[1]][ind[1]-1]
        }
      }
      d=density(na.omit(lms.list))
      den.lms=0
      for(i in 1:length(lms.list)){
        if(is.na(lms.list[i]))
          den.lms[i]=NA
        else
          den.lms[i]=d$y[which(abs(d$x-lms.list[i])==min(abs(d$x-lms.list[i])))]
      }
      m.den.lms=mean(den.lms, na.rm=T)
      ind1=which(den.lms>m.den.lms/4 & den.lms<m.den.lms/3)
      ind2=which(den.lms<m.den.lms/4)
      if(length(ind1)!=0)
        confidence[ind1]=confidence[ind1]*0.8
      if(length(ind2)!=0)
        confidence[ind2]=confidence[ind2]*0.7
      if(length(which(is.na(lms.list)))!=0)
        confidence[which(is.na(lms.list))]=confidence[which(is.na(lms.list))]*0.6
      
      #lms.frame=data.frame(lms=lms.list, ind=c(1:length(lms.list)))
      #lms.frame=lms.frame[do.call(order, c(lms.frame["lms"], decreasing=F)), ]
      #diff=unlist(lms.frame['lms'])[2:length(lms.list)]-unlist(lms.frame['lms'])[1:(length(lms.list)-1)]
      #bw=2*mean(na.omit(diff))
      #den=0
      #for(i in 1:length(lms.list)){
      #  if(is.na(lms.list[i]))
      #     den[i]=NA
      #  else
      #    den[i]=length(which(lms.list>lms.list[i]-bw & lms.list<lms.list[i]+bw))/length(lms.list)
      #}
    }
  }
  confidence
}

############## Landmark finding functions############

## returns the peaks (local maxima's) in the kernel density estimate of data.
landmarker <- function(data, channel.name, max.lms, span=3){
  d=data[,channel.name]
  A=density(na.omit(d))
  d=A$y
  pks=c()
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  for( i in 1:(length(d)-span)){
    if(!is.na(d[i+span%/%2]) & (d[i+span%/%2]==max(d[c(i:(i+span))], na.rm=T)))
      if(!is.na(d[i]) & !is.na(d[i+span]) & d[i+span%/%2]!=d[i] & d[i+span%/%2]!=d[i+span])
        pks=append(pks, i+span%/%2)
    
  }
  return(A$x[pks])
}


## returns the max.lms top score landmarks.
## the score of a landmarks is a funciton of its sharpness and density value
filter.lms <- function(lms, data, max.lms, peak.density.thr, peak.distance.thr){
  filtered=list()
  if(length(lms) == 0){
    filtered$lms=vector()
    filtered$score=vector()    
    return(filtered)
  }
  filtered$score=score.lms(lms, data, max.lms, peak.density.thr, peak.distance.thr)
  lms.score=data.frame(score=filtered$score, ind=c(1:length(lms)))
  lms.score=lms.score[do.call(order, c(lms.score["score"], decreasing=T)), ]
  ind=which(lms.score$score>0)
  if(length(ind)==0){
    filtered$lms=vector()
    filtered$score=vector()        
    return(filtered)
  }
  lms.score.ind=lms.score$ind[ind]
  if(length(lms.score.ind)<max.lms)
    filtered$lms=sort(lms[lms.score.ind], decreasing=F)
  else
    filtered$lms=sort(lms[lms.score.ind[c(1:max.lms)]], decreasing=F)
  return(filtered) 
}


## assigns a score to each landmark. the score of a landmarks is a funciton of its sharpness and density value.
## the peaks with density value less than peak.density.thr*maximum.peak.density are discarded.
## of the peaks with distance less than peak.distance.thr*range.data only one is considered.

score.lms <- function(lms, data, max.lms, peak.density.thr, peak.distance.thr){    
  bw=64
  score=vector()
  height.cutoff=peak.density.thr
  if(length(lms) == 0)
    return(score)
  A=density(na.omit(data))
  bw=min(64, length(A$x)/10)
  lms.max.height=max(returny(A, lms), na.rm=T)
  MIN.LMS.DIST=(max(A$x, na.rm=T)-min(A$x, na.rm=T))*peak.distance.thr
  last.lms=-1
  last.lms.i=-1
  last.lms.score=0
  for(i in 1:length(lms)){
    lms.ind=which(na.omit(abs(A$x-lms[i]))==min(na.omit(abs(A$x-lms[i])) ) )
    ind=(max(lms.ind-bw%/%2, 1)):(min(lms.ind+bw%/%2, length(A$x)))
    if(length(ind)==0)
      ind=1

    if(A$y[lms.ind] <   height.cutoff*lms.max.height){    
      w=0
    }
    else{
      ## computing the sharpness
      w=A$y[lms.ind]-A$y[ind]
      w[which(w<0)]=3*w[which(w<0)]
    }
    ## computing final score
    score[i]=sum(w, na.rm=T)*A$y[lms.ind]
    if(score[i]<0)
      score[i]=0
    if(last.lms<0){
      last.lms=lms[i]
      last.lms.i=i

    }
    else{
      ##If two lms's are very close only choose one with the better score.
      if(lms[i]-last.lms < MIN.LMS.DIST){
        if(score[i]>score[last.lms.i]){
          last.lms=lms[i]
          score[last.lms.i]=0
          last.lms.i=i          
        }
        else{
          score[i]=0
        }
      }
      else{
        last.lms=lms[i]
        last.lms.i=i
      }
    }
  }
  return(score)  
}

#Returns the indices of each flowFrame with boundary values that should be NA. In contrast to remBoundary, doesn't copy the data and doesn't return the NA'd flowSet
remBoundaryNC<-function(ncflowset,channel.names){
	ranges<-fsApply(ncflowset,range);
	index=list()
	res=list()
	for(i in 1:length(ncflowset)){
		d=exprs(ncflowset[[i]])
		index[[i]]=vector()
		for(p in channel.names){
			ind=which(d[,p]==ranges[[i]][1,p])
			index[[i]]=append(index[[i]],ind)
			#d[ind,p]=NA
			ind=which(d[,p]==ranges[[i]][2,p])
			index[[i]]=append(index[[i]],ind)
			#d[ind,p]=NA
		}
		#res[[i]]=new('flowFrame', exprs=d)
		#parameters(res[[i]])=parameters(flowset[[i]])
	}
	#names(res)=sampleNames(flowset)
	#res=as(res, 'flowSet')
	#phenoData(res)=phenoData(flowset)
	return(list(index=index))   
}

#setting the boundary values to NA. Also returns the indices of the NA values for recovery.
remBoundary <- function(flowset, channel.names){
  ranges <- fsApply(flowset, range)
  index=list()
  res=list()
  for(i in 1:length(flowset)){
    d=exprs(flowset[[i]])
    index[[i]]=vector()
    for(p in channel.names){
      ind=which(d[,p]==ranges[[i]][1,p])
      index[[i]]=append(index[[i]],ind)
      d[ind,p]=NA
      ind=which(d[,p]==ranges[[i]][2,p])
      index[[i]]=append(index[[i]],ind)
      d[ind,p]=NA
    }
    res[[i]]=new('flowFrame', exprs=d)
    parameters(res[[i]])=parameters(flowset[[i]])
  }
  names(res)=sampleNames(flowset)
  res=as(res, 'flowSet')
  phenoData(res)=phenoData(flowset)
  return(list(index=index, res=res))   
}

#Restore the boundary values that were NA in the ncdfFlowSet
restoreBoundaryNC <- function(org.flowset, normalized.flowset,remb.flowset,channel.names){
	for(i in 1:length(org.flowset)){
		for(p in channel.names){
			exprs(normalized.flowset[[i]])[remb.flowset$index[[i]], p]=exprs(org.flowset[[i]])[remb.flowset$index[[i]], p]
		}
	}
}

#Restore the boundary values that were set to NA.
restoreBoundary <- function(org.flowset, remb.flowset, channel.names){
  for(i in 1:length(org.flowset)){
    for(p in channel.names){
      exprs(remb.flowset$res[[i]])[remb.flowset$index[[i]], p]=exprs(org.flowset[[i]])[remb.flowset$index[[i]], p]
    }
  }
}

