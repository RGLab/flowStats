#Sequentially normalize a gating set by normalizing channels at each step of gating
#return a normalized gatingset
#target is the target sample to normalize against
#skip is a list of gate indices for which you want to skip normalization (i.e. dump gates that may be difficult to normalize, or FSC/SSC)
getNormGateList<-function(x){
	flowCore:::checkClass(x,"GatingSet")
	bfsgates<-lapply(x,function(y)which(sapply(RBGL::bfs(y@tree),function(x)!flowWorkspace:::.isBooleanGate.graphNEL(y,x))))
	bfsgates<-unique(do.call(rbind,bfsgates))
	return(data.frame(gate=flowWorkspace:::getNodes(x[[1]])[bfsgates],index=t(bfsgates)))
}
plotAllGates<-function(gh,cex=2,gsubset=NULL){
  nodes<-RBGL::bfs(gh@tree)
  nodes<-nodes[!sapply(nodes,function(x)flowWorkspace:::.isBooleanGate.graphNEL(gh,x))][-1L]
  if(is.null(gsubset))
    l<-length(nodes)
  else{
    nodes<-nodes[gsubset]
    l<-length(nodes)
  }
  rs<-floor(sqrt(l))
  cs<-l%/%rs
  while(rs*cs<l){
    if((rs+1)*cs<rs*(cs+1))
      rs<-rs+1
    else
      cs<-cs+1
  }
  grid<-expand.grid(1:rs,1:cs)
  sapply(1:(length(nodes)-1),function(g){
    print(flowWorkspace:::plotGate(gh,nodes[g],lwd=2,pch='.',cex=cex,main=strsplit(nodes[g],"\\.")[[1]][2]),split=c(grid[g,1],grid[g,2],rs,cs),more=T)
  })
  print(flowWorkspace:::plotGate(gh,nodes[length(nodes)],lwd=2,pch='.',cex=cex,main=strsplit(nodes[length(nodes)],"\\.")[[1]][2]),split=c(grid[length(nodes),1],grid[length(nodes),2],rs,cs),more=F)
}
plotSameGate<-function(gs,cex=2,gsubset=NULL,names=NULL){
	if(is.null(gsubset))
		stop("Must specify a gate to plot")
	l<-length(gs)
	rs<-floor(sqrt(l))
  	cs<-l%/%rs
  	while(rs*cs<l){
    	if((rs+1)*cs<rs*(cs+1))
      		rs<-rs+1
    	else
      		cs<-cs+1
  	}
  	grid<-expand.grid(1:cs,1:rs)
	g2<-NULL;
	if(length(gsubset)==2)
		g2<-gsubset[2]
	g1<-gsubset[1]
	if(is.null(names)){
		names<-flowWorkspace:::getSamples(gs)
	}else if(length(names)!=length(gs)){
		stop("names must be same length as gating set")
	}
	for(i in 1:(length(gs)-1)){
		nodes<-flowWorkspace:::getNodes(gs[[i]])
		print(flowWorkspace:::plotGate(gs[[i]],nodes[g1],lwd=2,pch='.',cex=cex,main=names[i]),split=c(grid[i,1],grid[i,2],cs,rs),more=T)
		if(!is.null(g2)){
			flowWorkspace:::plotGate(gs[[i]],nodes[g2],lwd=2,add=T)
		}
	}
	nodes<-flowWorkspace:::getNodes(gs[[length(gs)]])
	print(flowWorkspace:::plotGate(gs[[length(gs)]],nodes[g1],lwd=2,pch='.',cex=cex,main=names[length(gs)]),split=c(grid[length(gs),1],grid[length(gs),2],cs,rs),more=F)
	if(!is.null(g2)){
		flowWorkspace:::plotGate(gs[[length(gs)]],nodes[g2],lwd=2,add=T)
	}
}
#updateFlowFrameRange<-function(x){
#	flowCore:::checkClass(x,"GatingSet")
#	nc<-flowWorkspace:::getNcdf(x[[1]]);
#	s<-sampleNames(nc)
#	sapply(s,function(ss){
#		myframe<-nc[[ss]]
#		r<-apply(exprs(myframe),2,range)
#		pars<-myframe@parameters
#		data<-pars@data
#		data$minRange<-
#	})
#}
comparativeNormalizationPlot<-function(x,y,g,s,g2=NULL){
	flowCore:::checkClass(x,"GatingSet")
	flowCore:::checkClass(y,"GatingSet")
	flowCore:::checkClass(g,"numeric")
	flowCore:::checkClass(s,"numeric")
	if(g>length(flowWorkspace::getNodes(x[[1]]))){
		stop("Gate index out of bounds")
	}
	if(s>length(x)){
		stop("Sample index out of bounds")
	}
	if(length(x)!=length(y)){
		stop("The two gating sets should be of the same size")
	}
	if(!is.null(g2)){
		G2<-flowWorkspace::getGate(y[[s]],g2)@boundaries
	}else{
		G2<-NULL
	}
	par<-flowWorkspace::getParent(x[[1]],g)
	dims<-flowWorkspace::getDimensions(x[[1]],flowWorkspace::getNodes(x[[1]])[g])
	form<-sapply(dims,function(f)as.formula(paste("~`",f,"`",sep="")))
	print(densityplot(form[[1]],flowWorkspace::getData(x,par),main="Raw"),split=c(1,1,3,2),more=TRUE)
	print(densityplot(form[[1]],flowWorkspace::getData(y,par),main="Normalized"),split=c(1,2,3,2),more=TRUE)
	
	print(densityplot(form[[2]],flowWorkspace::getData(x,par),main="Raw"),split=c(2,1,3,2),more=TRUE)
	print(densityplot(form[[2]],flowWorkspace::getData(y,par),main="Normalized"),split=c(2,2,3,2),more=TRUE)
	
	print(flowWorkspace::plotGate(x[[s]],g,main="Raw"),split=c(3,1,3,2),more=TRUE)
	if(!is.null(G2)){
	trellis.focus(highlight=FALSE)
	panel.polygon(G2,border="red")
	trellis.unfocus()
	}
	print(flowWorkspace::plotGate(y[[s]],g,main="Normalized"),split=c(3,2,3,2),more=FALSE)
	if(!is.null(G2)){
	trellis.focus()
	panel.polygon(G2,border="red")
	trellis.unfocus()
	}
}




setMethod("normalize",c("GatingSetInternal","missing"),function(data,x="missing",...){
			.normalizeGatingSetInternal(x=data,...)			
		})
setMethod("normalize",c("GatingSet","missing"),function(data,x="missing",...){
			.normalizeGatingSet(x=data,...)			
		})

.normalizeGatingSet<-function(x,target=NULL,skipgates=NULL,skipdims=c("FSC-A","SSC-A","FSC-H","SSC-H","Time"),subsample=NULL,chunksize=10,nPeaks=list(),bwFac=2,...){
	samples<-flowWorkspace:::getSamples(x)
	valid<-target%in%samples
	if(!is.null(target)){
		if(!valid){
			stop("target ",target," not in the GatingSet")
		}
	}

	#Get all the non-boolean gates, breadth first traversal
	#Do a breadth-first traversal
#	browser()
	bfsgates<-lapply(x,function(y)
							RBGL::bfs(y@tree)[which(sapply(RBGL::bfs(y@tree)
																	,function(x)
																	!flowWorkspace:::.isBooleanGate.graphNEL(y,x)
														)
													)
												]
						)
	gates<-lapply(x,function(y)getNodes(y))
	
	#reorder the bfsgates so that indices match the order in getNodes
	for(i in seq_along(bfsgates)){
		bfsgates[[i]]<-match(bfsgates[[i]],gates[[i]])
	}
	bfsgates<-unique(do.call(rbind,bfsgates))
	np<-vector("list",max(bfsgates))
	for(p in seq_along(nPeaks)){
		np[[p]]<-nPeaks[[p]]
	}
	#gate-specific channel list to track normalization
	parentgates<-list();	
	
	#Keep a vector of normalized dimensions
	unnormalized<-NULL

	#Initialize master channel list
	unnormalized<-colnames(getData(x[[1]]))

		#Create a new gating set for the normalized data.
		#So we need a new data environment in a new gatingset object
		#x will be our new gatingset object and will be returned as a copy
		#which environments need to be reassigned?? all of them?
		#Copy the metadata, group, and data environments to new environments
		#Then update the data.
		#For ALL gates.. I should only need to do this once.
		#should cloning the gating set copy the ncdf files by default?
	if(flowWorkspace:::isNcdf(x[[1]]))
		x<-flowWorkspace:::cloneGatingSet(x,clone.data=TRUE,clone.gating=TRUE)	

	#for each gate, grab the dimensions and check if they are normalized.
	#Normalize what hasn't been normalized yet, then do the gating.
	#Set a target sample by name
	for(i in bfsgates){
		gate<-flowWorkspace:::getGate(x[[1]],i)

		#Root node
		if(class(gate)=="logical"){ #implies gate is NA, since NA is class logical, otherwise some flowCore gate class.
			next;
			}else{
				#check which dimensions are normalized already
				dims<-flowWorkspace:::getDimensions(x[[1]],flowWorkspace:::getNodes(x[[1]])[i])
				dims<-setdiff(dims,skipdims)
				#Get the PARENT gate (since we'll be gating the data using gate i)
				g<-flowWorkspace:::getParent(x[[1]],i)
				#initialize gate-specific normalization list
				if(is.null(parentgates[[as.character(g)]])){
					parentgates[[as.character(g)]]<-list();
					parentgates[[as.character(g)]]$unnormalized<-unnormalized
					parentgates[[as.character(g)]]$normalized<-NULL
				}
				#Data will be subset at gate g (parent) and normalized on dims of i (child)
				#keep a list of normalized and unnormalized channels for each parent gate..
				#Check the gate being normalized.. make sure 74 is done correctly
				if(!i%in%skipgates){
					wh.dim<-dims%in%parentgates[[as.character(g)]]$unnormalized
					parentgates[[as.character(g)]]$normalized<-c(parentgates[[as.character(g)]]$normalized,dims[wh.dim]);
					parentgates[[as.character(g)]]$unnormalized<-setdiff(parentgates[[as.character(g)]]$unnormalized,dims)
					stains<-dims[wh.dim]
					if(length(stains)!=0&gateHasSufficientData(x,g,...)){
						npks<-np[[i]]
						#TODO code for non ncdfFlowSet
						result<-warpSetGS(x,stains=stains,gate=g,target=target,subsample=subsample,chunksize=chunksize,peakNr=npks,bwFac=bwFac,...)
						if(flowWorkspace:::isNcdf(x[[1]])){
					 		sapply(sampleNames(result),function(s)ncdfFlow:::updateIndices(result,s,NA))
					
							#assign result to the Gating set data environment
							lfile<-flowWorkspace:::getNcdf(x)@file
							lapply(x,function(gh)assign("ncfs",result,graph:::nodeDataDefaults(gh@tree,"data")[["data"]]))		
						}else{
							browser()
							#TODO assign flowset to gatingset.
							for(j in seq_along(x)){
								odat<-getData(x[[j]])
								inds<-flowWorkspace:::getIndices(x[[j]],getNodes(x[[j]])[g])
								odat@exprs[inds,]<-result[[j]]@exprs
								assign("data",odat,graph:::nodeDataDefaults(x[[j]]@tree,"data"))
							}
						}
					}
						flowWorkspace:::recomputeGate(x,i);			
						for( j in seq_along(x)){
							flowWorkspace:::writeGatesToNetCDF(x[[j]])
						}		
				}
			}
			
			#flowWorkspace:::writeGatingSetGatesToNetCDFParallel(x,isNew=FALSE)
		}
		x
}

.normalizeGatingSetInternal<-function(x,target=NULL,skipgates=NULL,skipdims=c("FSC-A","SSC-A","FSC-H","SSC-H","Time"),subsample=NULL,chunksize=10,nPeaks=list(),bwFac=2,...){

#	browser()
	samples<-getSamples(x)
	valid<-target%in%samples
	if(!is.null(target)){
		if(!valid){
			stop("target ",target," not in the GatingSet")
		}
	}
	
	#Get all the non-boolean gates, breadth first traversal
	#Do a breadth-first traversal
	nodelist<-getNodes(x[[1]],order="bfs",prefix=TRUE)
	
	bfsgates<-unlist(lapply(nodelist,function(curNode){
#							browser()
				if(curNode=="root")
					return(0)
				else
				{
					curNodeInd<-as.integer(strsplit(curNode,split="\\.")[[1]][1])+1
					if(class(getGate(x[[1]],curNodeInd))!="BooleanGate")
						return(curNodeInd)
				}
			}))
	

	np<-vector("list",length(bfsgates))
	
	names(np)<-bfsgates
	for(p in seq_along(nPeaks)){
		np[[p]]<-nPeaks[[p]]
	}
	#gate-specific channel list to track normalization
	parentgates<-list();	
	
	#Keep a vector of normalized dimensions
	unnormalized<-NULL
	
	#Initialize master channel list
	unnormalized<-colnames(getData(x[[1]]))
	
	message("cloning the getingSet...")
		
	x<-clone(x)	
#	browser()
	#for each gate, grab the dimensions and check if they are normalized.
	#Normalize what hasn't been normalized yet, then do the gating.
	#Set a target sample by name
	for(i in bfsgates){
		
		gate<-getGate(x[[1]],i)
		
		#Root node
		if(class(gate)=="logical"){ #implies gate is NA, since NA is class logical, otherwise some flowCore gate class.
			next;
		}else{
			#check which dimensions are normalized already
			dims<-parameters(gate)
			dims<-setdiff(dims,skipdims)
			if(length(dims)>0)
			{
				
				#Data will be subset at gate g (parent) and normalized on dims of i (child)
				#keep a list of normalized and unnormalized channels for each parent gate..
				#Check the gate being normalized.. make sure 74 is done correctly
				if(!i%in%skipgates){
#					browser()
					#Get the PARENT gate (since we'll be gating the data using gate i)
					g<-getParent(x[[1]],i)
					#initialize gate-specific normalization list
					if(is.null(parentgates[[as.character(g)]])){
						parentgates[[as.character(g)]]<-list();
						parentgates[[as.character(g)]]$unnormalized<-unnormalized
						parentgates[[as.character(g)]]$normalized<-NULL
					}
					
					wh.dim<-dims%in%parentgates[[as.character(g)]]$unnormalized
					parentgates[[as.character(g)]]$normalized<-c(parentgates[[as.character(g)]]$normalized,dims[wh.dim]);
					parentgates[[as.character(g)]]$unnormalized<-setdiff(parentgates[[as.character(g)]]$unnormalized,dims)
					stains<-dims[wh.dim]
#					browser()	
					if(length(stains)!=0&gateHasSufficientData(x,g,...)){
						#choose the np element by name
						npks<-np[[as.character(i)]]

						result<-warpSetGS(x,stains=stains,gate=g,target=target,subsample=subsample,chunksize=chunksize,peakNr=npks,bwFac=bwFac,...)
											
						if(flowWorkspace:::isNcdf(x[[1]])){
							sapply(sampleNames(result),function(s)ncdfFlow:::updateIndices(result,s,NA))
							ncFlowSet(x)<-result
						}else{
							oldfs<-ncFlowSet(x)
							for(j in getSamples(x)){
								inds<-flowWorkspace::getIndices(x[[j]],g)
								oldfs[[j]]@exprs[inds,]<-result[[j]]@exprs
							}
							ncFlowSet(x)<-oldfs
						}
						recompute(x,i);	
					}
								
		
				}
			}
			
		}
		
	}
	x
}

#Function to check if each flowFrame at the given gate has enough data to normalize. 
#Set the threshold at 500 events for each flowFrame.
gateHasSufficientData<-function(x=NULL,g=NULL,minCountThreshold=500,...){
	#x could be flowSet or ncdfFlowSet
	res<-unlist(lapply(x,function(x)nrow(getData(x,g))>=minCountThreshold))
	if(all(res))
		return(TRUE)
	else
	{
		warning("not enough events to normalize: ",names(res[!res]))
		return(FALSE)
	}
}
#		hierarchy
#})
