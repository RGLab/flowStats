#Sequentially normalize a gating set by normalizing channels at each step of gating
#return a normalized gatingset
#target is the target sample to normalize against
#skip is a list of gate indices for which you want to skip normalization (i.e. dump gates that may be difficult to normalize, or FSC/SSC)
#getNormGateList<-function(x){
#	flowCore:::checkClass(x,"GatingSet")
#	bfsgates<-lapply(x,function(y)which(sapply(RBGL::bfs(y@tree),function(x)!flowWorkspace:::.isBooleanGate.graphNEL(y,x))))
#	bfsgates<-unique(do.call(rbind,bfsgates))
#	return(data.frame(gate=flowWorkspace::getNodes(x[[1]], showHidden = TRUE)[bfsgates],index=t(bfsgates)))
#}
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
    print(flowWorkspace::plotGate(gh,nodes[g],lwd=2,pch='.',cex=cex,main=strsplit(nodes[g],"\\.")[[1]][2]),split=c(grid[g,1],grid[g,2],rs,cs),more=T)
  })
  print(flowWorkspace::plotGate(gh,nodes[length(nodes)],lwd=2,pch='.',cex=cex,main=strsplit(nodes[length(nodes)],"\\.")[[1]][2]),split=c(grid[length(nodes),1],grid[length(nodes),2],rs,cs),more=F)
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
		names<-flowWorkspace::sampleNames(gs)
	}else if(length(names)!=length(gs)){
		stop("names must be same length as gating set")
	}
	for(i in 1:(length(gs)-1)){
		nodes<-flowWorkspace::getNodes(gs[[i]], showHidden = TRUE)
		print(flowWorkspace::plotGate(gs[[i]],nodes[g1],lwd=2,pch='.',cex=cex,main=names[i]),split=c(grid[i,1],grid[i,2],cs,rs),more=T)
		if(!is.null(g2)){
			flowWorkspace::plotGate(gs[[i]],nodes[g2],lwd=2,add=T)
		}
	}
	nodes<-flowWorkspace::getNodes(gs[[length(gs)]], showHidden = TRUE)
	print(flowWorkspace::plotGate(gs[[length(gs)]],nodes[g1],lwd=2,pch='.',cex=cex,main=names[length(gs)]),split=c(grid[length(gs),1],grid[length(gs),2],cs,rs),more=F)
	if(!is.null(g2)){
		flowWorkspace::plotGate(gs[[length(gs)]],nodes[g2],lwd=2,add=T)
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
	if(g>length(flowWorkspace::getNodes(x[[1]], showHidden = TRUE))){
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
	dims<-flowWorkspace::getDimensions(x[[1]],flowWorkspace::getNodes(x[[1]], showHidden = TRUE)[g])
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




setMethod("normalize",c("GatingSet","missing"),function(data,x="missing",...){
			.normalizeGatingSet(x=data,...)			
		})


.normalizeGatingSet <- function(x,target=NULL
                                    , populations = NULL
                                    , dims
                                    ,nPeaks=list(),ncdfFile = NULL, minCountThreshold = 500, ...){

#	browser()
	samples<-sampleNames(x)
	valid<-target%in%samples
	if(!is.null(target)){
		if(!valid){
			stop("target ",target," not in the GatingSet")
		}
	}
	
	#Get all the non-boolean gates, breadth first traversal
	#Do a breadth-first traversal
	
    nodes <- getNodes(x[[1]],order="bfs",showHidden = TRUE, path = "auto")
	
	#gate-specific channel list to track normalization
	parentgates<-list();	
	
	#Initialize master channel list
	unnormalized <- colnames(getData(x[[1]]))
	
    
	message("cloning the gatingSet...")
		
	x <- clone(x, ncdfFile = ncdfFile)	
	
	#for each gate, grab the dimensions and check if they are normalized.
	#Normalize what hasn't been normalized yet, then do the gating.
	#Set a target sample by name
	for(node in populations){
		
        
        if(node!="root")
          if(!flowWorkspace:::.isBoolGate(x[[1]],node))
        {
            
    		gate <- getGate(x[[1]],node)
    		
    		dims <- intersect(parameters(gate), dims)
            
        
    		if(length(dims)>0)
    		{
    			
    			#Data will be subset at gate g (parent) and normalized on dims of i (child)
    			#keep a list of normalized and unnormalized channels for each parent gate..
    			#Check the gate being normalized.. make sure 74 is done correctly
    			 
                message("Normalize ", node)
  					
  				#Get the PARENT gate (since we'll be gating the data using gate i)
  				parent <- getParent(x[[1]], node)
  				#initialize gate-specific normalization list
  				if(is.null(parentgates[[as.character(parent)]])){
  					parentgates[[as.character(parent)]]<-list();
  					parentgates[[as.character(parent)]]$unnormalized<-unnormalized
  					parentgates[[as.character(parent)]]$normalized<-NULL
  				}
                  #check which dimensions are normalized already
  				wh.dim<-dims%in%parentgates[[as.character(parent)]]$unnormalized
  				parentgates[[as.character(parent)]]$normalized<-c(parentgates[[as.character(parent)]]$normalized,dims[wh.dim]);
  				parentgates[[as.character(parent)]]$unnormalized<-setdiff(parentgates[[as.character(parent)]]$unnormalized,dims)
  				stains<-dims[wh.dim]
  #					browser()	
  				if(length(stains)!=0&&gateHasSufficientData(x, parent, minCountThreshold = minCountThreshold, ...))
                {
  					#choose the np element by name
				  npks <- nPeaks[[node]]

  					result <- warpSetGS(x,stains = stains
                                              ,node = parent
                                              ,target = target
                                              ,peakNr = npks
                                              ,...)
  										
  					if(flowWorkspace::isNcdf(x)){
  						sapply(sampleNames(result),function(s)ncdfFlow::updateIndices(result,s,NA))
  						flowData(x)<-result
  					}else{
  						oldfs<-flowData(x)
  						for(j in sampleNames(x)){
  							inds<-flowWorkspace::getIndices(x[[j]],parent)
  							oldfs[[j]]@exprs[inds,]<-result[[j]]@exprs
  						}
  						flowData(x)<-oldfs
  					}
  					recompute(x,node);	
  				}
    							
    			
    		}
          }
	}
	x
}

#Function to check if each sample at the given gate has enough data to normalize. 
#Set the threshold at 500 events for each flowFrame.
gateHasSufficientData<-function(x=NULL,g=NULL,minCountThreshold=500,...){
  
	res <- unlist(flowWorkspace::lapply(x,function(y){
              nrow(getData(y,g))>=minCountThreshold
                  }))
	if(all(res))
		return(TRUE)
	else
	{
		warning("not enough events to normalize: ",paste(names(res[!res]),collapse="\n"))
		return(FALSE)
	}
}
#		hierarchy
#})
