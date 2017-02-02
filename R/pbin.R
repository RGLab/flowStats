
proBin<-function(m,minEvents=500,channels=NULL)
{ 
  if(!is(m,"flowFrame"))
	stop("m needs to be an object of class flowFrame")
	
  timeCol <- flowCore:::findTimeChannel(m,strict=FALSE)  
  cNames <-colnames(m) 
  if(length(timeCol)==0){
		m <- exprs(m)
  }else{
	    m <-exprs(m)[, cNames[!(cNames==timeCol)]]
  }
  if(is.null(channels)){
		channels <- colnames(m)
  }

  nodeIndx=1
  node =data.frame(dataIndx=1,visited=FALSE,parent=0,left=0,right=0,stringsAsFactors=FALSE)
  data=new.env(parent=emptyenv())
  data[[as.character(1)]]=m
    
  axisLimits=list()
  axisLimits[[1]]=data.frame(axisMin=rep(-Inf,ncol(m)),axisMax=rep(Inf,ncol(m)))
    
  splitPars=list()
  splitPars[[1]]=data.frame(splitCol=0,splitMed=0)		
    
  while (1)
  {
    if( node$left[nodeIndx] !=0 && node$visited[node$left[nodeIndx]]==FALSE)
    {
      nodeIndx= node$left[nodeIndx]
    
    }
    else if( node$right[nodeIndx] !=0 && node$visited[node$right[nodeIndx]]==FALSE)
    {
      nodeIndx= node$right[nodeIndx]
    
    }
    else if(node$visited[nodeIndx]==FALSE)
    { 
      node$visited[nodeIndx]=TRUE
      
      if(nrow(data[[as.character(nodeIndx)]])>minEvents)
      {
		y<-apply(data[[as.character(nodeIndx)]][,channels,drop=FALSE],2,var)
		tName <- which.max(y)
		t <- which(colnames(data[[as.character(nodeIndx)]])== attr(tName,"names"))
        #t<-which.max(y)
        x<-median(data[[as.character(nodeIndx)]][,t])
        s1=data[[as.character(nodeIndx)]][,t] > x
        if(!(all(s1) || all(!s1)))
        { 
       
        m1<-data[[as.character(nodeIndx)]][s1,,drop=FALSE]      
        m2<-data[[as.character(nodeIndx)]][!s1,,drop=FALSE]  
	        
	temp=nrow(node)+1; 
        node$left[nodeIndx]=temp;
        node$right[nodeIndx]=temp+1;
        data[[as.character(temp)]]=m1;
        data[[as.character(temp+1)]]=m2;
        splitPars[[temp]]=data.frame(splitCol=t,splitMed=x)
	splitPars[[temp+1]]=data.frame(splitCol=t,splitMed=x)
		
        axisLimits[[temp]]=axisLimits[[nodeIndx]]
        axisLimits[[temp+1]]=axisLimits[[nodeIndx]]
        axisLimits[[temp]]$axisMin[t]=x
        axisLimits[[temp+1]]$axisMax[t]=x

        node[temp,]=list(dataIndx=as.numeric(temp),visited=FALSE,parent=nodeIndx,left=0,right=0)
        node[temp+1,]=list(dataIndx=as.numeric(temp+1),visited=FALSE,parent=nodeIndx,left=0,right=0)
      }
       }
    }
    else
    {
      nodeIndx=node$parent[nodeIndx]
      if(node$parent[nodeIndx]==0)
      {
        if(node$visited[node$right[nodeIndx]]==TRUE)
        {
          break;
        }
      
      }
    }
          
  }
  return(list(table=node, data=data,limits=axisLimits,splitPars=splitPars))
}


binByRef<-function(binRes,data)  #output of proBin function, flowFrame
{
   data<-exprs(data)
   lmts=binRes$table$dataIndx[binRes$table$left==0]
   storeIndx=1;
   binned=new.env(parent=emptyenv())
   for (i in lmts)
   { 
      cStoreIndx<-as.character(storeIndx)
      len=nrow(binRes$limits[[i]])
      binned[[cStoreIndx]]<-data
      for(j in seq_len(len))
      { 
	indx<-!is.na(cut( binned[[cStoreIndx]][,j],binRes$limits[[i]][j,],labels=FALSE))
	binned[[cStoreIndx]]<- matrix(binned[[cStoreIndx]][indx],ncol=ncol(binned[[cStoreIndx]]))
      }
      storeIndx<-storeIndx+1
   } 
   return(binned)
   
}
   
   
plotBins<-function(binRes,data,channels=c("FSC-H","SSC-H"),title="",residuals=NULL,shadeFactor=0.6) 
{
   if(!all(channels %in% colnames(data)))
   stop("Invalid channel(s)")
  
  timeCol <- flowCore:::findTimeChannel(data,strict=FALSE)
  cNames <-colnames(data) 
  if(length(timeCol)==0){
		data <- exprs(data)
  }else{
		data<- exprs(data)[, cNames[!(cNames==timeCol)]]
  }
 
   plotColX<-which(colnames(data)==channels[1])
   plotColY<-which(colnames(data)==channels[2])

   plot(x=data[,plotColX], y=data[,plotColY],
   xlab=channels[1],  ylab=channels[2],
   main=title, col="#0080ff30")

   fill<-"transparent"
   if(!is.null(residuals))
   {
   res <- (((residuals - min(residuals)) / max(residuals-min(residuals)))*shadeFactor)
   shade  <- rgb(0,0,0, res)
   }

   x.limits=range(data[,plotColX])
   y.limits=range(data[,plotColY])
   limits=binRes$limits[binRes$table$left==0]
   pars=binRes$splitPars[binRes$table$left==0]
        
   for(i in seq_along(limits))
   {
      x1=limits[[i]]$axisMin[plotColX]
      x2=limits[[i]]$axisMax[plotColX]
      y1=limits[[i]]$axisMin[plotColY]
      y2=limits[[i]]$axisMax[plotColY]
	  
      if(pars[[i]]$splitCol==plotColX)
      {
      ##plot vertical lines
      y1=max(limits[[i]]$axisMin[plotColY], y.limits[1])    
      y2=min(limits[[i]]$axisMax[plotColY], y.limits[2])
      x1=max(limits[[i]]$axisMin[plotColX], x.limits[1])    
      x2=min(limits[[i]]$axisMax[plotColX], x.limits[2])
      }
      else if(pars[[i]]$splitCol==plotColY)  
      {
      x1=max(limits[[i]]$axisMin[plotColX], x.limits[1])   
      x2=min(limits[[i]]$axisMax[plotColX], x.limits[2])
      y1=max(limits[[i]]$axisMin[plotColY], y.limits[1])    
      y2=min(limits[[i]]$axisMax[plotColY], y.limits[2])
      }
            
      if(!is.null(residuals))
      {
	fill<-shade[i]
      }
  
      rect(x1,y1,x2,y2,border="red", col=fill)
   }

}
   

calcPearsonChi<-function(ctrlRes,sampRes)
{  
      binId<-ctrlRes$table$dataIndx[ctrlRes$table$left==0]
    binCount<-length(binId)
    input=matrix(0,2,binCount)
    for( i in seq_len(binCount))
    { 
        input[1,i]=nrow(ctrlRes$data[[as.character(binId[i])]])     #expected
        input[2,i]=nrow(sampRes[[as.character(i)]])  #observed
    }
    chisq.test(input)
}

calcPBChiSquare<-function(ctrlRes,sampRes,ctrlCount,sampCount)
{
  result<-list()
  binId<-ctrlRes$table$dataIndx[ctrlRes$table$left==0]
  binCount<-length(binId)
  for( i in seq_len(binCount))
  { 
    ctrl= (nrow(ctrlRes$data[[as.character(binId[i])]]))/ctrlCount
    samp=(nrow(sampRes[[as.character(i)]]))/sampCount
    result$chiSq[i]<-( (samp-ctrl)^2)/(ctrl+samp)
    result$expected[i]<-ctrl*ctrlCount
    result$observed[i]<-samp*sampCount
    result$residuals[i]<-(samp-ctrl)/(sqrt(ctrl+samp))
  }
  result$pbStat<-(2*ctrlCount*sampCount*(sum(result$chiSq))/(ctrlCount+sampCount)
     -(binCount-1))/(sqrt(2*(binCount-1)))
    
  return(result)
  }
  
