#' Generate a prior specification based on a flowClust model This function
#' generates a prior specification based on a flowClust fit object It can be
#' passed to a second round of flowClust() with usePrior="yes" The prior could
#' be estimated from a single sample, for example, and then used to speed up
#' the convergence for other samples.
#' 
#' Generate a prior specification based on a flowClust model This function
#' generates a prior specification based on a flowClust fit object It can be
#' passed to a second round of flowClust() with usePrior="yes" The prior could
#' be estimated from a single sample, for example, and then used to speed up
#' the convergence for other samples.
#' 
#' 
#' @param x a flowClust fit object
#' @param kappa is the fraction of equivalent observations by which to weight
#' this prior relative to the flowClust model.
#' @param Nt the number of total equivalent observation
#' @param addCluster not currently supported
#' @export 
flowClust2Prior<-function(x,kappa,Nt=NULL,addCluster=NULL){
  if(is.null(Nt)){
    Nt<-nrow(x@z)
  }
    p<-ncol(x@mu)
    K<-x@K
    
    nu0<-Ng<-x@w*Nt
    if(all((nu0*kappa-p-1)>0)){
        Lambda0<-x@sigma;
        for(i in 1:K){
            Lambda0[i,,]<-Lambda0[i,,]*(kappa*nu0[i]-p-1)
        }
    }else{
        stop("Can't proceed. Prior nu0 is negative for cluster(s) ",paste(which((nu0-p-1)>0),collapse=","),"\n(p-1) = ",p-1,": Try increasing kappa")
    }
    
    Omega0<-array(0,c(K,p,p))
    for(i in 1:K){
        #Make precision of prior means spherical.
        Omega0[i,,]<-diag(1,p)
        if(p==1){
            dS<-x@sigma[i,,]
            dO<-Omega0[i,,]
        }else{
            dS<-det(x@sigma[i,,])
            dO<-det(Omega0[i,,])
        }
        k<-(dO/dS)^(1/p)
        #rescale Omega0
        Omega0[i,,]<-Omega0[i,,]*k
        #Multiply by Ng*kappa
        Omega0[i,,]<-solve(Omega0[i,,]*Ng[i]*kappa)
    }
    nu0<-nu0*kappa
    Mu0<-x@mu
    
    lambda<-x@lambda
    #w0 are dirichlet prior parameters
    #Should make them vague
    w0<-x@w*Nt
    
    if(!is.null(addCluster)){
      for(i in (K+1):(K+addCluster)){
           S<-stats::cov(Mu0)
           Lam<-array(0,c(K+1,p,p))
           om<-array(0,c(K+1,p,p))
           Mu0<-rbind(Mu0,colMeans(Mu0))
           for(i in 1:K){
               om[i,,]<-Omega0[i,,]
               Lam[i,,]<-Lambda0[i,,]
           }
           om[K+1,,]<-diag(1,p)
           Lam[K+1,,]<-diag(1,p)
           diag(Lam[K+1,,])<-diag(S)
           diag(om[K+1,,])<-diag(S)
           if(p==1){
               dS<-Lam[K+1,,]
               dO<-om[K+1,,]
           }else{
               dS<-det(Lam[K+1,,])
               dO<-det(om[K+1,,])
           }
           k<-(dO/dS)^(1/p)
           om[K+1,,]<-om[K+1,,]*k
           #om[K+1,,]<-om[K+1,,]
           Omega0<-om
           Lambda0<-Lam
           nu0<-c(nu0,p+2)
           w0<-c(w0,1)
           K<-K+1
      }
    }
    prior<-list(Mu0=Mu0,Lambda0=Lambda0,Omega0=Omega0,w0=w0,nu0=nu0,nu=x@nu,lambda=x@lambda,K=K)
    class(prior)<-"flowClustPrior"
    attr(prior,"lambda")<-x@lambda
    prior;
}

#' Generate a flowClust prior specification
#' 
#' Generate a flowClust prior specification from gates and data
#' 
#' Construct a prior specification. Generally not called by the user.
#' 
#' @aliases mkPrior mkPrior,list,flowSet,missing,missing-method
#' @param gate A list of flowCore gates. The gates should represent the SAME
#' population gated across multiple samples.
#' @param data A flowSet of the same size as the number of gates above. Each
#' flowFrame in the flowSet should contain the events representing the
#' population in its corresponding gate. i.e. it should be the gated data.
#' @param nu0 The nu0 hyperparameter. For estimation from data, it should be
#' nu0=NA.
#' @param Omega0 The Omega0 hyperparameter. For estimation from data it can be
#' missing.
#' @param \dots Not currently used.
#' @return Return values depend on the specific method called. Not meant for
#' user consumption.
#' @author Greg Finak \email{gfinak@@fhcrc.org}
#' @references \url{http://www.rglab.org}
#' @examples
#' 
#
#' ## The function is currently defined as
#' @rdname mkPrior
#' @export 
setGeneric("mkPrior",function(gate,data,nu0,Omega0,...){
standardGeneric("mkPrior");	
})
# ===============================================================================
# = Generate a prior from a polygonGate and a flowFrame. Provide nu0 and Omega0 =
# ===============================================================================
#' @rdname mkPrior
#' @export 
setMethod("mkPrior",signature("polygonGate","flowFrame","numeric","matrix")
                ,function(gate,data,nu0,Omega0){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);	
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-stats::cov(data)*(nu0-length(dims)-1)
	Mu0<-colMeans(data)
	if(length(dims)!=ncol(Omega0)){
		stop("Invalid dimensions of \"Omega0\". Received ",dim(Omega0)[1]," x ",dim(Omega0)[2]," but expecting ",length(dims)," x ",length(dims));
	}
	prior=list(Mu0=Mu0,Lambda0=Lambda0,nu0=nu0,Omega0=Omega0)
	#message("Making prior from gate and data");
	prior;
})
#hyperparemeters not specified. Default nu0=500, Omega0=diag(1/1000,D)

# =========================================================
# = Nu0 and Omega0 specified. RectangleGate and flowFrame =
# =========================================================
#' @rdname mkPrior
setMethod("mkPrior",signature("rectangleGate","flowFrame","numeric","matrix")
              ,function(gate,data,nu0,Omega0){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-stats::cov(data)*(nu0-length(dims)-1)
	Mu0<-colMeans(data)
	if(length(dims)!=ncol(Omega0)){
		stop("Invalid dimensions of \"Omega0\". Received ",dim(Omega0)[1]," x ",dim(Omega0)[2]," but expecting ",length(dims)," x ",length(dims));
	}
	prior=list(Mu0=Mu0,Lambda0=Lambda0,nu0=nu0,Omega0=Omega0)
	prior;
})

# =========================================================================
# = hyperparemeters not specified. 
# = Returns an  incomplete prior specification. Shouldn't be called       =
# = Directly by the user.												  =
# =========================================================================
#' @rdname mkPrior
setMethod("mkPrior",signature("rectangleGate","flowFrame","missing","missing")
          ,function(gate,data,nu0=NA,Omega0=NA){
	gc(reset=TRUE)
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-stats::cov(data) #This is the covariance matrix, not the real Lambda0
	Mu0<-colMeans(data)
	n<-dim(data)[1];
	prior<-list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})
# ===================================================================================
# = Generate a prior from a polygonGate and a flowFrame with nu0 and Omega0 missing
# = Returns and incomplete prior specification. Should not be called by the user    =
# ===================================================================================
#' @rdname mkPrior
setMethod("mkPrior",signature("polygonGate","flowFrame","missing","missing")
        ,function(gate,data,nu0=NA,Omega0=NA){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-stats::cov(data)
	Mu0<-colMeans(data)
	n<-dim(data)[1]
	prior=list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})


# =========================================================================
# = Plot a prior given some data (a flowFrame) and a prior specification. =
# = The prior specification is for a single cluster. 					  =
# =========================================================================



#' Plots a flowClust prior over some data.
#' 
#' Plots a flowClust prior overlaid on data.
#' 
#' Generates a plot of a "flowClustPrior" or "flowClustPriorList" object
#' overlaid on some data. Plots the prior means (Mu0), prior covariance of the
#' means (Omega0), and prior sample covariance (Lambda0).
#' 
#' @param data On object of class "flowFrame". The data to be plotted.
#' @param prior An object of class "flowClustPrior", or "flowClustPriorList",
#' returned by a call to \code{mkPrior}.
#' @param dims A character vector of the dimensions to be included in the plot.
#' The dimension names should match column names in the prior and in the
#' flowFrame.
#' @param \dots Additional arguments to plotting functions, such as
#' \code{smooth=TRUE/FALSE}
#' @return Silently returns zero.
#' @author Greg Finak <gfinak@@fhcrc.org>
#' @keywords aplot dplot
#' @export 
plotPrior<-function(data,prior,dims=NULL,...){
	if(!"smooth"%in%names(list(...))){
		sm=FALSE
	}else
	  stop('smooth argument is no longer supported!')
	if(class(prior)=="flowClustPriorList"){
		prior<-prior[[1]];
	}
	if(!is.null(dims)){
		if(all(dims%in%colnames(prior$Mu0))){
			dims<-dims[na.omit(match(colnames(prior$Mu0),dims))]
			dim.inds<-match(dims,colnames(prior$Mu0))
		}else{
			stop("Can't find ",dims[which(!(dims%in%colnames(prior$Mu0)))]," in the prior.");
		}
	}else{
		dims<-colnames(prior$Mu0)
		dim.inds<-match(dims,colnames(prior$Mu0))
	}
	k<-nrow(prior$Mu0);
	nd<-ncol(prior$Mu0)
	if(nd>1){
		
			plot(exprs(data[,dims]), pch = ".", ...);
				
		for(i in 1:k){
			points(t(as.matrix(prior$Mu0[i,dim.inds])),pch=20,col="red")
			lines(ellipse(solve(prior$Omega0[i,dim.inds,dim.inds]),centre=prior$Mu0[i,dim.inds]),col="green",lwd=2,lty=2)
			lines(ellipse(prior$Lambda0[i,dim.inds,dim.inds]/(prior$nu0[i]-nd-1),centre=prior$Mu0[i,dim.inds]),col="red",lwd=2,lty=2)
		}
	}else{
		for(i in 1:k){
			
			hist(exprs(data[,dims]),breaks=256,...);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds])),col="red")
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.975)*sqrt(1/prior$Omega0[i,dim.inds,dim.inds]),col="red",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.025)*sqrt(1/prior$Omega0[i,dim.inds,dim.inds]),col="red",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.975)*sqrt((prior$Lambda0[i,dim.inds,dim.inds])/(prior$nu0-2)),col="green",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.025)*sqrt((prior$Lambda0[i,dim.inds,dim.inds])/(prior$nu0-2)),col="green",lty=2);
		}
	}
	invisible(0);
}
# ===========================================================================================
# = Pad the prior to K clusters. K must be greater than the number of clusters in the prior =
# ===========================================================================================
.padPriorToK<-function(prior,k,env){
	#env is an environment with the data to be used to pad the prior.
	flag<-0;
	if(inherits(prior,"flowClustPriorList")){
		prior<-prior[[1]];
		flag<-1;
	}else if(inherits(prior,"flowClustPrior")){
		prior<-prior;
	}else{
		stop("prior must be of class \"flowClustPrior\" or \"flowClustPriorList\"");
	}
	nd<-ncol(prior$Mu0);
	nc<-nrow(prior$Mu0);
	if(nc>=k){
		stop("Prior cannot be padded to fewer than or equal to ", k, " clusters. It already has ",nc," clusters.");
	}
	Mu0<-matrix(NA,k,nd);
	Omega0<-array(NA,c(k,nd,nd));
	Lambda0<-array(NA,c(k,nd,nd));
	nu0<-numeric(k);
	colnames(Mu0)<-colnames(prior$Mu0);
	for(i in 1:nc){
		Mu0[i,]<-prior$Mu0[i,];
		Omega0[i,,]<-prior$Omega0[i,,]
		Lambda0[i,,]<-prior$Lambda0[i,,]
		nu0[i]<-prior$nu0[i]
	}
	for(i in (nc+1):k){
		Mu0[i,]<-rep(0,nd)##Set to vague values
		Omega0[i,,]<-diag(0,nd)## Set to vague values
		Lambda0[i,,]<-diag(0,nd)## Set to vague values
		nu0[i]<-nd+1; #vague
	}
	pprior<-list(Mu0=Mu0,Omega0=Omega0,Lambda0=Lambda0,nu0=nu0);
	class(pprior)<-"flowClustPrior";
	if(flag){
		pprior<-list(pprior);
		class(pprior)<-"flowClustPriorList";
	}
	return(pprior);
}


# ============================================================
# = Get the parent node of a given node in a graphNEL object =
# ============================================================
.graph_parent<-function(g,n){
	return(setdiff(unlist(graph::adj(ugraph(g),n)),unlist(graph::adj(g,n))))
}

# ==================================================================================================
# = Estimate the hyperparameters of the prior given the means and covariances of multiple samples. =
# ==================================================================================================
#' @importFrom corpcor cov.shrink
#' @importFrom stats cov var optimize qf qchisq quantile kmeans optim mahalanobis dist  na.omit qnorm
#' @importFrom graphics hist curve stripchart lines points abline title contour image
#' @importFrom grDevices gray heat.colors rainbow terrain.colors topo.colors cm.colors
#' @noRd
.estimateHyperParameters<-function(priors,model.means,model.cov,nu0){
	#Empirical Bayes Estimation of hyperparameters (except nu0)
	d<-dim(priors[[1]]$Lambda0)[1]
	MuG0<-NA;LambdaG0<-NA;OmegaG0<-NA;nuG0<-NA;
	MuG0<-colMeans(do.call(rbind,lapply(priors,function(x)x$Mu0)))
	if(model.means=="full"){
		#What if n<p
		cm<-do.call(rbind,lapply(priors,function(x)x$Mu0))
		if(nrow(cm)==1){
			
			# TODO Test this more thoroughly, i.e. case where n=1, and using covariance of the means at 1% of the covariance of the data.
			OmegaG0<-solve(diag(diag(priors[[1]]$Lambda*0.01)))
		}else if(nrow(cm)<ncol(cm)){
			OmegaG0<-solve(corpcor::cov.shrink(cm,verbose=FALSE));
		}else{
			OmegaG0<-try(solve(stats::cov(cm)),silent=TRUE);
			if(inherits(OmegaG0,"try-error")){
				OmegaG0<-solve(corpcor::cov.shrink(cm,verbose=FALSE))
			}
		}
		#Sometimes the estimate is unstable and negative. In that case use a diagonal covariance parameterization
		if(OmegaG0[1,1]<0){
			OmegaG0<-solve(diag(apply(cm,1,var)))
		}
		
	
	}else if(model.means=="DU"){
		# TODO code to handle nrow = 1
		d<-dim(priors[[1]]$Lambda0)[1]
		OmegaG0<-solve(diag(diag(stats::cov(do.call(rbind,lapply(priors,function(x)x$Mu0)))),d))
#		if(OmegaG0[1,1]<0){
#		}
	}else if(model.means=="DE"){
		# TODO code to handle nrow=1
		d<-dim(priors[[1]]$Lambda0)[1]
		OmegaG0<-solve(diag((det(stats::cov(do.call(rbind,lapply(priors,function(x)x$Mu0)))))^(1/d),d))
	}
	nuG0<-nu0	
	if(model.cov=="full"){
		#This estimate of nu0 is not correct for this model, working on it.
		if(is.na(nu0)){
			X<-lapply(priors,function(x)x$Lambda0)
			d<-dim(priors[[1]]$Lambda0)[1]		
			nu0<-stats::optimize(f=.LIW,interval=c(d+2,20000),X=X)$minimum
			nuG0<-nu0
		}
		LambdaG0<-Reduce("+",lapply(priors,function(x)x$Lambda0))/length(priors)
		LambdaG0<-LambdaG0*(nu0) #ML estimate of LambdaG0
	}else if(model.cov=="DE"){
		#Estimate nu0
		if(is.na(nu0)){
			d<-dim(priors[[1]]$Lambda0)[1]
	nu0<-(solve(diag(diag(Reduce("+",lapply(priors,function(x)(x$Lambda0)))),d))*d)/(sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0)))))*(d-1))+diag(1/(d-1),d)
			nu0<-det(nu0)^(-1/d)
			nuG0<-nu0
		}
		#Or use the given value	
			d<-dim(priors[[1]]$Lambda0)[1]		
			LambdaG0<-diag((length(priors)*dim(priors[[1]]$Lambda0)[1]*nu0)/sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0))))),d)
	}else if(model.cov=="DU"){
		#This estimate of nu0 is not correct for this model.. working on it.
			if(is.na(nu0)){
				d<-dim(priors[[1]]$Lambda0)[1]	
				X<-lapply(priors,function(x)x$Lambda0)	#nu0<-(solve(diag(diag(Reduce("+",lapply(priors,function(x)(x$Lambda0)))),d))*d)/(sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0)))))*(d-1))+1/(d-1)
				# nu0<-det(nu0)^(-1/d)
				nu0<-stats::optimize(f=.LIW,interval=c(d+2,20000),X=X)$minimum
				nuG0<-nu0
			}
			d<-dim(priors[[1]]$Lambda0)[1]	
		LambdaG0<-diag(diag(Reduce("+",lapply(priors,function(x)x$Lambda0))/length(priors)),dim(priors[[1]]$Lambda0)[1])*(nuG0-d-1)
	}
	n<-list(unlist(lapply(priors,function(x)x$n)))# number of events in each sample. To be used to compute wi later.
	#include the number of events for dirichlet
	prior=list(Mu0=MuG0,Omega0=OmegaG0,Lambda0=LambdaG0,nu0=nuG0,n=n);
	# Mu0 must be a matrix for flowClust
	prior$Mu0<-t(as.matrix(prior$Mu0)) 
	#Omega0 and Lambda0 must be K,py,py arrays
	prior$Lambda0<-array(prior$Lambda0,c(1,ncol(prior$Lambda0),ncol(prior$Lambda0)))
	prior$Omega0<-array(prior$Omega0,c(1,ncol(prior$Omega0),ncol(prior$Omega0)))
	class(prior)<-"flowClustPrior";
	prior<-list(prior);
	class(prior)<-"flowClustPriorList";
	return(prior);
}
# ======================================================================================
# = We should also be able to generate a prior from a single gate and multiple samples =
# ======================================================================================
#' @rdname mkPrior 
setMethod("mkPrior",signature("list","flowSet",nu0="missing","missing")
      ,function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
 mkPrior(gate,data,nu0=NA)
 })

# ==============================
# = Prior from a flowSet alone =
# ==============================
#' @rdname mkPrior
#' @param model.cov,model.means model names used for cov and means. one of c("full","DE","DU").
#'                                          "full" is the default.
setMethod("mkPrior",signature("missing","flowSet",nu0="ANY","missing")
        ,function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
	priors<-list();
	method=match.arg(model.cov,c("full","DE","DU"));
	model=match(model.means,c("full","DU","DE"))
	repflag<-FALSE;
	estflag<-FALSE;
	if(!(is(nu0,"numeric")|is.na(nu0))){
		stop("nu0 must be a numeric vector of length >= 1, or NA for estimation");
	}
	if(length(nu0)>1){
		repflag<-TRUE;
	}
	cat(".")
	for(i in 1:length(data)){
		if(nrow(data[[i]])>5){
			priors[[i]]<-mkPrior(data=data[[i]])
		}else{
			#Ignore the gate if there's not enough sample points
			next;
		}
	}
	priors<-priors[unlist(lapply(priors,function(x)!is.null(x)))]
	if(length(priors)==0){
		#No gate has enough data to construct a prior.
		return(NA);
	}

	prior<-.estimateHyperParameters(priors,model.means=model.means,model.cov=model.cov,nu0=nu0);
	return(prior);
	
})
# ================================================================================
# = Construct a prior from a flowFrame alone. Not meant to be called by the user =
# ================================================================================
#' @rdname mkPrior
setMethod("mkPrior",signature("missing","flowFrame",nu0="missing","missing"),function(gate,data,nu0,Omega0){
	gc(reset=TRUE)
	##Use all the dimensions, since they're not specified.
	data<-exprs(data)
	if(ncol(data)>=nrow(data)){
		Lambda0<-corpcor::cov.shrink(data,verbose=FALSE)
		class(Lambda0)<-"matrix"
	}else{
		#Lambda0<-cov.shrink(data,verbose=FALSE) #This is the covariance matrix, not the real Lambda0
		#class(Lambda0)<-"matrix"
		Lambda0<-stats::cov(data) #This is the covariance matrix, not the real Lambda0
		
	}
	Mu0<-colMeans(data)
	n<-dim(data)[1];
	prior<-list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})
# ====================================================
# = Multiple gates (same gates) and multiple samples =
# = Calls mkPrior with the "missing" hyperparamter   =
# = Signatures										 =
# ====================================================
#' @rdname mkPrior
setMethod("mkPrior",signature("list","flowSet",nu0="ANY","missing"),function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
	priors<-list();
	method=match.arg(model.cov,c("full","DE","DU"));
	model=match(model.means,c("full","DU","DE"))
	repflag<-FALSE;
	estflag<-FALSE;
	if(!all(unlist(lapply(gate,function(x)class(x)=="polygonGate"|class(x)=="rectangleGate"),use.names=FALSE))){
		stop("All elements of \"gate\" must be class polygonGate or rectangleGate");
	}
	if(!(is(nu0,"numeric")|is.na(nu0))){
		stop("nu0 must be a numeric vector of length >= 1, or NA for estimation");
	}
	if(length(nu0)>1){
		repflag<-TRUE;
	}
	cat(".")
	#Deal only with gates where there are more than 5 samples
	for(i in 1:length(data)){
		sub<-Subset(data[[i]],filter(data[[i]],gate[[i]]))
		if(nrow(sub)>5){
			priors[[i]]<-mkPrior(gate[[i]],sub)
		}else{
			#Ignore the gate if there's not enough sample points
			next;
		}
	}
	priors<-priors[unlist(lapply(priors,function(x)!is.null(x)))]
	if(length(priors)==0){
		#No gate has enough data to construct a prior.
		return(NA);
	}

	prior<-.estimateHyperParameters(priors,model.means=model.means,model.cov=model.cov,nu0=nu0);
	return(prior);
})

# From MCMCpack - want to remove this dependency.
.ddirichlet <- function(x, alpha) 
{
	dirichlet1 <- function(x, alpha) {
		logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
		s <- sum((alpha - 1) * log(x))
		exp(sum(s) - logD)
	}
	if (!is.matrix(x)) 
		if (is.data.frame(x)) 
			x <- as.matrix(x)
		else x <- t(x)
		if (!is.matrix(alpha)) 
			alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
								 byrow = TRUE)
		if (any(dim(x) != dim(alpha))) 
			stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
		pd <- vector(length = nrow(x))
		for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
																				 ])
		pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
		pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
		return(pd)
}

# ====================================
# = Log-likelihood for the dirichlet =
# ====================================
.LD<-function(alpha,x){
	
	-sum(log(.ddirichlet(x,alpha)))
}

# =======================================================
# = Log-Likelihood for the Inverse Wishart Distribution =
# =======================================================
.LIW<-function (nu, X) 
{
    n <- length(X)
    d <- nrow(X[[1]])
    Sinv <- lapply(X, function(x) t(solve(x)))
    sumSinv <- Reduce("+", lapply(Sinv, function(x) t(x)))
    S = solve(sumSinv/(nu * n))
    A <- n * d * (d - 1)/4 * log(pi) - n * sum((sapply(1:d, function(j) lgamma(nu/2 + 
        (1 - j)/2))))
    B <- n * nu/2 * log(det(S))
    C <- -Reduce("+", lapply(X, function(x) log(det(x)))) * (nu - 
        d - 1)/2
    D <- -Reduce("+", lapply(Sinv, function(x) sum(diag(S %*% 
        x))/2))
    E <- -n * nu * d/2 * log(2)
    return(-sum(c(A, B, C, D, E)))
}
# ======================
# = Combine two priors =
# ======================
.mergePriors<-function(priors){
	l<-length(priors)
	nd<-dim(priors[[1]]$Mu0)[2]
	Omega0<-array(NA,c(l,nd,nd));
	Lambda0<-array(NA,c(l,nd,nd));
	nu0<-numeric(l);
	Mu0<-matrix(NA,l,nd);
	w0<-numeric(l);
	for(i in 1:l){
		Omega0[i,,]<-priors[[i]]$Omega0[1,,];
		Lambda0[i,,]<-priors[[i]]$Lambda0[1,,];
		Mu0[i,]<-priors[[i]]$Mu0[1,];
		nu0[i]<-priors[[i]]$nu0[1];	
	}
	#estimate w0 from dirichlet likelihood. Numerical optimization. Initialized at unity.
	x<-do.call(cbind,unlist(lapply(priors,function(x)x$n),recursive=F))
	for(i in 1:nrow(x)){
		x[i,]<-x[i,]/sum(x[i,])
	}
	wstar<-optim(par=rep(1,length(w0)),x=x,fn=.LD,method="L-BFGS-B",lower=rep(1,length(w0)))
	if(wstar$convergence!=0){
		w0<-rep(1,length(w0));
	}
	else{
		w0<-wstar$par;
	}
	
	colnames(Mu0)<-colnames(priors[[1]]$Mu0)
	mprior<-list(Mu0=Mu0,Omega0=Omega0,Lambda0=Lambda0,nu0=nu0,w0=w0);
	class(mprior)<-"flowClustPrior";
	return(mprior)
}
# ============================================================================================================
# = Determine if the priors for the children of a given node should be combined into a single specification. =
# ============================================================================================================
.combinePriorsForChildNodes<-function(x,n){
	children<-unlist(adj(x,n),use.names=F);
	priors<-nodeData(x,children,"prior");
	
	#Which children have a  prior specification?
	havePriors<-unlist(lapply(priors,function(x)inherits(x,"flowClustPrior")))
	#Compare only those with priors
	priors<-priors[havePriors];
	if(length(priors)>0){
	dims<-lapply(priors,function(x)colnames(x$Mu0))
	groups<-matrix(0,nrow=length(dims),ncol=length(dims))
	colnames(groups)<-names(dims);
	rownames(groups)<-names(dims);
	unassigned<-names(dims);
	assigned<-list();
	while(length(unassigned)!=0){
		i<-unassigned[1];
		k<-length(assigned)+1
		assigned[[k]]<-i;
		unassigned<-setdiff(unassigned,i)
		for(j in setdiff(unassigned,i)){
			if(identical(dims[[i]],dims[[j]])){
				assigned[[k]]<-c(assigned[[k]],j)
				unassigned<-setdiff(unassigned,j)
			}
		}		
	}
	mergedPriors<-list()
	#assigned is a list containing the names of identical nodes.
	## TODO fix this it doesn't appear to be correct ordering
	for(i in 1:length(assigned)){
		mergedPriors[[i]]<-.mergePriors(priors[assigned[[i]]])
		attr(mergedPriors[[i]],"popnames")<-assigned[[i]];
	}
	class(mergedPriors)<-"flowClustPriorList"
	}else{
		mergedPriors<-NA
	}
	mergedPriors;
}

# =========================================================================================================
# = Map function for flowClust model with bayes prior. Assigns events to clusters based on quantile. =
# = Assigns outliers conditional on cluster membership, not globally. Important for correclty assigning   =
# = Rare cell populations																				  =
# =========================================================================================================
.fcbMap<-function(x,quantile){
	## TODO This needs to make use of the prior nu0.
	p<-length(x@varNames)
	qq<-qf(quantile,p,x@nu)
	qq<-(x@nu+p)/(x@nu+p*qq)
	(apply(x@z*apply(x@u,2,function(y)ifelse((y<qq),0,1)),1,function(x)ifelse(all(x==0),0,which.max(x))))
}
##Need a wrapper for "det" to handle 1D numeric 1x1 matrix and return the correct result
.det<-function(x){
	if(class(x)=="numeric"){
		x<-as.matrix(x);
	}
	as.matrix(det(x));
}
