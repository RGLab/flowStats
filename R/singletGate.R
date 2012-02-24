singletGate<-function(x, area,height,sidescatter,lower = NULL, 
    upper = NULL,maxit=100,nsd=5) {
    flowCore:::checkClass(x, "flowFrame")
    flowCore:::checkClass(area, "character")
    flowCore:::checkClass(height, "character")
    if (length(area)+length(height) != 2) {
        stop("stains must be of length 2, specifying the FSC area and FSC height channels")
    }
    if(!(is.null(lower)&is.null(upper))){
	flowCore:::checkClass(lower, "numeric")
    	flowCore:::checkClass(upper, "numeric")

    	if (lower > 1 | lower < 0 | upper > 1 | upper < 0) {
        	stop("lower and upper must be in the range [0,1]")
    	}
	lower<-qnorm(lower)
	upper<-qnorm(upper)
    }
    
	if(!is.null(nsd)){
	    flowCore:::checkClass(nsd,"numeric")
	    lower<- -abs(nsd);
	    upper<-abs(nsd);
    }
    x <- exprs(x[, c(area,height,sidescatter)])
    newcols <- c("A", "H","SSC")
    oldcols <- colnames(x)
    colnames(x) <- newcols
    #Model the forward scatter height as a function of area + side scatter + side scatter / area (of the models tested, this gave the best R-squared, around 0.78)
    form <- as.formula("H~A+SSC+I(SSC/A)")
    model <- rlm(form, as.data.frame(x),maxit=100)
	#TODO throw a warning or something if the model doesn't converge.
    est <- huber(resid(model))
    #Filter outliers based on the residuals. Threshold taken from lower/upper, either number of sd, or a quantile. 
    indices <- findInterval(resid(model) ,c(est$mu + lower * est$s,  est$mu + 
        upper * est$s))==1
    serr<-sum(resid(model)^2);sreg<-sum((fitted(model)-huber(x[,"H"])$mu)^2)
    stot<-sreg+serr
    R<-serr/stot
    retme<-list(indices=indices,Rsquared=R)
    class(retme)<-"singletFilter"
 	return(retme)
}
