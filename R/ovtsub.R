##' This function computes an Overton-like subtraction of two densities. It calculates the proportion of the reference density that is above a reference
##'
##' The test can be one-sided or two-sided. If one sided, it tests the region of the test density that is above the mode of the reference density. If two-sided it will look at the regions on either side of the mode of the reference density.
##' Densities are computed on a grid of 1024, and appropriately normalized. 
##' @title Overton-like subtraction of densities.
##' @param ref The reference channel specified as a \code{vector}
##' @param test The test (potentially positive) channel specified as a \code{vector}
##' @param twosided \code{boolean} flag testing whether the area of the density of the test curve above the reference curve will be calculated on both sides of the mode of the test curve (TRUE) or only on the positive side of the mode (FALSE, default). 
##' @return \code{numeric} value representing the proportion of the area of the test density above the reference density. 
##' @author Greg Finak
##' @examples
##' A = rnorm(10000,mean=1,sd=0.5)
##' B = rnorm(10000,mean=2,sd=0.5)
##' overton_like(A,B)
##' 
overton_like = 
function(ref, test,twosided=FALSE) {
	from = pmin(range(ref), range(test))[1]
	to = pmax(range(ref), range(test))[2]
	ref = density(ref, from = from, to = to,n=1024)
	test = density(test, from = from, to = to,n=1024)
	muA = sum(ref$y * ref$x * diff(ref$x)[1])
	ABnorm = sum(pmax(ref$y, test$y) * diff(ref$x)[1])
	bpart = test$y * diff(test$x)[1] / ABnorm
	apart = ref$y * diff(ref$x)[1] / ABnorm
	aboverlap = pmin(ref$y, test$y) * diff(ref$x)[1] / ABnorm
	bpospart = bpart
	aboverlappos = aboverlap
	if(!twosided){
		bpospart[test$x < muA & test$y / ref$y < 1] = 0 #zero out where test < ref and less than mu
		bpospart[test$x > muA & test$y / ref$y > 1] = 0 #zero out where test > ref and more than mu
		aboverlappos[test$x < muA & test$y / ref$y < 1] = 0
		aboverlappos[test$x > muA & test$y / ref$y > 1] = 0
	}else{
		bpart[test$x<muA]=-bpart[test$x<muA]
		aboverlap[test$x<muA]  = -aboverlap[test$x<muA]
	}
	if(!twosided){
		res = (sum(bpart) - sum(aboverlap) - sum(bpospart) + sum(aboverlappos)) *
	ABnorm
	}else{
		res=(sum(bpart)-sum(aboverlap))*ABnorm
	}
	if (res < 0 & res > -1) {
		res = res
	} else if (res > 1) {
		res = 1
	} else if (res < -1){
		res = -1
	}
	res
}
