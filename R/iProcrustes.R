iProcrustes <- function(x, xbar, rotation.only = FALSE, scalling=TRUE, 
                   translate=FALSE)
{  
   ## Apply SVD to get the rotation matrix (Q) and scaling factor (scal).
   ## Inputs:
   ## The columns of x represent the coordinates of the landmarks, and
   ## the number of row is the labels of the landmarks. x
   ## xbar is the reference landmark. 
   ## Outputs:
   ## Q is the rotation matrix, and scal the scaler (must be positive).

   ## check valid dimension
   if (!all(dim(x) == dim(xbar)))
      stop("The dimension of matrix x must equal to that of matrix xbar.")

   ## 1. translation
   if (translate) {
       ## translate the centroids of the points in x to that of points in xbar
       T <-  colMeans(x, na.rm=TRUE) 
       I <- matrix(1, nrow=nrow(x), ncol=1)
       x <- x - (I %*% T)
       xbar <- xbar - (I %*% colMeans(xbar, na.rm=TRUE))
   }

   ## 2. get rotation/reflection matrix
   m <- t(xbar) %*% x
   s <- svd(m)
   Q <- s$v %*% t(s$u) 
 
   ## force Q to be a rotation matrix
   if (det(Q) < 0 & rotation.only) {
       cat("iProcrustes: Q is a reflection matrix. Fix it to a rotation matrix...\n")
       Q <- .rotateOnly(Q, x)
   }

   ## 3. get scalling factor
   if (scalling) {
       scal <- sum(diag(t(xbar) %*% x %*% Q)) / 
               sum(diag(t(x) %*% x))
       scal <- abs(scal)
   }
   else 
       scal <- 1

   # returning results
   if (translate)
       return(list(T=T, Q=Q, scal=scal, T.xbar=xbar))
   else
       return(list(Q=Q, scal=scal))
}


## anti-reflection, rotation only matrix
## now only works for the reflection across Y-axis

.rotateOnly <- function(Q, x)
{
    ## this function only 2-D
    if (det(Q) == 1)
    {
        cat("Q is already a rotational matrix. Return its original value.")
        r2psi <- Q
    } 
    else 
    {   ## reflect across Y-axis only, for now... #degugging
        refX<- matrix(c(1, 0, 0 ,-1), nrow=2, ncol=2)
    	
	A <-  x %*% Q
    	xhat <- x %*% Q %*% refX
    	norm <- sqrt(xhat[1,1]^2+xhat[1,2]^2)
    	psi = acos(sum(xhat[2,]*c(1,0))/norm)

        ## rotate counterclockwise 
    	r2psi = matrix(c(cos(2*psi), -sin(2*psi), sin(2*psi), cos(2*psi)),
      	       ncol=2, nrow=2)

   	xtil = xhat %*% r2psi ## for debugging
   }
   return(r2psi)
}

## it will be better to write  a function to seperate rotation and reflection
## matrices. 
