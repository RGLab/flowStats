useSVD <- function(x, xbar, rotation.Only = FALSE)
{  
   ## Apply SVD to get the rotation matrix (Q) and scaling factor (scal).
   ## Inputs:
   ## The columns of x represent the coordinates of the landmarks, and
   ## the number of row is the labels of the landmarks. x
   ## xbar is the reference landmark. 
   ## Outputs:
   ## Q is the rotation matrix, and scal the scaler (must be positive).

   m <- t(xbar) %*% x
   s <- svd(m)
   Q <- s$v %*% t(s$u) 
 
   ## force Q to be a rotation matrix
   if (det(Q) < 0 & rotation.Only) {
       cat("useSVD: Q is a reflection matrix. Fix it to a rotation matrix...\n")
       Q <- .rotateOnly(Q, x)
   }

   scal <- sum(diag(t(xbar) %*% x %*% Q)) / 
           sum(diag(t(x) %*% x))
   scal <- abs(scal)

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
