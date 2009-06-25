iProcrustes <- function(x, xbar, rotation.only = TRUE, scalling=TRUE, 
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
   if (missing(x))
       stop("iProcrustes: the first argument x is missing.")
   if (missing(xbar))
       stop("iProcrustes: the second arugement xbar is mission.")
   if (!all(dim(x) == dim(xbar)))
       stop("iProcrustes: the dimension of matrix x must equal to that of
            matrix xbar.")    

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
       Q <- .rotateOnly(s)
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


## 
## get non-reflection transformation matrix
##
.rotateOnly <- function(s)
{   ## s: reture values of svd, including U, V and sigma
    ## To lose the reflection component of V*t(V), reverse the sign in the
    ## right most column of either U or V correspoinding to the smallest
    ## sigular value. 

    minsv <- which.min(s$d)
    s$u[, minsv] <- -s$u[, minsv]
    Q <- s$v %*% t(s$u)

    return(Q)
}


