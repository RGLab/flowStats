#' Random Generation from a t Mixture Model with Box-Cox Transformation
#' 
#' This function can be used to generate a sample from a multivariate \eqn{t}
#' mixture model with Box-Cox transformation.
#' 
#' 
#' @param N The number of observations.
#' @param w A vector of length \eqn{K}, containing the \eqn{K} cluster
#' proportions.
#' @param mu A matrix of size \eqn{K \times P}{K x P}, where \eqn{K} is the
#' number of clusters and \eqn{P} is the dimension, containing the \eqn{K} mean
#' vectors.
#' @param sigma An array of dimension \eqn{K \times P \times P}{K x P x P},
#' containing the \eqn{K} covariance matrices.
#' @param nu The degrees of freedom used for the \eqn{t} distribution.
#' @param lambda The Box-Cox transformation parameter.  If missing, the
#' conventional \eqn{t} distribution without transformation will be used.
#' @return A matrix of size \eqn{N \times P}{N x P}.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}
#' @keywords datagen
#' @examples
#' ### Number of components
#' K <- 5
#' ### Dimension
#' p <- 2
#' ### Number of observations
#' n <- 200
#' Mu <- matrix(runif(K*p, 0, 20), K, p)
#' Sigma <- array(0, c(K, p, p))
#' 
#' for (k in 1:K)
#' {
#'     Sigma[k,,][outer(1:p, 1:p, ">")] <- runif(p*(p-1)/2,-.1,.1)
#'     diag(Sigma[k,,]) <- runif(p,0,1)
#'     ### Make sigma positive definite
#'     Sigma[k,,] <- Sigma[k,,] %*% t(Sigma[k,,])
#' }
#' 
#' ### Generate the weights
#' w <- rgamma(K,10,1)
#' w <- w/sum(w)
#' 
#' y <- SimulateMixture(n, w, Mu, Sigma, nu=4)
#' @export 
#' @importFrom mnormt rmt
SimulateMixture <- function(N, w, mu, sigma, nu=4, lambda)
{
    # Number of clusters
    K <- length(w)
    if (K==1) {
        mu <- matrix(mu, 1)
        sigma <- array(sigma, c(1, ncol(mu), ncol(mu)))
    } else if (length(mu)==K) {
        mu <- matrix(mu, K, 1)
        sigma <- array(sigma, c(K, 1, 1))
    }
	# dimension
    py <- ncol(mu)
	y <- matrix(0, N, py)
    nu <- rep(nu,K)
    if (!missing(lambda))
        lambda <- rep(lambda,K)

    label <- sample(1:K, N, replace=T, prob=w)
    count <- table(c(label,1:K))-1
    for (k in 1:K) if (count[k]>0) {    
        y[label==k,] <- rmt(count[k], mu[k,], sigma[k,,], nu[k])
        if (!missing(lambda)) y[label==k,] <- rbox(y[label==k,], lambda[k])    
    }
    y
}
