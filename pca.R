
# This function does principal component analysis for input data.class.
# Input:   X - nDim * nSam
#          r - Reduced dimension
# Output:  rX - Obtained data with reduced dimension
#          eigvector - Principal components
#          eigvalue - Eigenvalues of associated eigenvectors
#          
# Coded by Miao Cheng
# Date: 2015-11-5

# library(base)


pca_mc <- function(X, r = 0)
{
  nSam = length(X[1,])
  nDim = length(X[,1])
  X <- as.matrix(X)
  
  meanx <- rowMeans(X)
  cX <- X - rep(meanx, nSam)
  
  if(nDim < nSam)
  {
    XX <- cX*t(cX)
    res <- eigen(XX)
    eigenvector <- res$vectors
    eigenvalue <- res$values
  }
  else
  {
    res <- svd(cX)
    D <- diag(res$d)
    D <- D*D
    eigenvalue <- D
    eigenvector <- res$u
  }
  
  if(r != 0 && r > 0)
  {
    eigenvector <- eigenvector[,1:r]
    eigenvalue <- eigenvalue[,1:r]
  }
  
  output <- list(pcaV=eigenvector, pcaD=eigenvalue)
  
  return(output)
  
}




