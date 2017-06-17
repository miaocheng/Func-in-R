# 
# This function implements the algorithm of locality preserving projections (LPP).
# Input: X - nDim * nSam
#        para - parameters:
#           r - Reduced dimension
#           k - Parameter for construction of nearest-neighbor graph
# Output: rX - Obtained data with reduced dimension
#         V - Projection vectors of LPP
#         D - Eigenvalues of associated eigenvectors
# Coded by Miao Cheng
#################################################################################

source("./eudist_mc.R")


lpp_mc <- function(X, para){
  
  r <- para$rDim
  k <- para$k
  
  nDim <- dim(X)[1]
  nSam <- dim(X)[2]
  
  D <- eudist_mc(X, X)
  for (ii in 1:nSam)
  {
    D[ii, ii] = 1e10
  }
  
  wets <- NULL
  W <- NULL
  for (i in 1:nSam)
  {
    tmp <- D[i, ]
    ind <- order(tmp)
    ind <- ind[1:k]
    
#     tmp <- array(0, c(k, 3))
#     tmp[, 1] <- i
#     tmp[, 2] <- ind
#     wets <- rbind(wets, tmp)
    
    tmp <- array(0, c(nSam, k))
    tmp[i, ] <- 1
    for (j in 1:k)
    {
      idx <- ind[j]
      tmp[idx, j] <- -1
    }
    W <- cbind(W, tmp)
    
#     wets <- as.matrix(wets)
    W <- as.matrix(W)
    
  }
  
  L <- W %*% t(W)
  
  S <- X %*% L
  S <- S %*% t(X)
  
  SX <- svd(S, nu = nDim, nv = nDim)
  U <- SX$u
  d <- SX$d
#   rW <- length(U[1,])
  V <- subset(U, select = c(nDim:1))
  D <- rev(d)
  
  if(r != 0 && r < nDim)
  {
    V <- V[,1:r]
    D <- D[1:r]
  }
  
  rX <- crossprod(V, X)
  
  output <- list(rX=rX, V=V, D=D)
  return(output)
  
}





