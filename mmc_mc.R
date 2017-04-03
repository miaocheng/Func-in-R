#
# This function implements the efficient version of maximum margin criterion (MMC),
# in which high-dimensional data can be handled.
# Input:   X - nDim * nSam
#          L - Labels of training set
#          r - Reduced dimension
# Output:  rX - Obtained data with reduced dimension
#          ldaV - Projection vectors of MMC
#          wD - Eigenvalues of associated eigenvectors
#          
# Coded by Miao Cheng



# library(base)


mmc_mc <- function(X, L, r = 0)
{
  nSam = length(X[1,])
  nDim = length(X[,1])
  X <- as.matrix(X)
  L <- as.matrix(L)
  meanx <- rowMeans(X)
  
  xL <- unique(L)
  nL <- length(xL)
  B <- NULL
  W <- NULL
  for (i in 1:nL)      #seq(from=1, to=nL, by=1))
  {
    # ml <- L %in% xL(i)
    # ml <- match(L, xL(i), nomatch=0)
    ml <- which(L==xL[i], arr.ind=TRUE)
    nC <- length(ml)
    cX <- subset(X, select = ml)
    
    # Constructing between and within scatter
    meanc <- rowMeans(cX)
    B <- cbind(B, meanc)
    wX <- cX - rep(meanc, nC)
    W <- cbind(W, wX)
  }
  B <- B - rep(meanx, nL)
  
  if (nDim > 1500)
  {
    B <- t(X) %*% B
    W <- t(X) %*% W
  }
  
  Sb <- B %*% t(B)
  Sw <- W %*% t(W)
  Sb <- Sb / nL
  Sw <- Sw / nSam
  St <- Sb - Sw
  
  resT = svd(St)
  U <- resT$u
  S <- resT$d
  
  if (nDim > 1500)
  {
    ss <- S^(-0.5)
    ss <- diag(ss)
    tmp <- X %*% U
    bV <- tmp %*% ss
    
    resV <- svd(bV)
    U <- resV$u
    S <- resV$d
    
    S <- S^(-2)
    
    S <- rev(S)
    rT <- length(U[1,])
    U <- subset(U, select = c(rT:1))
  }
  
  if(r != 0 && r < rT)
  {
    U <- U[,1:r]
    S <- S[1:r]
  }
  
  
  rX <- crossprod(U, X)
  
  output <- list(rX=rX, V=U, D=S)
  
  return(output)
  
  
}




