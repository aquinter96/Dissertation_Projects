## B_inits.R

B_inits <- function(X, s){
  
  #define initial values for B0, B, and Phi1
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  Xcent <- scale(X, scale = FALSE)
  B0 <- matrix(0, nrow = q, ncol = 1, byrow = T)
  B0new <-  matrix(colMeans(X), nrow = q, ncol = 1, byrow = T)
  B <- matrix(0, nrow = q, ncol = s)
  Bnew <- as.matrix((svd(Xcent)$v%*%diag(svd(Xcent)$d)/sqrt(n))[,1:s])
  Phi1 <- matrix(0, nrow = q, ncol = q)
  Phi1new <- diag(q)
  
  ##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
  B31 <- Bnew%*%(eigen(1/(ncol(X))*crossprod(Bnew))$vectors)
  B51 <- B31%*%qr.Q(qr(t(B31[(1:s),(1:s)])))
  if(s == 1){
    Bnew <- B51 %*% solve(B51[(1:s),(1:s)])
    Bnew[1,1] <- 1
  }else{
    Bnew <- B51 %*% solve(diag(diag(B51[(1:s),(1:s)])))
    Bnew[(1:s),(1:(s))][upper.tri(Bnew[(1:s),(1:(s))], diag = F)] <- 0
    diag(Bnew[(1:s),(1:(s))]) <- 1
  }
  
  inits <- list()
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(((norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001)) & (niter < 5000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    invmodv <-Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(s) - t(B)%*%invmodv%*%B
    
    EW <- t(B)%*%invmodv%*%(t(X)-matrix(B0, nrow = ncol(X), ncol = n, byrow = F))
    B0new <- (1/n)*as.matrix(colSums(X - t(B%*%EW)))
    Bnew <- (crossprod(X,t(EW)) - matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(EW))%*%solve(n*condvar + tcrossprod(EW))
    Phi1new <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0new, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(Bnew) + matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0new, nrow = q, ncol = n, byrow = F)) + 
                                 2*matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(Bnew) + n*Bnew%*%condvar%*%t(Bnew) + Bnew%*%EW%*%t(EW)%*%t(Bnew)))
    niter <- niter + 1
  }
  B0 <- B0new
  Phi1 <- Phi1new
  Bem <- Bnew
  B3 <- Bem%*%(eigen(1/(ncol(X))*crossprod(Bem))$vectors)
  B5 <- B3%*%qr.Q(qr(t(B3[(1:s),(1:s)])))
  if(s==1){
    B <- B5 %*% solve(B5[(1:s),(1:s)])
    B[1,1] <- 1
  }
  else{
    B <- B5 %*% solve(diag(diag(B5[(1:s),(1:s)])))
    B[(1:s),(1:(s))][upper.tri(B[(1:s),(1:(s))], diag = F)] <- 0
    diag(B[(1:s),(1:(s))]) <- 1
  }
  
  inits$B0 <- B0
  inits$B <- B
  inits$Phi1 <- Phi1
  
  return(inits)
}

## end of code