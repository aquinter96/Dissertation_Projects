## C_inits.R

C_inits <- function(Data, Mest, t){
  
  ############################################################################ 
  ## Obtain estimate of the matrix product A*Gamma using method of moments
  ############################################################################ 
  
  if(is.null(dim(Mest$X_estimates$X_final_pars$B))){
    B <- as.matrix(Mest$X_estimates$X_final_pars$B)
  }else{
    B <- Mest$X_estimates$X_final_pars$B
  }
  if(is.null(dim(Mest$M_estimates$M_final_pars$A))){
    A <- as.matrix(Mest$M_estimates$M_final_pars$A)
  }else{
    A <- Mest$M_estimates$M_final_pars$A
  }
  Gamma <- Mest$M_estimates$M_final_pars$Gamma
  
  Psi <- Mest$X_estimates$X_final_pars$Psi
  Phi2 <- Mest$M_estimates$M_final_pars$Phi2
  Phi3 <- Mest$M_estimates$M_final_pars$Phi3
  
  X <- Data$X
  M <- Data$M
  Y <- Data$Y
  
  n <- nrow(X)
  r <- ncol(Y)
  m <- ncol(M)
  
  inits <- list()
  
  C0init <- as.matrix(colMeans(Y))
  
  XCov1 <- solve(Psi%*%t(B)%*%B%*%Psi)%*%Psi%*%t(B)%*%cov(X,Y)
  MCov1 <- solve(t(A)%*%A)%*%t(A)%*%cov(M, Y)
  
  CDelta <- t(solve(Phi2)%*%(MCov1 - Gamma%*%Psi%*%XCov1))
  Dinit <- t(XCov1 - t(Gamma)%*%t(CDelta))
  
  ############################################################################ 
  ## Use sample covariance matrix and estimate of A*Gamma to get the initial
  ## estimates of A/Phi2 using the nlminb() function
  ############################################################################ 
  
  Cobj <-function(c){
    cmat <- matrix(0, nrow = r, ncol = t)
    cmat[(1:t),(1:t)][lower.tri(cmat[(1:t),(1:t)], diag = T)] <- c[1:(t*(t+1)/2)]
    cmat[(t+1):r,] <- c[(1+t*(t+1)/2):(r*t-t*(t-1)/2)]
    phi5mat <- diag(c[(r*t-t*(t-1)/2)+(1:r)])
    norm((var(Y) - CDelta%*%(Gamma%*%Psi%*%t(Gamma) + Phi2)%*%t(CDelta) - cmat%*%t(cmat) - Dinit%*%Psi%*%t(Dinit) - 2*CDelta%*%Gamma%*%Psi%*%t(Dinit) - phi5mat), type = "F")
  }
  
  Cinit <- matrix(0, nrow = r, ncol = t)
  if(t==1){
    Cinit[1,1] <- 1
  }
  else{
    diag(Cinit[(1:t),(1:t)]) <- 1
  }
  
  
  Phi5init <- diag(r)
  Cinitvec <- c(Cinit[(1:t),(1:t)][lower.tri(Cinit[(1:t),(1:t)], diag = T)], Cinit[(t+1):r,])
  Phi5initvec <- diag(Phi5init)
  Cnewvec <- nlminb(c(Cinitvec, Phi5initvec), Cobj)
  Cinit <- matrix(0, nrow = r, ncol = t)
  Cinit[(1:t),(1:t)][lower.tri(Cinit[(1:t),(1:t)], diag = T)] <- Cnewvec$par[1:(t*(t+1)/2)]
  Cinit[(t+1):r,] <- Cnewvec$par[(1+t*(t+1)/2):(r*t-t*(t-1)/2)]
  
  Phi5init <- diag(Cnewvec$par[(r*t-t*(t-1)/2)+(1:r)])
  C3 <- Cinit%*%(eigen(1/n*t(Cinit)%*%Cinit)$vectors)
  C5 <- C3%*%qr.Q(qr(t(C3[(1:t),(1:t)])))
  
  ############################################################################ 
  ## Use estimates of A/A*Gamma/Phi2 to get initial estimates of Gamma and Phi3
  ############################################################################ 
  
  if(t==1){
    Phi4init <- C5[1,1]^2
    Cinit <- C5 %*% solve(C5[(1:t),(1:t)])
    Cinit[1,1] <- 1
  }
  else{
    Phi4init <- diag(diag(C5[(1:t),(1:t)]))%*%diag(diag(C5[(1:t),(1:t)]))
    Cinit <- C5 %*% solve(diag(diag(C5[(1:t),(1:t)])))
    Cinit[(1:t),(1:t)][upper.tri(Cinit[(1:t),(1:t)], diag = F)] <- 0
    diag(Cinit[(1:t),(1:t)]) <- 1
  }
  CCinv <- Matrix::solve(t(Cinit)%*%Cinit)
  # Deltainit <- CCinv%*%t(Cinit)%*%CDelta
  Deltainit <- Matrix::solve(t(C5)%*%C5)%*%t(C5)%*%CDelta
  
  inits$C0 <- C0init
  inits$C <- Cinit
  inits$D <- Dinit
  inits$Delta <- Deltainit
  inits$Phi4 <- Phi4init
  inits$Phi5 <- Phi5init
  
  return(inits)
}

## end of code