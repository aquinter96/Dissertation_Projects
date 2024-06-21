# EMAlgBAdLassoCV_Parallel.R

EMAlgBAdLassoCV <- function(X, xk_sig, Binit, Phi1init, tuningpB = seq(5.1,5.5,0.1), nfolds=10, weights){
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  B0 <- matrix(0, nrow = q, ncol = 1, byrow = T)
  B0new <-  matrix(colMeans(X), nrow = q, ncol = 1, byrow = T)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Bnew <- Binit
  Phi1 <- matrix(0, nrow = q, ncol = q)
  Phi1new <- Phi1init
  
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (niter < 15000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    ## Here grid search is used to find optimal tuning parameters for B in the LASSO part. 
    #for the first iteration of the algorithm ONLY, use grid search to calculate the optimal tuning parameter combination
    #for B
    if(niter == 0){
      
      cv.start <- Sys.time()
      
      ## Grid search for the optimal tuning parameter is performed using the parSapply() function for parallelization.
      ## Each tuning parameter being checked is applied to the tunesearchB function which calculates the corresponding
      ## BIC and stored in the BICs vector
      BICs <- parSapply(cl, tuningpB, tunesearchB, X, B0, B, Phi1, weights)
      
      if(length(BICs) > 1){
        sigmaNAs <- is.na(BICs)
        tuningpB <- tuningpB[!sigmaNAs]
        BICs <- BICs[!sigmaNAs]
      }
      if(length(BICs) == 0){
        if(length(tuningpBinit) == 1){
          tuningpBhat <- tuningpBinit
        }
        if(length(tuningpBinit) > 1){
          tuningpBhat <- min(tuningpBinit)
        }
      }else{
        tuningpBhat <- tuningpB[which(BICs == min(BICs))]
      }
      
      #once BIC has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest BIC (i.e. the optimal tuning parameter)
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    invmodv <- Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
    
    EW <- crossprod(B,invmodv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
    
    for(i in 1:xk_sig){
      for(j in 1:xk_sig){
        deltaijB <- 0
        thetaijB <- 0
        omegaijB <- 0
        if(j < i){
          deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
          thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i,], n)))
          omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%Bnew[i,-j]
          bbar <- ((thetaijB-omegaijB))/(deltaijB)
          if(weights[i,j] == 0){
            lambdaB <- (tuningpBhat*Phi1[i,i])/deltaijB*(1/1e-5)
          }
          else{
            lambdaB <- (tuningpBhat*Phi1[i,i])/deltaijB*(1/abs(weights[i,j]))
          }
          Bnew[i,j] <- lasso(bbar, lambdaB)
        }
      }
    }
    for(i in (xk_sig+1):q){
      for(j in 1:xk_sig){
        deltaijB <- 0
        thetaijB <- 0
        omegaijB <- 0
        deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
        thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i,], n)))
        omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%Bnew[i,-j]
        bbar <- ((thetaijB-omegaijB))/(deltaijB)
        if(weights[i,j] == 0){
          lambdaB <- (tuningpBhat*Phi1[i,i])/deltaijB*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpBhat*Phi1[i,i])/deltaijB*(1/abs(weights[i,j]))
        }
        Bnew[i,j] <- lasso(bbar, lambdaB)
      }
    }
    
    B0new <- (1/n)*as.matrix(colSums(X - t(Bnew%*%EW)))
    Phi1new <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0new, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(Bnew) + matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0new, nrow = q, ncol = n, byrow = F)) +
                                 2*matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(Bnew) + n*Bnew%*%condvar%*%t(Bnew) + Bnew%*%EW%*%t(EW)%*%t(Bnew)))
    
    niter <- niter + 1
  }
  B0 <- B0new
  B <- Bnew
  Phi1 <- Phi1new
  
  invmodv <- Matrix::solve(Phi1 + tcrossprod(B))
  condvar <- diag(xk_sig) - crossprod(B,invmodv)%*%B
  EW <- crossprod(B,invmodv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
  
  log.lik <- -(n/2)*(sum(diag(invmodv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1 +tcrossprod(B))))
  
  BICopt <- log(n)*(sum(B != 0) - xk_sig) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter, "optimal lambda" = tuningpBhat, "X" = X, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}

#Function for calculating the initial values of the Y model for JSFRM
initvalcalc <- function(X, Y, Binit, s, m, n){
  
  AGamma <- t(solve(t(Binit)%*%Binit)%*%t(Binit)%*%cov(X, Y))
  
  p <- ncol(Y)
  
  Aobj <-function(a){
    amat <- matrix(0, nrow = p, ncol = m)
    amat[(1:m),(1:m)][lower.tri(amat[(1:m),(1:m)], diag = T)] <- a[1:(m*(m+1)/2)]
    amat[(m+1):p,] <- a[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
    phi2mat <- diag(a[(p*m-m*(m-1)/2)+(1:p)])
    norm((var(Y) - AGamma%*%t(AGamma) - amat%*%t(amat) - phi2mat), type = "F")
  }
  
  Ainit <- matrix(0, nrow = p, ncol = m)
  if(m==1){
    Ainit[1,1] <- 1
  }
  else{
    diag(Ainit[(1:m),(1:m)]) <- 1
  }
  Phi2init <- diag(p)
  Ainitvec <- c(Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)], Ainit[(m+1):p,])
  Phi2initvec <- diag(Phi2init)
  Anewvec <- nlminb(c(Ainitvec, Phi2initvec), Aobj)
  Ainit <- matrix(0, nrow = p, ncol = m)
  Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)] <- Anewvec$par[1:(m*(m+1)/2)]
  Ainit[(m+1):p,] <- Anewvec$par[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
  
  Phi2init <- diag(Anewvec$par[(p*m-m*(m-1)/2)+(1:p)])
  A3 <- Ainit%*%(eigen(1/n*t(Ainit)%*%Ainit)$vectors)
  A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))
  
  if(m==1){
    Phi3init <- Ainit[1,1]^2
    Ainit <- A5 %*% solve(A5[(1:m),(1:m)])
    Ainit[1,1] <- 1
  }
  else{
    Phi3init <- diag(diag(A5[(1:m),(1:m)]))%*%diag(diag(A5[(1:m),(1:m)]))
    Ainit <- A5 %*% solve(diag(diag(A5[(1:m),(1:m)])))
    Ainit[(1:m),(1:(m))][upper.tri(Ainit[(1:m),(1:(m))], diag = F)] <- 0
    diag(Ainit[(1:m),(1:(m))]) <- 1
  }
  AAinv <- Matrix::solve(t(Ainit)%*%Ainit)
  Gammainit <- AAinv%*%t(Ainit)%*%AGamma
  output <- list("Ainit" = Ainit, "Ginit" = Gammainit, "Phi2init" = Phi2init, "Phi3init" = Phi3init)
  return(output)
}

## end of code