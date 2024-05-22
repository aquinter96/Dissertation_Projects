
library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/coord_desc.cpp")

datgen <- function(ak, bk, gamma, n,varvecB,varvecA,varvecG){
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    W <- matrix(rnorm(n, mean = 0, sd = 1), nrow = 1, n) ## generate matrix U N(0,I) m x n
  }
  else{
    W <- matrix(rnorm(n*ncol(gamma), mean = 0, sd = 1), nrow = ncol(gamma), n) ## generate matrix U N(0,I) m x n
  }
  xlin.term  <-  bk %*% rbind(rep(1, ncol(W)),W)
  xlin.term <- t(xlin.term)
  xx <- matrix(0, nrow = n, ncol = nrow(bk))
  for(jj in 1: (ncol(xx))){
    xx[, jj] <- xlin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecB[jj])))
  }
  colnames(xlin.term) <- paste0('X', 1:(ncol(xx)))
  colnames(xx) <- paste0('X', 1:(ncol(xx)))
  X = xx
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    Z <- gamma*W + matrix(rnorm(n, mean = 0, sqrt(varvecG)), nrow = 1) ## generate latent matrix Z
  }
  else{
    Z <- matrix(0, nrow=nrow(gamma), ncol = n)
    for(i in 1:nrow(gamma)){
      Z[i,] <- gamma[i,]%*%W + matrix(rnorm(n, mean = 0, sd = sqrt(varvecG[i])), nrow = 1) ## generate latent matrix Z
    }
  }
  ylin.term  <-  ak %*% rbind(rep(1, ncol(Z)),Z)
  ylin.term <- t(ylin.term)
  yy <- matrix(0, nrow = n, ncol = nrow(ak))
  for(jj in 1: (ncol(yy))){
    yy[, jj] <- ylin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecA[jj])))
  }
  colnames(ylin.term) <- paste0('Y', 1:(ncol(yy)))
  colnames(yy) <- paste0('Y', 1:(ncol(yy)))
  Y = yy
  output <- list("X" = X, "Y" = Y, "Z" = Z)
}


#Function for basic EM-based factor model for initial estimates of B0, B, and Phi1
EMAlgB <- function(X, xk_sig){
  
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  Xcent <- scale(X, scale = FALSE)
  B0 <- matrix(0, nrow = q, ncol = 1, byrow = T)
  B0new <-  matrix(colMeans(X), nrow = q, ncol = 1, byrow = T)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Bnew <- as.matrix((svd(Xcent)$v%*%diag(svd(Xcent)$d)/sqrt(n))[,1:xk_sig])
  Phi1 <- matrix(0, nrow = q, ncol = q)
  Phi1new <- diag(q)
  
  B31 <- Bnew%*%(eigen(1/(ncol(X))*crossprod(Bnew))$vectors)
  B51 <- B31%*%qr.Q(qr(t(B31[(1:xk_sig),(1:xk_sig)])))
  if(xk_sig == 1){
    Bnew <- B51 %*% solve(B51[(1:xk_sig),(1:xk_sig)])
    Bnew[1,1] <- 1
  }else{
    Bnew <- B51 %*% solve(diag(diag(B51[(1:xk_sig),(1:xk_sig)])))
    Bnew[(1:xk_sig),(1:(xk_sig))][upper.tri(Bnew[(1:xk_sig),(1:(xk_sig))], diag = F)] <- 0
    diag(Bnew[(1:xk_sig),(1:(xk_sig))]) <- 1
  }
  
  while(((norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001)) & (niter < 5000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    invmodv <-Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
    
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
  B5 <- B3%*%qr.Q(qr(t(B3[(1:xk_sig),(1:xk_sig)])))
  if(xk_sig==1){
    B <- B5 %*% solve(B5[(1:xk_sig),(1:xk_sig)])
    B[1,1] <- 1
  }
  else{
    B <- B5 %*% solve(diag(diag(B5[(1:xk_sig),(1:xk_sig)])))
    B[(1:xk_sig),(1:(xk_sig))][upper.tri(B[(1:xk_sig),(1:(xk_sig))], diag = F)] <- 0
    diag(B[(1:xk_sig),(1:(xk_sig))]) <- 1
  }
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter)
  return(output)
}





#  EM Algorithm with only R code
EMAlgBAdLassoCV <- function(X, xk_sig, Binit, Phi1init, tuningpB = seq(5.1,5.5,0.1), nfolds=10, weights){
  alg.start <- Sys.time()
  
  tuningpBinit <- tuningpB
  
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
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  varX <- cov(X)
  
  while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (niter < 15000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
    if(niter == 0){
      
      deltaijB <- 0 #initialize coordinate descent starting values
      thetaijB <- 0
      omegaijB <- 0
      lambdaB <- NULL
      bbar <- NULL
      
      B0CV <- rep(0, q)
      BCV <- matrix(0, nrow = n, ncol = q)
      Phi1CV <- matrix(0, nrow = q, ncol = q)
      B0start <- B0
      Bstart <- B
      Phi1start <- Phi1
      
      BICmat <- rep(0,length(tuningpB))
      
      cv.start <- Sys.time()
      maxit <- 0
      
      for(a in 1:length(tuningpB)){
        
        B0CV <- B0start
        BCV <- Bstart
        Phi1CV <- Phi1start
        
        B0 <- rep(0, q)
        B <- matrix(0, nrow = q, ncol = xk_sig)
        Phi1 <- matrix(0, nrow = q, ncol = q)
        
        maxit <- 0
        
        while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (maxit < 15000)){
          
          B0 <- B0CV
          B <- BCV
          Phi1 <- Phi1CV
          
          invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
          condvar <- diag(xk_sig) - crossprod(B,invmodcv)%*%B
          EW <- crossprod(B,invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
          
          BCVupdate <- matrix(0, nrow = q, ncol = xk_sig)
          
          while(norm(BCV - BCVupdate, type = "F") > 0.001){
            
            BCVupdate <- BCV
            
            for(i in 1:xk_sig){
              for(j in 1:xk_sig){
                deltaijB <- 0
                thetaijB <- 0
                omegaijB <- 0
                #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
                if(j < i){
                  deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
                  thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0CV[i,], n)))
                  omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%BCV[i,-j]
                  bbar <- ((thetaijB-omegaijB))/(deltaijB)
                  if(weights[i,j] == 0){
                    lambdaB <- (tuningpB[a]*Phi1CV[i,i])/deltaijB*(1/1e-5)
                  }
                  else{
                    lambdaB <- (tuningpB[a]*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
                  }
                  BCV[i,j] <- lasso(bbar, lambdaB)
                }
              }
            }
            #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
            for(i in (xk_sig+1):q){
              for(j in 1:xk_sig){
                deltaijB <- 0
                thetaijB <- 0
                omegaijB <- 0
                deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
                thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0CV[i,], n)))
                omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%BCV[i,-j]
                bbar <- ((thetaijB-omegaijB))/(deltaijB)
                if(weights[i,j] == 0){
                  lambdaB <- (tuningpB[a]*Phi1CV[i,i])/deltaijB*(1/1e-5)
                }
                else{
                  lambdaB <- (tuningpB[a]*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
                }
                BCV[i,j] <- lasso(bbar, lambdaB)
              }
            }
          }
          #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
          #to check for convergence
          B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
          Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                                      2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
          maxit <- maxit + 1
        }
        #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
        invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
        log.lik <- -(n/2)*(sum(diag(invmodcv%*%(varX+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                             q*log((2*pi))+log(det(Phi1CV +tcrossprod(BCV))))
        
        BICmat[a] <- log(n)*(sum(BCV != 0) - xk_sig) - 2*log.lik
        
        #calculate CV MSE for each tuning parameter
        #using a running sum over each fold 
      }
      if(length(BICmat) > 1){
        sigmaNAs <- is.na(BICmat)
        tuningpB <- tuningpB[!sigmaNAs]
        BICmat <- BICmat[!sigmaNAs]
      }
      if(length(BICmat) == 0){
        if(length(tuningpBinit) == 1){
          tuningpBhat <- tuningpBinit
        }
        if(length(tuningpBinit) > 1){
          tuningpBhat <- min(tuningpBinit)
        }
      }else{
        tuningpBhat <- tuningpB[which(BICmat == min(BICmat))]
      }
      
      B0 <- B0start
      B <- Bstart
      Phi1 <- Phi1start
      
      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    invmodv <- Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
    
    EW <- crossprod(B,invmodv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
    
    Bnewupdate <- matrix(0, nrow = q, ncol = xk_sig)
    
    while(norm(Bnew - Bnewupdate) > 0.00001){
      
      Bnewupdate <- Bnew
      
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
  
  log.lik <- -(n/2)*(sum(diag(invmodv%*%(varX+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1 +tcrossprod(B))))
  
  BICopt <- log(n)*(sum(B != 0) - xk_sig) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter, "optimal lambda" = tuningpBhat, "X" = X, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}




#  EM ALgorithm with C++ code for coordinate descent
EMAlgBAdLassoCV2 <- function(X, xk_sig, Binit, Phi1init, tuningpB = seq(5.1,5.5,0.1), nfolds=10, weights){
  alg.start <- Sys.time()
  
  tuningpBinit <- tuningpB
  
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
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  varX <- cov(X)
  
  while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (niter < 15000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
    if(niter == 0){
      
      deltaijB <- 0 #initialize coordinate descent starting values
      thetaijB <- 0
      omegaijB <- 0
      lambdaB <- NULL
      bbar <- NULL
      
      B0CV <- rep(0, q)
      BCV <- matrix(0, nrow = n, ncol = q)
      Phi1CV <- matrix(0, nrow = q, ncol = q)
      B0start <- B0
      Bstart <- B
      Phi1start <- Phi1
      
      BICmat <- rep(0,length(tuningpB))
      
      cv.start <- Sys.time()
      maxit <- 0
      
      for(a in 1:length(tuningpB)){
        
        B0CV <- B0start
        BCV <- Bstart
        Phi1CV <- Phi1start
        
        B0 <- rep(0, q)
        B <- matrix(0, nrow = q, ncol = xk_sig)
        Phi1 <- matrix(0, nrow = q, ncol = q)
        
        maxit <- 0
        
        while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (maxit < 15000)){
          
          B0 <- B0CV
          B <- BCV
          Phi1 <- Phi1CV
          
          invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
          condvar <- diag(xk_sig) - crossprod(B,invmodcv)%*%B
          EW <- crossprod(B,invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
          
          BCV <- coordesc(n, B, B0, X, EW, condvar, Phi1, tuningpB[a], weights)
          
                    #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
          #to check for convergence
          B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
          Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                                      2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
          maxit <- maxit + 1
        }
        #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
        invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
        log.lik <- -(n/2)*(sum(diag(invmodcv%*%(varX+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                             q*log((2*pi))+log(det(Phi1CV +tcrossprod(BCV))))
        
        BICmat[a] <- log(n)*(sum(BCV != 0) - xk_sig) - 2*log.lik
        
        #calculate CV MSE for each tuning parameter
        #using a running sum over each fold 
      }
      if(length(BICmat) > 1){
        sigmaNAs <- is.na(BICmat)
        tuningpB <- tuningpB[!sigmaNAs]
        BICmat <- BICmat[!sigmaNAs]
      }
      if(length(BICmat) == 0){
        if(length(tuningpBinit) == 1){
          tuningpBhat <- tuningpBinit
        }
        if(length(tuningpBinit) > 1){
          tuningpBhat <- min(tuningpBinit)
        }
      }else{
        tuningpBhat <- tuningpB[which(BICmat == min(BICmat))]
      }
      
      B0 <- B0start
      B <- Bstart
      Phi1 <- Phi1start
      
      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    invmodv <- Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
    
    EW <- crossprod(B,invmodv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
    
    Bnewupdate <- matrix(0, nrow = q, ncol = xk_sig)
    
    Bnew <- coordesc(n, B, B0, X, EW, condvar, Phi1, tuningpBhat, weights)
    
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
  
  log.lik <- -(n/2)*(sum(diag(invmodv%*%(varX+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1 +tcrossprod(B))))
  
  BICopt <- log(n)*(sum(B != 0) - xk_sig) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter, "optimal lambda" = tuningpBhat, "X" = X, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}

#Function that uses CV to select the optimal X model across a range of possible values for the number of latent factors s
OverallBAlg <- function(X, tuningpB = seq(0,15,0.1), s_seq = 1:4, nfolds=10){
  Blist <- replicate(length(s_seq),NA,simplify=F)
  BIClist <- rep(0, length(s_seq))
  for(i in 1:length(s_seq)){
    Binit <- EMAlgB(X, s_seq[i])
    Blist[[i]] <- EMAlgBAdLassoCV(X, s_seq[i], Binit$B, Binit$Phi1, tuningpB, weights = Binit$B)
    BIClist[i] <- Blist[[i]]$BICopt
  }
  sopt <- s_seq[which(BIClist == min(BIClist))]
  optB <- Blist[[which(BIClist == min(BIClist))]]
  return(list("Bopt" = optB, "optimal s" = sopt, "BIC" = optB$BICopt, "lambda" = optB$`optimal lambda`))
}

#Function that uses CV to select the optimal Y model across a range of possible values for the number of latent factors m given the optimal X model
set.seed(2435)

params1A0 <- rnorm(20)
params1B0 <- rnorm(20)
params1Bs <- matrix(c(1, -0.65, 0.65, rep(0,7), rep(0.70, 10), 0, 1, 0.70, rep(-0.85, 7), rep(0, 10), 0, 0, 1, rep(0, 7), rep(-0.85, 10)), nrow = 20, ncol = 3)
params1As <- matrix(c(1, 0.45, -0.80, rep(0.80, 7), rep(0, 10), 0, 1, 0.70, rep(0, 7), rep(-0.95, 10)), nrow = 20, ncol = 2)
params1G <- matrix(c(0.60, -0.75, -0.70, 0.65, 0.80, -0.55), nrow = 2, ncol = 3)

q <- 20
p <- 20
s <- 3
m <- 2

params1Astotal <- cbind(params1A0, params1As)
params1Bstotal <- cbind(params1B0, params1Bs)

varvecB <- rep(1,q)
varvecA <- rep(1,p)
varvecG <- rep(1,m)

datset <- 1
sampsize <- 400

set.seed(datset)

dat <- datgen(params1Astotal, params1Bstotal, params1G, sampsize, varvecB, varvecA, varvecG)

Binit <- EMAlgB(dat$X, 3)


# R code only version
Bres1 <- EMAlgBAdLassoCV(dat$X, 3, Binit$B, Binit$Phi1, tuningpB = seq(0, 15, 1), weights = Binit$B)

#R/C++ version
Bres2 <- EMAlgBAdLassoCV2(dat$X, 3, Binit$B, Binit$Phi1, tuningpB = seq(0, 15, 1), weights = Binit$B)

