args <- commandArgs(TRUE)
datset <- as.numeric(args[1]) 
sampsize <- as.numeric(args[2])

library(parallel)

#Function for generating data for given A, A0, B, B0, Gamma, and error variances
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

#Function for fitting the Y model of JSFRM. Also performs CV to select the optimal tuning for parameter for adaptive LASSO for estimation of the B matrix
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
  
  while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (niter < 15000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
    if(niter == 0){
      
      cv.start <- Sys.time()
      
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
      
      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
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


#Function for fitting the Y model of JSFRM. Also performs CV to select the optimal tuning for parameter for adaptive LASSO for estimation of the A matrix
EMAlgAGammaAdLassoCV <- function(X, Y, xk_sig, k_sig, B, B0, Phi1, Ainit, Ginit, Phi2init, Phi3init, tuningpA = seq(0.1,10,0.1), nfolds=10, weights){
  
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  
  if(xk_sig == 1){
    B <- as.matrix(B)
    Phi3init <- as.matrix(Phi3init)
  }
  
  A0 <- matrix(0, nrow = p, ncol = 1, byrow = T) #use PC for initial estimates of A and B
  A0new <-  matrix(colMeans(Y), nrow = p, ncol = 1, byrow = T)
  A <- matrix(0, nrow = p, ncol = k_sig)
  Anew <- Ainit
  
  Gamma <- matrix(0, nrow = k_sig, ncol = xk_sig)
  Gammanew <- Ginit
  
  Phi2 <- matrix(0, nrow = p, ncol = p)
  Phi2new <- Phi2init
  Phi3 <- matrix(0, nrow = k_sig, ncol = k_sig)
  Phi3new <- Phi3init
  
  while(((norm(A0 - A0new) > 0.00001) | (norm(A - Anew) > 0.00001) | (norm(Gamma - Gammanew) > 0.00001) | (norm(Phi2 - Phi2new) > 0.00001) | (norm(Phi3 - Phi3new) > 0.00001)) & (niter < 15000)){
    
    A0 <- A0new
    A <- Anew
    Gamma <- Gammanew
    Phi2 <- Phi2new
    Phi3 <- Phi3new
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
    if(niter == 0){
      
      deltaijA <- 0 #initialize coordinate descent starting values
      thetaijA <- 0
      omegaijA <- 0
      lambdaA <- NULL
      abar <- NULL
      A0CV <- rep(0, p)
      ACV <- matrix(0, nrow = p, ncol = k_sig)
      GammaCV <- matrix(0, nrow = k_sig, ncol = xk_sig)
      Phi2CV <- matrix(0, nrow = p, ncol = p)
      Phi3CV <- matrix(0, nrow = k_sig, ncol = k_sig)
      A0start <- A0
      Astart <- A
      Gammastart <- Gamma
      Phi2start <- Phi2
      Phi3start <- Phi3
      
      cv.start <- Sys.time()

      BICs <- parSapply(cl, tuningpA, tunesearchA, X, B0, B, Y, A0, A, Gamma, Phi1, Phi2, Phi3, weights)

      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      if(length(BICs) > 1){
        sigmaNAs <- is.na(BICs)
        tuningpA <- tuningpA[!sigmaNAs]
        BICs <- BICs[!sigmaNAs]
      }
      if(length(BICs) == 0){
        if(length(tuningpAinit) == 1){
          tuningpAhat <- tuningpAinit
        }
        if(length(tuningpAinit) > 1){
          tuningpAhat <- min(tuningpAinit)
        }
      }else{
        tuningpAhat <- tuningpA[which(BICs == min(BICs))]
      }
      
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
    
    sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                   (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
    
    sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
    
    inv11 <- Matrix::solve(sigma11)
    condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
    condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m
    EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
    EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
    EW <-  matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
    
    #now that the optimal tuning parameter has been chosen, run the EM algorithm as normal, using the optimal tuning parameter for all subsequent iterations
    
    for(i in 1:k_sig){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        if(j < i){
          deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
          thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
          omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%Anew[i,-j]
          abar <- ((thetaijA-omegaijA))/(deltaijA)
          if(weights[i,j] == 0){
            lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
          }
          Anew[i,j] <- lasso(abar, lambdaA)
        }
      }
    }
    for(i in (k_sig+1):p){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%Anew[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
        }
        Anew[i,j] <- lasso(abar, lambdaA)
      }
    }
    
    #A0 Update
    A0new <- (1/n)*as.matrix(colSums(Y - t(Anew%*%EZ)))
    
    #Phi2 Update
    Phi2new <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0new, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(Anew) + matrix(A0new, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0new, nrow = p, ncol = n, byrow = F)) + 
                                 2*matrix(A0new, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(Anew) + n*Anew%*%condvar%*%t(Anew) + Anew%*%EZ%*%t(EZ)%*%t(Anew)))
    
    #conditional covariance matrix of W and Z
    condvarWZ <- matrix(cbind(rbind(diag(xk_sig), Gamma), rbind(t(Gamma),(Phi3 + Gamma%*%t(Gamma)))), nrow = xk_sig+k_sig, ncol = xk_sig+k_sig) -
      matrix(cbind(rbind(t(B), Gamma%*%t(B)), rbind(t(Gamma)%*%t(Anew),(Phi3 + Gamma%*%t(Gamma))%*%t(Anew))), nrow = xk_sig+k_sig, ncol = q+p)%*%solve(matrix(cbind(rbind((Phi1 + tcrossprod(B)), Anew%*%Gamma%*%t(B)), rbind(B%*%t(Anew%*%Gamma),
                                                                                                                                                                                                                              (Phi2new + Anew%*%(Phi3 + tcrossprod(Gamma))%*%t(Anew)))), nrow = q+p, ncol = q+p))%*%t(matrix(cbind(rbind(t(B), Gamma%*%t(B)), rbind(t(Gamma)%*%t(Anew),(Phi3 + Gamma%*%t(Gamma))%*%t(Anew))), nrow = xk_sig+k_sig, ncol = q+p))
    
    #Gamma Update
    Gammanew <- as.matrix((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:xk_sig),(1:xk_sig)] + EW%*%t(EW)))
    
    #Phi3 Update
    if(k_sig == 1){
      Phi3new <- as.matrix((1/n)*(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(Gammanew) -
                                    2*EZ%*%t(EW)%*%t(Gammanew) + n*Gammanew%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(Gammanew) + Gammanew%*%EW%*%t(EW)%*%t(Gammanew)))
    }else{
      Phi3new <- as.matrix((1/n)*diag(diag((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(Gammanew) -
                                              2*EZ%*%t(EW)%*%t(Gammanew) + n*Gammanew%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(Gammanew) + Gammanew%*%EW%*%t(EW)%*%t(Gammanew)))))
    }
    niter <- niter + 1
  }
  
  A <- Anew
  A0 <- A0new
  Gamma <- Gammanew
  Phi2 <- Phi2new
  Phi3 <- Phi3new
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                 (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
  
  inv11 <- Matrix::solve(sigma11)
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  BICopt <- log(n)*(sum(A != 0) - k_sig) - 2*log.lik
  output <- list("A" = A, "A0" = A0, "Gamma" = Gamma, "Phi2" = Phi2, "Phi3" = Phi3, "iterations" = niter, "optimal lambda" = tuningpAhat, "Y" = Y, "EZ" = EZ, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
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
OverallAGAlg <- function(X, Y, Bobj,tuningpA = seq(0,15,0.1), m_seq = 1:4, nfolds=10){
  Ainit <- NA
  s <- Bobj$`optimal s`
  B <- Bobj$Bopt$B
  B0 <- Bobj$Bopt$B0
  Phi1 <- Bobj$Bopt$Phi1
  Alist <- replicate(length(m_seq),NA,simplify=F)
  BIClist <- rep(0, length(m_seq))
  for(i in 1:length(m_seq)){
    inits <- initvalcalc(X, Y, B, s, m_seq[i], nrow(Y))
    Aweight <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, 0, weights = matrix(rep(0,ncol(Y)*m_seq[i]),ncol=m_seq[i]))
    Alist[[i]] <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, tuningpA, weights = Aweight$A)
    BIClist[i] <- Alist[[i]]$BICopt
  }
  mopt <- m_seq[which(BIClist == min(BIClist))]
  optA <- Alist[[which(BIClist == min(BIClist))]]
  return(list("model results" = optA, "X" = X, "B" = B, "B0" = B0, "Phi1" = Phi1, "optimal s BIC" = s, "optimal m BIC" = mopt, "ABIC" = optA$BICopt, "lambdaA" = optA$`optimal lambda`, "BBIC" = Bobj$BIC, "lambdaB" = Bobj$lambda))
}

set.seed(2435)

params1A0 <- rnorm(100)
params1B0 <- rnorm(100)
params1Bs <- matrix(0, nrow = 100, ncol = 3)
diag(params1Bs[(1:3),(1:3)]) <- 1
params1Bs[(1:3),(1:3)][lower.tri(params1Bs[(1:3),(1:3)], diag = F)] <- rbinom(3,1,0.5)*(2*rbinom(3,1,0.5)-1)*runif(3, 0.35, 1)
params1Bs[4:100,] <- rbinom(291,1,0.5)*(2*rbinom(291,1,0.5)-1)*runif(291, 0.35, 1)
params1As <- matrix(0, nrow = 100, ncol = 2)
diag(params1As[(1:2),(1:2)]) <- 1
params1As[(1:2),(1:2)][lower.tri(params1As[(1:2),(1:2)], diag = F)] <- rbinom(1,1,0.5)*(2*rbinom(1,1,0.5)-1)*runif(1, 0.35, 1)
params1As[3:100,] <- rbinom(196,1,0.5)*(2*rbinom(196,1,0.5)-1)*runif(196, 0.35, 1)
params1G <- matrix(c(0.60, -0.75, -0.70, 0.65, 0.80, -0.55), nrow = 2, ncol = 3)
q <- 100
p <- 100
s <- 3
m <- 2

params1Astotal <- cbind(params1A0, params1As)
params1Bstotal <- cbind(params1B0, params1Bs)

varvecB <- rep(1,q)
varvecA <- rep(1,p)
varvecG <- rep(1,m)

set.seed(datset)

dat <- datgen(params1Astotal, params1Bstotal, params1G, sampsize, varvecB, varvecA, varvecG)

lasso <- function(x,y){ #initialize soft-thresholding function
  result <- NULL
  if(abs(x) <= y){
    result <- 0
  } else{
    result <- x - y*sign(x)
  }
  return(result)
}

tunesearchB <- function(tuningpB, X, B0, B, Phi1, weights){
  
  B0CV <- B0
  BCV <- B
  Phi1CV <- Phi1
  
  q <- nrow(B)
  xk_sig <- ncol(B)
  n <- nrow(X)
  
  B0 <- rep(0, q)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Phi1 <- matrix(0, nrow = q, ncol = q)
  
  deltaijB <- 0
  thetaijB <- 0
  omegaijB <- 0
  
  maxit <- 0
  
  while(( (norm(B0 - B0CV, type = "F") > 0.00001) | (norm(B - BCV, type = "F") > 0.00001) | (norm(Phi1 - Phi1CV, type = "F") > 0.00001) ) & (maxit < 15000)){
    
    B0 <- B0CV
    B <- BCV
    Phi1 <- Phi1CV
    
    invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - crossprod(B,invmodcv)%*%B
    EW <- crossprod(B,invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
    
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
            lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/1e-5)
          }
          else{
            lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
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
          lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
        }
        BCV[i,j] <- lasso(bbar, lambdaB)
      }
    }
    #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
    #to check for convergence
    B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
    Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                                2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
    maxit <- maxit + 1
  }
  
  invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
  log.lik <- -(n/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1CV +tcrossprod(BCV))))
  
  BICs <- log(n)*(sum(BCV != 0) - xk_sig) - 2*log.lik
  
  BICs
  
}

tunesearchA <- function(tuningpA, X, B0, B, Y, A0, A, Gamma, Phi1, Phi2, Phi3, weights){
  
  ACV <- A
  A0CV <- A0
  Phi2CV <- Phi2
  Phi3CV <- Phi3
  GammaCV <- Gamma
  
  q <- nrow(B)
  p <- nrow(A)
  xk_sig <- ncol(Gamma)
  k_sig <- nrow(Gamma)
  n <- nrow(Y)
  
  A <- matrix(0, nrow = p, ncol = k_sig)
  A0 <- rep(0, p)
  Gamma <- matrix(0, nrow = k_sig, ncol = xk_sig)
  Phi2 <- matrix(0, nrow = p, ncol = p)
  Phi3 <- matrix(0, nrow = k_sig, ncol = k_sig)
  
  deltaijA <- 0
  thetaijA <- 0
  omegaijA <- 0
  
  maxit <- 0
  
  while(((norm(A - ACV) > 0.0001) | (norm(A0 - A0CV) > 0.0001) | (norm(Gamma - GammaCV) > 0.0001) | (norm(Phi2 - Phi2CV) > 0.0001) | (norm(Phi3 - Phi3CV) > 0.0001)) & (maxit < 15000)){
    
    A0 <-A0CV
    A <- ACV
    Gamma <- GammaCV
    Phi2 <- Phi2CV
    Phi3 <- Phi3CV
    
    sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
    
    sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                   (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
    
    sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
    
    inv11 <- Matrix::solve(sigma11)
    condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
    condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
    
    EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
    EW <- matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
    EZ <- matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
    
    for(i in 1:k_sig){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
        if(j < i){
          deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
          thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0CV[i,], n)))
          omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
          abar <- ((thetaijA-omegaijA))/(deltaijA)
          if(weights[i,j] == 0){
            lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
          }
          ACV[i,j] <- lasso(abar, lambdaA)
        }
      }
    }
    for(i in (k_sig+1):p){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0CV[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
        }
        ACV[i,j] <- lasso(abar, lambdaA)
      }
    }
    
    print(ACV)
    #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
    #to check for convergence
    A0CV <- (1/n)*as.matrix(colSums(Y - t(ACV%*%EZ)))
    
    Phi2CV <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(ACV) + matrix(A0, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0, nrow = p, ncol = n, byrow = F)) + 
                                2*matrix(A0, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(ACV) + n*ACV%*%condvar%*%t(ACV) + ACV%*%EZ%*%t(EZ)%*%t(ACV)))
    
    GammaCV <- as.matrix((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:xk_sig),(1:xk_sig)] + EW%*%t(EW)))
    
    #Phi3 Update
    if(k_sig == 1){
      Phi3CV <- (1/n)*as.matrix(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(GammaCV) -
                                  2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))
    }else{
      Phi3CV <- (1/n)*as.matrix(diag(diag(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(GammaCV) -
                                            2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))))
    }
    maxit <- maxit + 1
  }
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(GammaCV)%*% t(ACV)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(GammaCV%*%t(B),(Phi3 + GammaCV%*%t(GammaCV))%*%t(ACV)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = GammaCV W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), ACV%*%GammaCV%*%t(B)), rbind(B%*%t(ACV%*%GammaCV),
                                                                                     (Phi2CV + ACV%*%(Phi3 + tcrossprod(GammaCV))%*%t(ACV)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(GammaCV) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( GammaCV  , Phi3 + GammaCV%*%t(GammaCV) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
  
  inv11 <- Matrix::solve(sigma11)
  condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
  condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
  
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  BICs <- log(n)*(sum(ACV != 0) - k_sig) - 2*log.lik
  
  BICs
  
}

tottime <- system.time({
  
  cl <- makeCluster(4)

  clusterExport(cl, c("lasso", "tunesearchB", "tunesearchA"))
  
  Bres <- OverallBAlg(dat$X)
  
  Results <- OverallAGAlg(dat$X, dat$Y, Bres)
  
  stopCluster(cl)
  
})

saveRDS(list(Results, tottime), file=paste("parapars100rep", datset, "nod4cpu1n", sampsize, ".rds", sep=""))