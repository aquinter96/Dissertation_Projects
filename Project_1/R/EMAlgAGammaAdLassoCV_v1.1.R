#' Penalized A EM Algorithm
#' 
#' An EM algorithm that can impose sparsity on A using coordinate descent.
#' The optimal tuning parameter for A is dtermined via grid search and BIC.
#' The estimates of A0, A, Gamma, Phi2, and Phi3 that correspond to the optimal tuning
#' parameter chosen are identified as the final estimates
#' 
#' #' THIS CODE IS FOR THE VERSION OF THE Y EM ALGORITHM WHERE ALL CHANGES WE MADE WERE
#' KEPT, WITH THE EXCEPTION OF REVERTING TO THE ORIGINAL METHOD OF ONLY CHECKING 
#' A FOR CONVERGENCE WHEN PERFORMING THE GRID SEARCH
#' 
#' @param X The X dataset
#' @param Y The Y dataset
#' @param xk_sig The number of latent factors for X
#' @param k_sig The number of latent factors for Y
#' @param B The final estimate of B from EMAlgBAdLassoCV
#' @param B0 The final estimate of B0 from EMAlgBAdLassoCV
#' @param Phi1 The final estimate of Phi1 from EMAlgBAdLassoCV
#' @param Ainit The initial estimate of A obtained from initvalcalc
#' @param Ginit The initial estimate of Gamma obtained from initvalcalc
#' @param Phi2init The inital estimate of Phi2 obtained from initvalcalc
#' @param Phi3init The inital estimate of Phi3 obtained from initvalcalc
#' @param tuningpA A vector of numbers containing the candidate tuning
#' parameter values to be compared via grid search
#' @param nfolds Deprecated parameter
#' @param weights The qxs matrix containing the weights to be used in the
#' adaptive lasso coordinate descent update. Note that element (i,j) of
#' weights is the weight used for the coordinate descent update of
#' element (i,j) of A. By default, the weights are the unpenalized MLE
#' estimates
#' 
#' @return Returns the final estimates of A0, A, Gamma, Phi2, and Phi3, which are the
#' the estimates obtained when using the optimal tuning parameter value
#' as determined by the grid search. Also returns the number of iterations
#' until convergence when using the optimal tuning parameter, the optimal
#' tuning parameter itself, the corresponding estimates of E(Z), the total
#' computation time and time of the grid search, and the BIC value
#' corresponding to the optimal tuning parameter
#' 
#' @examples
#' EMAlgAGammaAdLassoCV(data$X, data$Y, 3, 2, Best, B0est, Phi1est, Ainit, Ginit,
#' Phi2init, Phi3init, seq(0, 5, 0.1), weights = matrix(0, nrow = p, ncol = k_sig))
#' 
#' 
#' 
#' @export
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
  
  tuningpAinit <- tuningpA
  
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
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  
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
      ACV <- matrix(0, nrow = p, ncol = k_sig)
      Astart <- A
      
      BICmat <- rep(0, length(tuningpA))
      
      cv.start <- Sys.time()
      for(a in 1:length(tuningpA)){
        
        ACV <- Astart
        
        A <- matrix(0, nrow = p, ncol = k_sig)
        
        maxit <- 0
        
        while((norm(A - ACV) > 0.0001) & (maxit < 15000)){
          
          A <- ACV
          
          sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
          
          sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                         (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
          
          sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
          
          inv11 <- Matrix::solve(sigma11)
          condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
          condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
          
          EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
          EW <-  matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
          EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
          
          ACVupdate <- matrix(0, nrow = p, ncol = k_sig)
          
          while(norm(ACV - ACVupdate, type = "F") > 0.001){
            
            ACVupdate <- ACV
            
            for(i in 1:k_sig){
              for(j in 1:k_sig){
                deltaijA <- 0
                thetaijA <- 0
                omegaijA <- 0
                #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
                if(j < i){
                  deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
                  thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
                  omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
                  abar <- ((thetaijA-omegaijA))/(deltaijA)
                  if(weights[i,j] == 0){
                    lambdaA <- (tuningpA[a]*Phi2[i,i])/deltaijA*(1/1e-5)
                  }
                  else{
                    lambdaA <- (tuningpA[a]*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
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
                thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
                omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
                abar <- ((thetaijA-omegaijA))/(deltaijA)
                if(weights[i,j] == 0){
                  lambdaA <- (tuningpA[a]*Phi2[i,i])/deltaijA*(1/1e-5)
                }
                else{
                  lambdaA <- (tuningpA[a]*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
                }
                ACV[i,j] <- lasso(abar, lambdaA)
              }
            }
          }
        maxit <- maxit + 1
      }
      
      sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(ACV)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(ACV)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
      
      sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), ACV%*%Gamma%*%t(B)), rbind(B%*%t(ACV%*%Gamma),
                                                                                       (Phi2 + ACV%*%(Phi3 + tcrossprod(Gamma))%*%t(ACV)))), nrow = q+p, ncol = q+p)
      
      sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
      
      inv11 <- Matrix::solve(sigma11)
      condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
      condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
      
      log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
      
      BICmat[a] <- log(n)*(sum(ACV != 0) - k_sig) - 2*log.lik
      
    }
    if(length(BICmat) > 1){
      sigmaNAs <- is.na(BICmat)
      tuningpA <- tuningpA[!sigmaNAs]
      BICmat <- BICmat[!sigmaNAs]
    }
    if(length(BICmat) == 0){
      if(length(tuningpAinit) == 1){
        tuningpAhat <- tuningpAinit
      }
      if(length(tuningpAinit) > 1){
        tuningpAhat <- min(tuningpAinit)
      }
    }else{
      tuningpAhat <- tuningpA[which(BICmat == min(BICmat))]
    }
    
    A <- Astart
    
    #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
    #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
    if(length(tuningpAhat) > 1){
      tuningpAhat <- min(tuningpAhat)
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
  
  Anewupdate <- matrix(0, nrow = p, ncol = k_sig)
  
  while(norm(Anew - Anewupdate) > 0.001){
    
    Anewupdate <- Anew
    
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
