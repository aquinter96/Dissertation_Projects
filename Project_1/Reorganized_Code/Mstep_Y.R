## Mstep_Y.R

Mstep_Y <- function(X, Y, EW, EZ, condvarWZ, condvar, A0, A, Gamma, Phi2, Phi3, tuningpA, weights){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  xk_sig <- ncol(Gamma)
  k_sig <- ncol(A)
  
  ACV <- A
  
  ##################################              
  ## for each element in the first part of A: mxm  
  ## S is p*m
  ##################################              
  for(i in 1:k_sig){
    for(j in 1:k_sig){
      deltaijA <- 0
      thetaijA <- 0
      omegaijA <- 0
      #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
      if(j < i){
        deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%A[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
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
      omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%A[i,-j]
      abar <- ((thetaijA-omegaijA))/(deltaijA)
      if(weights[i,j] == 0){
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
      }
      else{
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
      }
      ACV[i,j] <- lasso(abar, lambdaA)
    }
  }
  
  #################################################################################
  ## back to the while loop to update A0/Gamma/Phi2/Phi3
  #################################################################################
  
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
    
  ests <- list("A0" = A0CV, "A" = ACV, "Gamma" = GammaCV, "Phi2" = Phi2CV, "Phi3" = Phi3CV)
  
  return(ests)
  
}
  
## end of code