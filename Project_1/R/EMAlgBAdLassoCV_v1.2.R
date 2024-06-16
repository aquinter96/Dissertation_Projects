#' Penalized B EM Algorithm
#' 
#' An EM algorithm that can impose sparsity on B using coordinate descent.
#' The optimal tuning parameter for B is dtermined via grid search and BIC.
#' The estimates of B0, B, and Phi1 that correspond to the optimal tuning
#' parameter chosen are identified as the final estimates
#'
#' THIS CODE IS FOR THE VERSION OF THE X EM ALGORITHM WHERE ALL CHANGES WE MADE WERE
#' KEPT, WITH THE EXCEPTION OF REVERTING TO THE ORIGINAL METHOD OF ONLY ITERATING
#' THE COORDINATE DESCENT ALGORITHM ONCE PER M STEP

#' @param X The X dataset used to estimate B0, B, and Phi1. Must be the same
#' X used in EMAlgB
#' @param xk_sig The desired number of latent factors of X, corresponding
#' to the number of W variables and the number of columns of B. Must be the
#' same xk_sig used in EMAlgB
#' @param Binit The initial value of B estimated in EMAlgB
#' @param Phi1init The initial value of Phi1 estimated in EMAlgB
#' @param tuningpB A vector of numbers containing the candidate tuning
#' parameter values to be compared via grid search
#' @param nfolds Deprecated parameter
#' @param weights The qxs matrix containing the weights to be used in the
#' adaptive lasso coordinate descent update. Note that element (i,j) of
#' weights is the weight used for the coordinate descent update of
#' element (i,j) of B
#' 
#' @return Returns the final estimates of B0, B, and Phi1, which are the
#' the estimates obtained when using the optimal tuning parameter value
#' as determined by the grid search. Also returns the number of iterations
#' until convergence when using the optimal tuning parameter, the optimal
#' tuning parameter itself, the corresponding estimates of E(W), the total
#' computation time and time of the grid search, and the BIC value
#' corresponding to the optimal tuning parameter.
#' 
#' @examples
#' EMAlgBAdLassoCV(data$X, 3, Binit, Phi1init, seq(0, 5, 0.1), weights = Binit)
#' 
#' 
#' 
#' @export
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
      
      B0Cv <- rep(0, q)
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
          #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
          #to check for convergence
          B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
          Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                                      2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
          maxit <- maxit + 1
        }
        #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
        invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
        log.lik <- -(n/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
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