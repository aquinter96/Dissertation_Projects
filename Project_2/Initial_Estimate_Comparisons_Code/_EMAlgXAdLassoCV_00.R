## _EMAlgXAdLassoCV_00.R

EMAlgXAdLassoCV <- function(Xpanel = list(), x_joint, x_indiv = c(), Ainits = list(), Binits = list(), vcovinits = list(), 
                            tuningpA = seq(5.1,5.5,0.1), tuningpB = seq(5.1,5.5,0.1), weightsA = list(), weightsB = list()){
  alg.start <- Sys.time()
  
  #define initial values for A1, A2, B1, B2, PhiR, PhiS1, PhiS2, Phi11, and Phi12
  niter <- 0
  n <- nrow(Xpanel[[1]])
  q1 <- ncol(Xpanel[[1]])
  q2 <- ncol(Xpanel[[2]])
  m <- x_joint
  u1 <- x_indiv[1]
  u2 <- x_indiv[2]
  A1 <- matrix(0, nrow = q1, ncol = m)
  A2 <- matrix(0, nrow = q2, ncol = m)
  B1 <- matrix(0, nrow = q1, ncol = u1)
  B2 <- matrix(0, nrow = q2, ncol = u2)
  PhiR <- matrix(0, nrow = m, ncol = m)
  PhiS1 <- matrix(0, nrow = u1, ncol = u1)
  PhiS2 <- matrix(0, nrow = u2, ncol = u2)
  Phi11 <- matrix(0, nrow = q1, ncol = q1)
  Phi12 <- matrix(0, nrow = q2, ncol = q2)
  A1new <- Ainits[[1]]
  A2new <- Ainits[[2]]
  B1new <- Binits[[1]]
  B2new <- Binits[[2]]
  PhiRnew <- vcovinits[[1]]
  PhiS1new <- vcovinits[[2]]
  PhiS2new <- vcovinits[[3]]
  Phi11new <- vcovinits[[4]]
  Phi12new <- vcovinits[[5]]
  
  ############################################################################ 
  ## assume sparsity for A1/A2//B1/B2, and thus, consider LASSO penalty  
  ## Here, we define a lasso function: a soft-thresholding function
  ############################################################################ 
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(((norm(A1 - A1new, type = "F") > 0.001) | (norm(A2 - A2new, type = "F") > 0.001) | (norm(B1 - B1new, type = "F") > 0.001) | (norm(B2 - B2new, type = "F") > 0.001) | (norm(PhiR - PhiRnew, type = "F") > 0.001) | (norm(PhiS1 - PhiS1new, type = "F") > 0.001) | (norm(PhiS2 - PhiS2new, type = "F") > 0.001) | (norm(Phi11 - Phi11new, type = "F") > 0.001) | (norm(Phi12 - Phi12new, type = "F") > 0.001)) & (niter < 10000)){
    A1 <- A1new
    A2 <- A2new
    B1 <- B1new
    B2 <- B2new
    PhiR <- PhiRnew
    PhiS1 <- PhiS1new
    PhiS2 <- PhiS2new
    Phi11 <- Phi11new
    Phi12 <- Phi12new
    
    ## Here grid search is used to find optimal tuning parameters for A1/A2 and B1/B2 in the LASSO part. 
    #for the first iteration of the algorithm ONLY, use grid search to calculate the optimal tuning parameter combination
    #for B and A
    
    if(niter == 0){
      
      epsilonijA1 <- 0 #initialize coordinate descent starting values
      thetaijA1 <- 0
      omegaijA1 <- 0
      tauijA1 <- 0
      epsilonijA2 <- 0
      thetaijA2 <- 0
      omegaijA2 <- 0
      tauijA2 <- 0
      epsilonijB1 <- 0
      thetaijB1 <- 0
      omegaijB1 <- 0
      tauijB1 <- 0
      epsilonijB2 <- 0
      thetaijB2 <- 0
      omegaijB2 <- 0
      tauijB2 <- 0
      lambdaA <- NULL
      lambdaB <- NULL
      a1bar <- NULL
      a2bar <- NULL
      b1bar <- NULL
      b2bar <- NULL
      A1CV <- matrix(0, nrow = q1, ncol = m)
      A2CV <- matrix(0, nrow = q2, ncol = m)
      B1CV <- matrix(0, nrow = q1, ncol = u1)
      B2CV <- matrix(0, nrow = q2, ncol = u2)
      PhiRCV <- matrix(0, nrow = m, ncol = m)
      PhiS1CV <- matrix(0, nrow = u1, ncol = u1)
      PhiS2CV <- matrix(0, nrow = u2, ncol = u2)
      Phi11CV <- matrix(0, nrow = q1, ncol = q1)
      Phi12CV <- matrix(0, nrow = q2, ncol = q2)
      A1start <- A1
      A2start <- A2
      B1start <- B1
      B2start <- B2
      PhiRstart <- PhiR
      PhiS1start <- PhiS1
      PhiS2start <- PhiS2
      Phi11start <- Phi11
      Phi12start <- Phi12
      maxit <- 0
      
      BICmat <- matrix(0, nrow = length(tuningpA), ncol = length(tuningpB))
      
      cv.start <- Sys.time()
      
      
      ## tuningpA = seq(5.1,5.5,0.1), tuningpB = seq(5.1,5.5,0.1),
      ## This means 5 * 5 = 25 sets of tuning parameters
      
      for(a in 1:length(tuningpA)){
        for(b in 1:length(tuningpB)){
          
          A1CV <- A1start
          A2CV <- A2start
          B1CV <- B1start
          B2CV <- B2start
          PhiRCV <- PhiRstart
          PhiS1CV <- PhiS1start
          PhiS2CV <- PhiS2start
          Phi11CV <- Phi11start
          Phi12CV <- Phi12start
          
          
          A1 <- matrix(0, nrow = q1, ncol = m)
          A2 <- matrix(0, nrow = q2, ncol = m)
          B1 <- matrix(0, nrow = q1, ncol = u1)
          B2 <- matrix(0, nrow = q2, ncol = u2)
          PhiR <- matrix(0, nrow = m, ncol = m)
          PhiS1 <- matrix(0, nrow = u1, ncol = u1)
          PhiS2 <- matrix(0, nrow = u2, ncol = u2)
          Phi11 <- matrix(0, nrow = q1, ncol = q1)
          Phi12 <- matrix(0, nrow = q2, ncol = q2)
          
          maxit <- 0
          
          ##########################################
          ## This while loop is for calculating estimates of A1/A2/B1/B2/PhiR/PhiS1/PhiS2/Phi11/Phi12
          ## for the current tuning parameter values a (for A) and b (for B)
          ##########################################
          while(((norm(A1 - A1CV, type = "F") > 0.001) | (norm(A2 - A2CV, type = "F") > 0.001) | (norm(B1 - B1CV, type = "F") > 0.001) | (norm(B2 - B2CV, type = "F") > 0.001) | (norm(PhiR - PhiRCV, type = "F") > 0.001) | (norm(PhiS1 - PhiS1CV, type = "F") > 0.001) | (norm(PhiS2 - PhiS2CV, type = "F") > 0.001) | (norm(Phi11 - Phi11CV, type = "F") > 0.001) | (norm(Phi12 - Phi12CV, type = "F") > 0.001)) & (maxit <= 10000)){
            
            A1 <- A1CV
            A2 <- A2CV
            B1 <- B1CV
            B2 <- B2CV
            PhiR <- PhiRCV
            PhiS1 <- PhiS1CV
            PhiS2 <- PhiS2CV
            Phi11 <- Phi11CV
            Phi12 <- Phi12CV
            
            sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2)), cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2)), cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2)))
            
            sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2)), cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12))
            
            sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2)),
                             cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2)),
                             cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2))
            
            condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
            condvarR <- condvar[(1:m), (1:m)]
            condvarS1 <- condvar[(m+(1:u1)), (m+(1:u1))]
            condvarS2 <- condvar[(m+u1+(1:u2)), (m+u1+(1:u2))]
            condvarRS1 <- condvar[(1:m), (m+(1:u1))]
            condvarRS2 <- condvar[(1:m), (m+u1+(1:u2))]
            
            ERS <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]))
            ER <- ERS[1:m,]
            ES1 <- ERS[m+(1:u1),]
            ES2 <- ERS[m+u1+(1:u2),]
            
            A1update <- matrix(0, nrow = q1, ncol = m)
            A2update <- matrix(0, nrow = q2, ncol = m)
            B1update <- matrix(0, nrow = q1, ncol = u1)
            B2update <- matrix(0, nrow = q2, ncol = u2)
            
            ##########################################
            ## This loop updates A1/A2/B1/B2 until convergence using coordinate descent for the current M-step
            ##########################################
            while((norm(A1update - A1CV, type = "F") > 0.005) | (norm(A2update - A2CV, type = "F") > 0.005) | (norm(B1update - B1CV) > 0.005) | (norm(B2update - B2CV) > 0.005)){
              
              A1update <- A1CV
              A2update <- A2CV
              B1update <- B1CV
              B2update <- B2CV
              
              ##################################              
              ## for each element in the first part of A1: mxm  
              ## A1 is q1*m
              ##################################              
              for(i in 1:m){
                for(j in 1:m){
                  epsilonijA1 <- 0
                  thetaijA1 <- 0
                  omegaijA1 <- 0
                  tauijA1 <- 0
                  
                  ## coordinate descent updates 
                  ## calculate coordinate descent updates for the first mxm submatrix of A, while ensuring that said mxm matrix is lower triangular
                  if(j < i){
                    epsilonijA1 <- (tcrossprod(ER) + n*condvarR)[j,j]
                    thetaijA1 <- crossprod(ER[j,],Xpanel[[1]][,i])
                    omegaijA1 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A1CV[i,-j]
                    for(k in 1:u1){
                      tauijA1 <- tauijA1 + t(ER%*%t(ES1) + n*condvarRS1)[k,j]%*%B1CV[i,k]
                    }
                    a1bar <- ((thetaijA1-omegaijA1-tauijA1))/(epsilonijA1)
                    if(weightsA[[1]][i,j] == 0){
                      lambdaA <- (tuningpA[a]*Phi11CV[i,i])/epsilonijA1*(1/1e-5)
                    }
                    else{
                      lambdaA <- (tuningpA[a]*Phi11CV[i,i])/epsilonijA1*(1/abs(weightsA[[1]][i,j]))
                    }
                    A1CV[i,j] <- lasso(a1bar, lambdaA)
                  }
                }
              }
              
              #calculate coordinate descent updates for remaining (q1-(m))xm submatrix of A1 
              for(i in (m+1):q1){
                for(j in 1:m){
                  epsilonijA1 <- 0
                  thetaijA1 <- 0
                  omegaijA1 <- 0
                  tauijA1 <- 0
                  epsilonijA1 <- (tcrossprod(ER) + n*condvarR)[j,j]
                  thetaijA1 <- crossprod(ER[j,],Xpanel[[1]][,i])
                  omegaijA1 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A1CV[i,-j]
                  for(k in 1:u1){
                    tauijA1 <- tauijA1 + t(ER%*%t(ES1) + n*condvarRS1)[k,j]%*%B1CV[i,k]
                  }
                  a1bar <- ((thetaijA1-omegaijA1-tauijA1))/(epsilonijA1)
                  if(weightsA[[1]][i,j] == 0){
                    lambdaA <- (tuningpA[a]*Phi11CV[i,i])/epsilonijA1*(1/1e-5)
                  }
                  else{
                    lambdaA <- (tuningpA[a]*Phi11CV[i,i])/epsilonijA1*(1/abs(weightsA[[1]][i,j]))
                  }
                  A1CV[i,j] <- lasso(a1bar, lambdaA)
                }
              }
              
              #calculate coordinate descent updates for the (q2)xm submatrix of A2 
              for(i in 1:q2){
                for(j in 1:m){
                  epsilonijA2 <- 0
                  thetaijA2 <- 0
                  omegaijA2 <- 0
                  tauijA2 <- 0
                  epsilonijA2 <- (tcrossprod(ER) + n*condvarR)[j,j]
                  thetaijA2 <- crossprod(ER[j,],Xpanel[[2]][,i])
                  omegaijA2 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A2CV[i,-j]
                  for(k in 1:u2){
                    tauijA2 <- tauijA2 + t(ER%*%t(ES2) + n*condvarRS2)[k,j]%*%B2CV[i,k]
                  }
                  a2bar <- ((thetaijA2-omegaijA2-tauijA2))/(epsilonijA2)
                  if(weightsA[[2]][i,j] == 0){
                    lambdaA <- (tuningpA[a]*Phi12CV[i,i])/epsilonijA2*(1/1e-5)
                  }
                  else{
                    lambdaA <- (tuningpA[a]*Phi12CV[i,i])/epsilonijA2*(1/abs(weightsA[[2]][i,j]))
                  }
                  A2CV[i,j] <- lasso(a2bar, lambdaA)
                }
              }
              
              # calculate coordinate descent updates for the (u1)xu1 submatrix of B1 
              # B1 is q1*u1 matrix
              for(i in 1:u1){
                for(j in 1:u1){
                  epsilonijB1 <- 0
                  thetaijB1 <- 0
                  omegaijB1 <- 0
                  tauijB1 <- 0
                  if(j < i){
                    epsilonijB1 <- (tcrossprod(ES1) + n*condvarS1)[j,j]
                    thetaijB1 <- crossprod(ES1[j,],Xpanel[[1]][,i])
                    omegaijB1 <- (tcrossprod(ES1) + n*condvarS1)[-j,j]%*%B1CV[i,-j]
                    for(k in 1:m){
                      tauijB1 <- tauijB1 + t(ES1%*%t(ER) + n*t(condvarRS1))[k,j]%*%A1CV[i,k]
                    }
                    b1bar <- ((thetaijB1-omegaijB1-tauijB1))/(epsilonijB1)
                    if(weightsB[[1]][i,j] == 0){
                      lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB1*(1/1e-5)
                    }
                    else{
                      lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB1*(1/abs(weightsB[[1]][i,j]))
                    }
                    B1CV[i,j] <- lasso(b1bar, lambdaB)
                  }
                }
              }
              
              #calculate coordinate descent updates for remaining (q1-(u1))*u1 submatrix of B1
              for(i in (u1+1):q1){
                for(j in 1:u1){
                  epsilonijB1 <- 0
                  thetaijB1 <- 0
                  omegaijB1 <- 0
                  tauijB1 <- 0
                  epsilonijB1 <- (tcrossprod(ES1) + n*condvarS1)[j,j]
                  thetaijB1 <- crossprod(ES1[j,],Xpanel[[1]][,i])
                  omegaijB1 <- (tcrossprod(ES1) + n*condvarS1)[-j,j]%*%B1CV[i,-j]
                  for(k in 1:m){
                    tauijB1 <- tauijB1 + t(ES1%*%t(ER) + n*t(condvarRS1))[k,j]%*%A1CV[i,k]
                  }
                  b1bar <- ((thetaijB1-omegaijB1-tauijB1))/(epsilonijB1)
                  if(weightsB[[1]][i,j] == 0){
                    lambdaB <- (tuningpB[b]*Phi11CV[i,i])/epsilonijB1*(1/1e-5)
                  }
                  else{
                    lambdaB <- (tuningpB[b]*Phi11CV[i,i])/epsilonijB1*(1/abs(weightsB[[1]][i,j]))
                  }
                  B1CV[i,j] <- lasso(b1bar, lambdaB)
                }
              }
              
              # calculate coordinate descent updates for the (u2)*u2 submatrix of B2
              # B1 is q2*u2 matrix
              for(i in 1:u2){
                for(j in 1:u2){
                  epsilonijB2 <- 0
                  thetaijB2 <- 0
                  omegaijB2 <- 0
                  tauijB2 <- 0
                  if(j < i){
                    epsilonijB2 <- (tcrossprod(ES2) + n*condvarS2)[j,j]
                    thetaijB2 <- crossprod(ES2[j,],Xpanel[[2]][,i])
                    omegaijB2 <- (tcrossprod(ES2) + n*condvarS2)[-j,j]%*%B2CV[i,-j]
                    for(k in 1:m){
                      tauijB2 <- tauijB2 + t(ES2%*%t(ER) + n*t(condvarRS2))[k,j]%*%A2CV[i,k]
                    }
                    b2bar <- ((thetaijB2-omegaijB2-tauijB2))/(epsilonijB2)
                    if(weightsB[[2]][i,j] == 0){
                      lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB2*(1/1e-5)
                    }
                    else{
                      lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB2*(1/abs(weightsB[[2]][i,j]))
                    }
                    B2CV[i,j] <- lasso(b2bar, lambdaB)
                  }
                }
              }
              
              #calculate coordinate descent updates for remaining (q2-(u2))*u2 submatrix of B2
              for(i in (u2+1):q2){
                for(j in 1:u2){
                  epsilonijB2 <- 0
                  thetaijB2 <- 0
                  omegaijB2 <- 0
                  tauijB2 <- 0
                  epsilonijB2 <- (tcrossprod(ES2) + n*condvarS2)[j,j]
                  thetaijB2 <- crossprod(ES2[j,],Xpanel[[2]][,i])
                  omegaijB2 <- (tcrossprod(ES2) + n*condvarS2)[-j,j]%*%B2CV[i,-j]
                  for(k in 1:m){
                    tauijB2 <- tauijB2 + t(ES2%*%t(ER) + n*t(condvarRS2))[k,j]%*%A2CV[i,k]
                  }
                  b2bar <- ((thetaijB2-omegaijB2-tauijB2))/(epsilonijB2)
                  if(weightsB[[2]][i,j] == 0){
                    lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB2*(1/1e-5)
                  }
                  else{
                    lambdaB <- (tuningpB[b]*Phi12CV[i,i])/epsilonijB2*(1/abs(weightsB[[2]][i,j]))
                  }
                  B2CV[i,j] <- lasso(b2bar, lambdaB)
                }
              }
            }
            
            #################################################################################
            ## back to the while loop to update PhiR/PhiS1/PhiS2/Phi11/Phi12
            #################################################################################
            
            PhiRCV <- diag(diag(condvarR + 1/n*ER%*%t(ER)))
            PhiS1CV <- diag(diag(condvarS1 + 1/n*ES1%*%t(ES1)))
            PhiS2CV <- diag(diag(condvarS2 + 1/n*ES2%*%t(ES2)))
            Phi11CV <- 1/n*diag(diag(t(Xpanel[[1]])%*%Xpanel[[1]] - 2*t(Xpanel[[1]])%*%t(A1CV%*%ER) - 2*t(Xpanel[[1]])%*%t(B1CV%*%ES1) + 2*n*A1CV%*%condvarRS1%*%t(B1CV)
                                     + 2*A1CV%*%ER%*%t(B1%*%ES1) + n*A1CV%*%condvarR%*%t(A1CV) + A1CV%*%ER%*%t(A1CV%*%ER) + n*B1CV%*%condvarS1%*%t(B1CV) + B1CV%*%ES1%*%t(B1%*%ES1)))
            Phi12CV <- 1/n*diag(diag(t(Xpanel[[2]])%*%Xpanel[[2]] - 2*t(Xpanel[[2]])%*%t(A2CV%*%ER) - 2*t(Xpanel[[2]])%*%t(B2CV%*%ES2) + 2*n*A2CV%*%condvarRS2%*%t(B2CV)
                                     + 2*A2CV%*%ER%*%t(B2%*%ES2) + n*A2CV%*%condvarR%*%t(A2CV) + A2CV%*%ER%*%t(A2CV%*%ER) + n*B2CV%*%condvarS2%*%t(B2CV) + B2CV%*%ES2%*%t(B2%*%ES2)))
            maxit <- maxit + 1
          }
          
          sigma11 <- rbind(cbind(A1CV%*%PhiRCV%*%t(A1CV) + B1CV%*%PhiS1CV%*%t(B1CV) + Phi11CV, A1CV%*%PhiRCV%*%t(A2CV)), cbind(A2CV%*%PhiRCV%*%t(A1CV), A2CV%*%PhiRCV%*%t(A2CV) + B2CV%*%PhiS2CV%*%t(B2CV) + Phi12CV))
          
          ## once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
          log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]],Xpanel[[2]]))+(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))))) + (q1+q2)*log((2*pi))+log(det(sigma11)))
          
          ## obtain the BIC under this set of tuning parameter (tuningpA[a], tuningpB[b])
          BICmat[a,b] <- log(n)*(sum(A1CV != 0) + sum(A2CV != 0) + sum(B1CV != 0) + sum(B2CV != 0) - (m+u1+u2)) - 2*log.lik
        }
      }
      
      ## To get the optimal tuning parameter ... 
      tuningpAhat <- tuningpA[which.min(apply(BICmat, MARGIN = 1, min))]
      tuningpBhat <- tuningpB[which.min(apply(BICmat, MARGIN = 2, min))]
      
      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
      
      A1 <- A1start
      A2 <- A2start
      B1 <- B1start
      B2 <- B2start
      PhiR <- PhiRstart
      PhiS1 <- PhiS1start
      PhiS2 <- PhiS2start
      Phi11 <- Phi11start
      Phi12 <- Phi12start
      
      print("CV done")
    }
    
    sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2)), cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2)), cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2)))
    
    sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2)), cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12))
    
    sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2)),
                     cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2)),
                     cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2))
    
    condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
    condvarR <- condvar[(1:m), (1:m)]
    condvarS1 <- condvar[(m+(1:u1)), (m+(1:u1))]
    condvarS2 <- condvar[(m+u1+(1:u2)), (m+u1+(1:u2))]
    condvarRS1 <- condvar[(1:m), (m+(1:u1))]
    condvarRS2 <- condvar[(1:m), (m+u1+(1:u2))]
    
    ERS <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]))
    ER <- ERS[1:m,]
    ES1 <- ERS[m+(1:u1),]
    ES2 <- ERS[m+u1+(1:u2),]
    
    #now that the optimal tuning parameter has been chosen, run the EM algorithm as normal, using the optimal tuning parameter for all subsequent iterations
    epsilonijA1 <- 0 #initialize coordinate descent starting values
    thetaijA1 <- 0
    omegaijA1 <- 0
    tauijA1 <- 0
    epsilonijA2 <- 0
    thetaijA2 <- 0
    omegaijA2 <- 0
    tauijA2 <- 0
    epsilonijB1 <- 0
    thetaijB1 <- 0
    omegaijB1 <- 0
    tauijB1 <- 0
    epsilonijB2 <- 0
    thetaijB2 <- 0
    omegaijB2 <- 0
    tauijB2 <- 0
    lambdaA <- NULL
    lambdaB <- NULL
    a1bar <- NULL
    a2bar <- NULL
    b1bar <- NULL
    b2bar <- NULL
    
    A1update <- matrix(0, nrow = q1, ncol = m)
    A2update <- matrix(0, nrow = q2, ncol = m)
    B1update <- matrix(0, nrow = q1, ncol = u1)
    B2update <- matrix(0, nrow = q2, ncol = u2)
    
    #############################################################################
    ## Here we get A1/A2/B1/B2 under the optimal tuning parameter 
    ## This part replicates the above code. Can be wrapped into a function ... 
    #############################################################################
    
    while((norm(A1update - A1new, type = "F") > 0.005) | (norm(A2update - A2new, type = "F") > 0.005) | (norm(B1update - B1new) > 0.005) | (norm(B2update - B2new) > 0.005)){
      
      A1update <- A1new
      A2update <- A2new
      B1update <- B1new
      B2update <- B2new
      
      for(i in 1:m){
        for(j in 1:m){
          epsilonijA1 <- 0
          thetaijA1 <- 0
          omegaijA1 <- 0
          tauijA1 <- 0
          #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
          if(j < i){
            epsilonijA1 <- (tcrossprod(ER) + n*condvarR)[j,j]
            thetaijA1 <- crossprod(ER[j,],Xpanel[[1]][,i])
            omegaijA1 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A1new[i,-j]
            for(k in 1:u1){
              tauijA1 <- tauijA1 + t(ER%*%t(ES1) + n*condvarRS1)[k,j]%*%B1new[i,k]
            }
            a1bar <- ((thetaijA1-omegaijA1-tauijA1))/(epsilonijA1)
            if(weightsA[[1]][i,j] == 0){
              lambdaA <- (tuningpAhat*Phi11new[i,i])/epsilonijA1*(1/1e-5)
            }
            else{
              lambdaA <- (tuningpAhat*Phi11new[i,i])/epsilonijA1*(1/abs(weightsA[[1]][i,j]))
            }
            A1new[i,j] <- lasso(a1bar, lambdaA)
          }
        }
      }
      
      for(i in (m+1):q1){
        for(j in 1:m){
          epsilonijA1 <- 0
          thetaijA1 <- 0
          omegaijA1 <- 0
          tauijA1 <- 0
          epsilonijA1 <- (tcrossprod(ER) + n*condvarR)[j,j]
          thetaijA1 <- crossprod(ER[j,],Xpanel[[1]][,i])
          omegaijA1 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A1new[i,-j]
          for(k in 1:u1){
            tauijA1 <- tauijA1 + t(ER%*%t(ES1) + n*condvarRS1)[k,j]%*%B1new[i,k]
          }
          a1bar <- ((thetaijA1-omegaijA1-tauijA1))/(epsilonijA1)
          if(weightsA[[1]][i,j] == 0){
            lambdaA <- (tuningpAhat*Phi11new[i,i])/epsilonijA1*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpAhat*Phi11new[i,i])/epsilonijA1*(1/abs(weightsA[[1]][i,j]))
          }
          A1new[i,j] <- lasso(a1bar, lambdaA)
        }
      }
      
      for(i in 1:q2){
        for(j in 1:m){
          epsilonijA2 <- 0
          thetaijA2 <- 0
          omegaijA2 <- 0
          tauijA2 <- 0
          epsilonijA2 <- (tcrossprod(ER) + n*condvarR)[j,j]
          thetaijA2 <- crossprod(ER[j,],Xpanel[[2]][,i])
          omegaijA2 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A2new[i,-j]
          for(k in 1:u2){
            tauijA2 <- tauijA2 + t(ER%*%t(ES2) + n*condvarRS2)[k,j]%*%B2new[i,k]
          }
          a2bar <- ((thetaijA2-omegaijA2-tauijA2))/(epsilonijA2)
          if(weightsA[[2]][i,j] == 0){
            lambdaA <- (tuningpAhat*Phi12new[i,i])/epsilonijA2*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpAhat*Phi12new[i,i])/epsilonijA2*(1/abs(weightsA[[2]][i,j]))
          }
          A2new[i,j] <- lasso(a2bar, lambdaA)
        }
      }
      
      for(i in 1:u1){
        for(j in 1:u1){
          epsilonijB1 <- 0
          thetaijB1 <- 0
          omegaijB1 <- 0
          tauijB1 <- 0
          if(j < i){
            epsilonijB1 <- (tcrossprod(ES1) + n*condvarS1)[j,j]
            thetaijB1 <- crossprod(ES1[j,],Xpanel[[1]][,i])
            omegaijB1 <- (tcrossprod(ES1) + n*condvarS1)[-j,j]%*%B1new[i,-j]
            for(k in 1:m){
              tauijB1 <- tauijB1 + t(ES1%*%t(ER) + n*t(condvarRS1))[k,j]%*%A1new[i,k]
            }
            b1bar <- ((thetaijB1-omegaijB1-tauijB1))/(epsilonijB1)
            if(weightsB[[1]][i,j] == 0){
              lambdaB <- (tuningpBhat%*%Phi11new[i,i])/epsilonijB1*(1/1e-5)
            }
            else{
              lambdaB <- (tuningpBhat%*%Phi11new[i,i])/epsilonijB1*(1/abs(weightsB[[1]][i,j]))
            }
            B1new[i,j] <- lasso(b1bar, lambdaB)
          }
        }
      }
      #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
      for(i in (u1+1):q1){
        for(j in 1:u1){
          epsilonijB1 <- 0
          thetaijB1 <- 0
          omegaijB1 <- 0
          tauijB1 <- 0
          epsilonijB1 <- (tcrossprod(ES1) + n*condvarS1)[j,j]
          thetaijB1 <- crossprod(ES1[j,],Xpanel[[1]][,i])
          omegaijB1 <- (tcrossprod(ES1) + n*condvarS1)[-j,j]%*%B1new[i,-j]
          for(k in 1:m){
            tauijB1 <- tauijB1 + t(ES1%*%t(ER) + n*t(condvarRS1))[k,j]%*%A1new[i,k]
          }
          b1bar <- ((thetaijB1-omegaijB1-tauijB1))/(epsilonijB1)
          if(weightsB[[1]][i,j] == 0){
            lambdaB <- (tuningpBhat%*%Phi11new[i,i])/epsilonijB1*(1/1e-5)
          }
          else{
            lambdaB <- (tuningpBhat%*%Phi11new[i,i])/epsilonijB1*(1/abs(weightsB[[1]][i,j]))
          }
          B1new[i,j] <- lasso(b1bar, lambdaB)
        }
      }
      
      for(i in 1:u2){
        for(j in 1:u2){
          epsilonijB2 <- 0
          thetaijB2 <- 0
          omegaijB2 <- 0
          tauijB2 <- 0
          if(j < i){
            epsilonijB2 <- (tcrossprod(ES2) + n*condvarS2)[j,j]
            thetaijB2 <- crossprod(ES2[j,],Xpanel[[2]][,i])
            omegaijB2 <- (tcrossprod(ES2) + n*condvarS2)[-j,j]%*%B2new[i,-j]
            for(k in 1:m){
              tauijB2 <- tauijB2 + t(ES2%*%t(ER) + n*t(condvarRS2))[k,j]%*%A2new[i,k]
            }
            b2bar <- ((thetaijB2-omegaijB2-tauijB2))/(epsilonijB2)
            if(weightsB[[2]][i,j] == 0){
              lambdaB <- (tuningpBhat%*%Phi12new[i,i])/epsilonijB2*(1/1e-5)
            }
            else{
              lambdaB <- (tuningpBhat%*%Phi12new[i,i])/epsilonijB2*(1/abs(weightsB[[2]][i,j]))
            }
            B2new[i,j] <- lasso(b2bar, lambdaB)
          }
        }
      }
      
      for(i in (u2+1):q2){
        for(j in 1:u2){
          epsilonijB2 <- 0
          thetaijB2 <- 0
          omegaijB2 <- 0
          tauijB2 <- 0
          epsilonijB2 <- (tcrossprod(ES2) + n*condvarS2)[j,j]
          thetaijB2 <- crossprod(ES2[j,],Xpanel[[2]][,i])
          omegaijB2 <- (tcrossprod(ES2) + n*condvarS2)[-j,j]%*%B2new[i,-j]
          for(k in 1:m){
            tauijB2 <- tauijB2 + t(ES2%*%t(ER) + n*t(condvarRS2))[k,j]%*%A2new[i,k]
          }
          b2bar <- ((thetaijB2-omegaijB2-tauijB2))/(epsilonijB2)
          if(weightsB[[2]][i,j] == 0){
            lambdaB <- (tuningpBhat%*%Phi12new[i,i])/epsilonijB2*(1/1e-5)
          }
          else{
            lambdaB <- (tuningpBhat%*%Phi12new[i,i])/epsilonijB2*(1/abs(weightsB[[2]][i,j]))
          }
          B2new[i,j] <- lasso(b2bar, lambdaB)
        }
      }
    }
    
    PhiRnew <- diag(diag(condvarR + 1/n*ER%*%t(ER)))
    PhiS1new <- diag(diag(condvarS1 + 1/n*ES1%*%t(ES1)))
    PhiS2new <- diag(diag(condvarS2 + 1/n*ES2%*%t(ES2)))
    Phi11new <- 1/n*diag(diag(t(Xpanel[[1]])%*%Xpanel[[1]] - 2*t(Xpanel[[1]])%*%t(A1new%*%ER) - 2*t(Xpanel[[1]])%*%t(B1new%*%ES1) + 2*n*A1new%*%condvarRS1%*%t(B1new)
                              + 2*A1new%*%ER%*%t(B1%*%ES1) + n*A1new%*%condvarR%*%t(A1new) + A1new%*%ER%*%t(A1new%*%ER) + n*B1new%*%condvarS1%*%t(B1new) + B1new%*%ES1%*%t(B1%*%ES1)))
    Phi12new <- 1/n*diag(diag(t(Xpanel[[2]])%*%Xpanel[[2]] - 2*t(Xpanel[[2]])%*%t(A2new%*%ER) - 2*t(Xpanel[[2]])%*%t(B2new%*%ES2) + 2*n*A2new%*%condvarRS2%*%t(B2new)
                              + 2*A2new%*%ER%*%t(B2%*%ES2) + n*A2new%*%condvarR%*%t(A2new) + A2new%*%ER%*%t(A2new%*%ER) + n*B2new%*%condvarS2%*%t(B2new) + B2new%*%ES2%*%t(B2%*%ES2)))
    
    niter <- niter + 1
  }
  
  A1 <- A1new
  A2 <- A2new
  B1 <- B1new
  B2 <- B2new
  PhiR <- PhiRnew
  PhiS1 <- PhiS1new
  PhiS2 <- PhiS2new
  Phi11 <- Phi11new
  Phi12 <- Phi12new
  
  sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2)), cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2)), cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2)))
  
  sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2)), cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12))
  
  sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2)),
                   cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2)),
                   cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2))
  
  condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
  condvarR <- condvar[(1:m), (1:m)]
  condvarS1 <- condvar[(m+(1:u1)), (m+(1:u1))]
  condvarS2 <- condvar[(m+u1+(1:u2)), (m+u1+(1:u2))]
  condvarRS1 <- condvar[(1:m), (m+(1:u1))]
  condvarRS2 <- condvar[(1:m), (m+u1+(1:u2))]
  
  ERS <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]))
  ER <- ERS[1:m,]
  ES1 <- ERS[m+(1:u1),]
  ES2 <- ERS[m+u1+(1:u2),]
  
  log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]],Xpanel[[2]]))+(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))))) + (q1+q2)*log((2*pi))+log(det(sigma11)))
  
  BICopt <- log(n)*(sum(A1 != 0) + sum(A2 != 0) + sum(B1 != 0) + sum(B2 != 0) - (m+u1+u2)) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("A1" = A1, "A2" = A2, "B1" = B1, "B2" = B2, "PhiR" = PhiR, "PhiS1" = PhiS1, "PhiS2" = PhiS2, "Phi11" = Phi11, "Phi12" = Phi12, "iterations" = niter, "optimal lambdaA" = tuningpAhat, "optimal lambdaB" = tuningpBhat, "ERS" = ERS, "BICopt" = BICopt)
  return(output)
}
