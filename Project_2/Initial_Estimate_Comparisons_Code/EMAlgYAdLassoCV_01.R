## _EMAlgYAdLassoCV_00.R
## _EMAlgYAdLassoCV_01.R  (to divide into multiple functions)
## EMAlgYAdLassoCV_01.R  (divided into multiple functions)

#########################################################################################################
## This function was used in both X-part and Y-part estimation. Will move it to another file later. 
#########################################################################################################
lasso <- function(x,y){ #initialize soft-thresholding function
  result <- NULL
  if(abs(x) <= y){
    result <- 0
  } else{
    result <- x - y*sign(x)
  }
  return(result)
}


#########################################################################################################
## E-step 
EMAlgYAdLasso_EStep <- function(Data = Data,  ## (X1/X2,Y1/Y2)
                                old_Model_Param = old_Model_Param) 
  ## Input: Data, Model_param
  ## Output: conditional covariance and ... 
{
  # calculate covariance of latent factors with observed data, i.e. Cov(Ei, Fi), where
  # E = [Ri, Si1, Si2, Vi, Wi1, Wi2] and F = [Xi1, Xi2, Yi1, Yi2]
  
  ## to hold all E-step results needed for M-step  
  Estep <- list()
  
  MyPar <- old_Model_Param
  ## get all meta parameters (dimensions)
  {
    m <- ncol(MyPar$PhiR)
    u1 <- ncol(MyPar$PhiS1)
    u2 <-  ncol(MyPar$PhiS2)   
    q1 <-  ncol(MyPar$Phi11)
    q2 <-  ncol(MyPar$Phi12)
    
    p1 <-  ncol(MyPar$Phi21)
    p2 <- ncol(MyPar$Phi22)
    j <- ncol(MyPar$Phi3)
    z1 <- ncol(MyPar$Phi41)
    z2 <- ncol(MyPar$Phi42)      
  }
  ## get all matrix parameters 
  {
    PhiR   <-   MyPar$PhiR
    PhiS1  <-   MyPar$PhiS1
    PhiS2  <-   MyPar$PhiS2
    Phi11  <-   MyPar$Phi11
    Phi12  <-   MyPar$Phi12
    A1     <-   MyPar$A1
    A2     <-   MyPar$A2
    B1     <-   MyPar$B1
    B2     <-   MyPar$B2
    
    Phi3   <-  MyPar$Phi3
    Phi41  <-  MyPar$Phi41
    Phi42  <-  MyPar$Phi42
    Phi21  <-  MyPar$Phi21
    Phi22  <-  MyPar$Phi22
    Theta  <-  MyPar$Theta
    Psi    <-  MyPar$Psi
    Pi1    <-  MyPar$Pi1
    Pi2    <-  MyPar$Pi2
    Omega1 <-  MyPar$Omega1
    Omega2 <-  MyPar$Omega2
    Gamma1 <-  MyPar$Gamma1
    Gamma2 <-  MyPar$Gamma2
    Delta1 <-  MyPar$Delta1
    Delta2 <-  MyPar$Delta2 
    
    PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
    PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
    PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))
  }
  
  Estep$sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2), PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1), PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2)),
                         cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS1S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS1S))),
                         cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS2S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS2S))),
                         cbind(t(A1%*%PhiR%*%t(Theta) + B1%*%PhiS1S%*%t(Psi)), t(A2%*%PhiR%*%t(Theta) + B2%*%PhiS2S%*%t(Psi)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Theta) + 
                                                                                                                                   (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Psi) + Gamma1%*%Phi3), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Theta) + 
                                                                                                                                                                                                          (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Psi) + Gamma2%*%Phi3)),
                         cbind(t(A1%*%PhiR%*%t(Pi1) + B1%*%PhiS1S%*%t(Omega1)), t(A2%*%PhiR%*%t(Pi1) + B2%*%PhiS2S%*%t(Omega1)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi1) + 
                                                                                                                                     (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega1) + Delta1%*%Phi41), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi1) + (Gamma2%*%Psi + 
                                                                                                                                                                                                                                                                   Delta2%*%Omega2)%*%PhiS%*%t(Omega1))),
                         cbind(t(A1%*%PhiR%*%t(Pi2) + B1%*%PhiS1S%*%t(Omega2)), t(A2%*%PhiR%*%t(Pi2) + B2%*%PhiS2S%*%t(Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi2) + 
                                                                                                                                     (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega2)), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega2) + 
                                                                                                                                                                                               Delta2%*%Phi42)))
  
  #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
  
  Estep$sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), 
                               A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                         cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), 
                               A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                         cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + 
                                                                                                                                      B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + 
                                 (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, 
                               (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + 
                                 Gamma1%*%Phi3%*%t(Gamma2)),
                         cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + 
                                                                                                                                      B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + 
                                                                                                                                                                                            (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), 
                               (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + 
                                 Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
  
  #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
  # 
  Estep$sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2), PhiR%*%t(Theta), PhiR%*%t(Pi1), PhiR%*%t(Pi2)),
                         cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2), PhiS1S%*%t(Psi), PhiS1S%*%t(Omega1), PhiS1S%*%t(Omega2)),
                         cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2, PhiS2S%*%t(Psi), PhiS2S%*%t(Omega1), PhiS2S%*%t(Omega2)),
                         cbind(Theta%*%PhiR, Psi%*%t(PhiS1S), Psi%*%t(PhiS2S), Theta%*%PhiR%*%t(Theta) + Psi%*%PhiS%*%t(Psi) + Phi3, Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1), 
                               Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)),
                         cbind(Pi1%*%PhiR, Omega1%*%t(PhiS1S), Omega1%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1)), Pi1%*%PhiR%*%t(Pi1) + Omega1%*%PhiS%*%t(Omega1) + 
                                 Phi41, Pi1%*%PhiR%*%t(Pi2) + Omega1%*%PhiS%*%t(Omega2)),
                         cbind(Pi2%*%PhiR, Omega2%*%t(PhiS1S), Omega2%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)), Pi2%*%PhiR%*%t(Pi1) + Omega2%*%PhiS%*%t(Omega1), 
                               Pi2%*%PhiR%*%t(Pi2) + Omega2%*%PhiS%*%t(Omega2) + Phi42))
  
  #calculate conditional covariance of latent factors given observed data, i.e. Cov(Ei|Fi)
  #and then extract the conditional covariances of each latent factor (Cov(Ri|Fi), Cov(Si1|Fi), etc)
  
  Estep$condvar <- Estep$sigma22 - Estep$sigma21%*%solve(Estep$sigma11)%*%t(Estep$sigma21)

  # calculate conditional expectations of latent factors given observed data
  # 
  # ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
  # Xpanel = list(t(MyData$X1), t(MyData$X2)), Ypanel = list(t(MyData$Y1), t(MyData$Y2))  
  Estep$ERSVW <- Estep$sigma21%*%solve(Estep$sigma11)%*%rbind( (Data$X1), (Data$X2), (Data$Y1), (Data$Y2))
  
  return(Estep)
}

##Example:  EStepRes <- EMAlgYAdLasso_EStep(Data = Data, old_Model_Param = old_Model_Param) 

#########################################################################################################
## M-step 
EMAlgYAdLasso_MStep <- function(Data = Data,  ## (X1/X2,Y1/Y2)
                                old_Model_Param = old_Model_Param, 
                                E_res = EStepRes,
                                tuningpGammaA = 0,
                                tuningpDeltaB = 0,
                                weightsGamma = list(old_Model_Param$Gamma1, old_Model_Param$Gamma2), 
                                weightsDelta = list(old_Model_Param$Delta1, old_Model_Param$Delta2) ) 
  ## Input: Data, Model_param, E_res
  ## Output: new_Model_para ... 
{
  
  ###########################################################################################
  #the coordinate descent update parameters for each matrix (epsilonijG1, thetaijG1, etc.)
  ###########################################################################################
  
  MyPar <- old_Model_Param
  ## The following should go to M-step ... 
  ## Meta parameters (dimensions)
  {
    m <- ncol(MyPar$PhiR)
    u1 <- ncol(MyPar$PhiS1)
    u2 <-  ncol(MyPar$PhiS2)   
    q1 <-  ncol(MyPar$Phi11)
    q2 <-  ncol(MyPar$Phi12)
    
    p1 <-  ncol(MyPar$Phi21)
    p2 <-   ncol(MyPar$Phi22)
    y_joint <- ncol(MyPar$Phi3) ## j 
    z1 <- ncol(MyPar$Phi41)
    z2 <- ncol(MyPar$Phi42) 
    
    n <- ncol(Data$X1) ## sample size 
  }
  
  ## conditional variances from E-step 
  condvar <- E_res$condvar
  ## The following should go to M-step ... 
  condvarR <- condvar[(1:m), (1:m)]
  condvarS <- condvar[(m+(1:(u1+u2))), (m+(1:(u1+u2)))]
  condvarV <- condvar[(m+u1+u2+(1:y_joint)), (m+u1+u2+(1:y_joint))]
  condvarW1 <- condvar[(m+u1+u2+y_joint+(1:z1)), (m+u1+u2+y_joint+(1:z1))]
  condvarW2 <- condvar[(m+u1+u2+y_joint+z1+(1:z2)), (m+u1+u2+y_joint+z1+(1:z2))]
  condvarRS <- condvar[(1:m), (m+(1:(u1+u2)))]
  condvarRW1 <- condvar[(1:m), (m+u1+u2+y_joint+(1:z1))]
  condvarRW2 <- condvar[(1:m), (m+u1+u2+y_joint+z1+(1:z2))]
  condvarSW1 <- condvar[(m+(1:(u1+u2))), (m+u1+u2+y_joint+(1:z1))]
  condvarSW2 <- condvar[(m+(1:(u1+u2))), (m+u1+u2+y_joint+z1+(1:z2))]
  condvarVR <- condvar[(m+u1+u2+(1:y_joint)), (1:m)]
  condvarVS <- condvar[(m+u1+u2+(1:y_joint)), (m+(1:(u1+u2)))]
  condvarVW1 <- condvar[(m+u1+u2+(1:y_joint)), (m+u1+u2+y_joint+(1:z1))]
  condvarVW2 <- condvar[(m+u1+u2+(1:y_joint)), (m+u1+u2+y_joint+z1+(1:z2))]
  
  ## conditional R/S/V/W from E-step 
  ERSVW <- E_res$ERSVW
  ## The following should go to M-step ... 
  ER <- ERSVW[1:m,]
  ES1 <- ERSVW[m+(1:u1),]
  ES2 <- ERSVW[m+u1+(1:u2),]
  ES <- ERSVW[m+(1:(u1+u2)),]
  EV <- ERSVW[m+u1+u2+(1:y_joint),]
  EW1 <- ERSVW[m+u1+u2+y_joint+(1:z1),]
  EW2 <- ERSVW[m+u1+u2+y_joint+z1+(1:z2),]

  {
    epsilonijG1 <- 0 #initialize coordinate descent starting values
    thetaijG1 <- 0
    omegaijG1 <- 0
    tauijG1 <- 0
    epsilonijG2 <- 0
    thetaijG2 <- 0
    omegaijG2 <- 0
    tauijG2 <- 0
    epsilonijD1 <- 0
    thetaijD1 <- 0
    omegaijD1 <- 0
    tauijD1 <- 0
    epsilonijD2 <- 0
    thetaijD2 <- 0
    omegaijD2 <- 0
    tauijD2 <- 0
    lambdaG <- NULL
    lambdaD <- NULL
    g1bar <- NULL
    g2bar <- NULL
    d1bar <- NULL
    d2bar <- NULL
  }  
  
  ## Main part of M-step: to update using GSD for ... 
  {
    ###########################################################            
    ## To update Gamma1/Gamma2/Delta1/Delta2 only 
    ###########################################################            
    
    ## while((norm(Gamma1update - Gamma1CV, type = "F") > 0.005) | (norm(Gamma2update - Gamma2CV, type = "F") > 0.005) | (norm(Delta1update - Delta1CV) > 0.005) | (norm(Delta2update - Delta2CV) > 0.005)){
      
      ## Is this problematic? 
      ##   Gamma1CV is updated when we update Gamma1CV? 
      Gamma1CV <- MyPar$Gamma1
      Gamma2CV <- MyPar$Gamma2
      Delta1CV <- MyPar$Delta1
      Delta2CV <- MyPar$Delta2
      Phi21CV <- MyPar$Phi21
      Phi22CV <- MyPar$Phi22
      
      ## why not using j instead of y_joint?
      for(i in 1:y_joint){
        for(j in 1:y_joint){
          #reset the coordinate descent parameters to 0 each time we start a new coordinate
          epsilonijG1 <- 0
          thetaijG1 <- 0
          omegaijG1 <- 0
          tauijG1 <- 0
          #calculate coordinate descent updates for the first jxj submatrix of Gamma, while ensuring that said sxs matrix is lower triangular
          if(j < i){
            epsilonijG1 <- (tcrossprod(EV) + n*condvarV)[j,j]
            thetaijG1 <- crossprod(EV[j,],t(Data$Y1)[,i])
            omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1CV[i,-j]
            for(k in 1:z1){
              tauijG1 <- tauijG1 + Delta1CV[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
            }
            g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
            if(weightsGamma[[1]][i,j] == 0){
              lambdaG <- (tuningpGammaA*Phi21CV[i,i])/epsilonijG1*(1/1e-5)
            }
            else{
              lambdaG <- (tuningpGammaA*Phi21CV[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
            }
            Gamma1CV[i,j] <- lasso(g1bar, lambdaG)
          }
        }
      }
      #calculate coordinate descent updates for remaining (p1-(j+1))xj submatrix of Gamma1
      for(i in (y_joint+1):p1){
        for(j in 1:y_joint){
          epsilonijG1 <- 0
          thetaijG1 <- 0
          omegaijG1 <- 0
          tauijG1 <- 0
          epsilonijG1 <- (tcrossprod(EV) + n*condvarV)[j,j]
          thetaijG1 <- crossprod(EV[j,],t(Data$Y1)[,i])
          omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1CV[i,-j]
          for(k in 1:z1){
            tauijG1 <- tauijG1 + Delta1CV[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
          }
          g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
          if(weightsGamma[[1]][i,j] == 0){
            lambdaG <- (tuningpGammaA*Phi21CV[i,i])/epsilonijG1*(1/1e-5)
          }
          else{
            lambdaG <- (tuningpGammaA*Phi21CV[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
          }
          Gamma1CV[i,j] <- lasso(g1bar, lambdaG)
        }
      }
      
      #coordiante descent for Gamma2
      for(i in 1:p2){
        for(j in 1:y_joint){
          epsilonijG2 <- 0
          thetaijG2 <- 0
          omegaijG2 <- 0
          tauijG2 <- 0
          epsilonijG2 <- (tcrossprod(EV) + n*condvarV)[j,j]
          thetaijG2 <- crossprod(EV[j,],t(Data$Y2)[,i])
          omegaijG2 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma2CV[i,-j]
          for(k in 1:z2){
            tauijG2 <- tauijG2 + Delta2CV[i,k]*t(EV%*%t(EW2) + n*condvarVW2)[k,j]
          }
          g2bar <- ((thetaijG2-omegaijG2-tauijG2))/(epsilonijG2)
          if(weightsGamma[[2]][i,j] == 0){
            lambdaG <- (tuningpGammaA*Phi22CV[i,i])/epsilonijG2*(1/1e-5)
          }
          else{
            lambdaG <-(tuningpGammaA*Phi22CV[i,i])/epsilonijG2*(1/abs(weightsGamma[[2]][i,j]))
          }
          Gamma2CV[i,j] <- lasso(g2bar, lambdaG)
        }
      }
      
      for(i in 1:z1){
        for(j in 1:z1){
          epsilonijD1 <- 0
          thetaijD1 <- 0
          omegaijD1 <- 0
          tauijD1 <- 0
          #calculate coordinate descent updates for the first z1xz1 submatrix of Delta1, while ensuring that said z1xz1 matrix is lower triangular
          if(j < i){
            epsilonijD1 <- (tcrossprod(EW1) + n*condvarW1)[j,j]
            thetaijD1 <- crossprod(EW1[j,],t(Data$Y1)[,i])
            omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1CV[i,-j]
            for(k in 1:y_joint){
              tauijD1 <- tauijD1 + Gamma1CV[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
            }
            d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
            if(weightsDelta[[1]][i,j] == 0){
              lambdaD <- (tuningpDeltaB*Phi21CV[i,i])/epsilonijD1*(1/1e-5)
            }
            else{
              lambdaD <- (tuningpDeltaB*Phi21CV[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
            }
            Delta1CV[i,j] <- lasso(d1bar, lambdaD)
          }
        }
      }
      #calculate coordinate descent updates for remaining (z1-(z1+1))xz1 submatrix of Delta1
      for(i in (z1+1):p1){
        for(j in 1:z1){
          epsilonijD1 <- 0
          thetaijD1 <- 0
          omegaijD1 <- 0
          tauijD1 <- 0
          epsilonijD1 <- (tcrossprod(EW1) + n*condvarW1)[j,j]
          thetaijD1 <- crossprod(EW1[j,],t(Data$Y1)[,i])
          omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1CV[i,-j]
          for(k in 1:y_joint){
            tauijD1 <- tauijD1 + Gamma1CV[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
          }
          d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
          if(weightsDelta[[1]][i,j] == 0){
            lambdaD <- (tuningpDeltaB*Phi21CV[i,i])/epsilonijD1*(1/1e-5)
          }
          else{
            lambdaD <- (tuningpDeltaB*Phi21CV[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
          }
          Delta1CV[i,j] <- lasso(d1bar, lambdaD)
        }
      }
      
      for(i in 1:z2){
        for(j in 1:z2){
          epsilonijD2 <- 0
          thetaijD2 <- 0
          omegaijD2 <- 0
          tauijD2 <- 0
          #calculate coordinate descent updates for the first z2xz2 submatrix of Delta2, while ensuring that said z2xz2 matrix is lower triangular
          if(j < i){
            epsilonijD2 <- (tcrossprod(EW2) + n*condvarW2)[j,j]
            thetaijD2 <- crossprod(EW2[j,],t(Data$Y2)[,i])
            omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2CV[i,-j]
            for(k in 1:y_joint){
              tauijD2 <- tauijD2 + Gamma2CV[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
            }
            d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
            if(weightsDelta[[2]][i,j] == 0){
              lambdaD <- (tuningpDeltaB*Phi22CV[i,i])/epsilonijD2*(1/1e-5)
            }
            else{
              lambdaD <- (tuningpDeltaB*Phi22CV[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
            }
            Delta2CV[i,j] <- lasso(d2bar, lambdaD)
          }
        }
      }
      
      #calculate coordinate descent updates for remaining (z2-(z2+1))xz2 submatrix of Delta2
      for(i in (z2+1):p2){
        for(j in 1:z2){
          epsilonijD2 <- 0
          thetaijD2 <- 0
          omegaijD2 <- 0
          tauijD2 <- 0
          epsilonijD2 <- (tcrossprod(EW2) + n*condvarW2)[j,j]
          thetaijD2 <- crossprod(EW2[j,],t(Data$Y2)[,i])
          omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2CV[i,-j]
          for(k in 1:y_joint){
            tauijD2 <- tauijD2 + Gamma2CV[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
          }
          d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
          if(weightsDelta[[2]][i,j] == 0){
            lambdaD <- (tuningpDeltaB*Phi22CV[i,i])/epsilonijD2*(1/1e-5)
          }
          else{
            lambdaD <- (tuningpDeltaB*Phi22CV[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
          }
          Delta2CV[i,j] <- lasso(d2bar, lambdaD)
        }
      }
    }
    
    ############################################################################ 
    ## This returns to the first while loop in the Cross-Validation part.           
    ############################################################################ 
    # Psi <- MyPar$Psi 
  
    ThetaCV <- (n*condvarVR + EV%*%t(ER) - n*(MyPar$Psi)%*%t(condvarRS) - (MyPar$Psi)%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    Pi1CV <- (n*t(condvarRW1) + EW1%*%t(ER) - n*MyPar$Omega1%*%t(condvarRS) - MyPar$Omega1%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    Pi2CV <- (n*t(condvarRW2) + EW2%*%t(ER) - n*MyPar$Omega2%*%t(condvarRS) - MyPar$Omega2%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    PsiCV <- (n*condvarVS + EV%*%t(ES) - n*ThetaCV%*%condvarRS - ThetaCV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    Omega1CV <- (n*t(condvarSW1) + EW1%*%t(ES) - n*Pi1CV%*%condvarRS - Pi1CV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    Omega2CV <- (n*t(condvarSW2) + EW2%*%t(ES) - n*Pi2CV%*%condvarRS - Pi2CV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    
    Phi21CV <- 1/n*diag(diag((Data$Y1)%*%t(Data$Y1) - 2*(Data$Y1)%*%t(Gamma1CV%*%EV) - 2*(Data$Y1)%*%t(Delta1CV%*%EW1) + 2*n*Gamma1CV%*%condvarVW1%*%t(Delta1CV)
                             + 2*Gamma1CV%*%EV%*%t(Delta1CV%*%EW1) + n*Gamma1CV%*%condvarV%*%t(Gamma1CV) + Gamma1CV%*%EV%*%t(Gamma1CV%*%EV) + n*Delta1CV%*%condvarW1%*%t(Delta1CV) + Delta1CV%*%EW1%*%t(Delta1CV%*%EW1)))
    Phi22CV <- 1/n*diag(diag((Data$Y2)%*%t(Data$Y2) - 2*(Data$Y2)%*%t(Gamma2CV%*%EV) - 2*(Data$Y2)%*%t(Delta2CV%*%EW2) + 2*n*Gamma2CV%*%condvarVW2%*%t(Delta2CV)
                             + 2*Gamma2CV%*%EV%*%t(Delta2CV%*%EW2) + n*Gamma2CV%*%condvarV%*%t(Gamma2CV) + Gamma2CV%*%EV%*%t(Gamma2CV%*%EV) + n*Delta2CV%*%condvarW2%*%t(Delta2CV) + Delta2CV%*%EW2%*%t(Delta2CV%*%EW2)))
    Phi3CV <- 1/n*diag(diag(n*condvarV + EV%*%t(EV) - 2*n*condvarVR%*%t(ThetaCV) - 2*EV%*%t(ThetaCV%*%ER) - 2*n*condvarVS%*%t(PsiCV) - 2*EV%*%t(PsiCV%*%ES) + n*ThetaCV%*%condvarR%*%t(ThetaCV) + ThetaCV%*%ER%*%t(ThetaCV%*%ER)
                            + 2*n*ThetaCV%*%condvarRS%*%t(PsiCV) + 2*ThetaCV%*%ER%*%t(PsiCV%*%ES) + n*PsiCV%*%condvarS%*%t(PsiCV) + PsiCV%*%ES%*%t(PsiCV%*%ES)))
    Phi41CV <- 1/n*diag(diag(n*condvarW1 + EW1%*%t(EW1) - 2*n*t(Pi1CV%*%condvarRW1) -2*EW1%*%t(Pi1CV%*%ER) -2*n*t(Omega1CV%*%condvarSW1) - 2*EW1%*%t(Omega1CV%*%ES) + n*Pi1CV%*%condvarR%*%t(Pi1CV) + Pi1CV%*%ER%*%t(Pi1CV%*%ER)
                             + 2*n*Pi1CV%*%condvarRS%*%t(Omega1CV) + 2*Pi1CV%*%ER%*%t(Omega1CV%*%ES) + n*Omega1CV%*%condvarS%*%t(Omega1CV) + Omega1CV%*%ES%*%t(Omega1CV%*%ES)))
    Phi42CV <- 1/n*diag(diag(n*condvarW2 + EW2%*%t(EW2) - 2*n*t(Pi2CV%*%condvarRW2) -2*EW2%*%t(Pi2CV%*%ER) -2*n*t(Omega2CV%*%condvarSW2) - 2*EW2%*%t(Omega2CV%*%ES) + n*Pi2CV%*%condvarR%*%t(Pi2CV) + Pi2CV%*%ER%*%t(Pi2CV%*%ER)
                             + 2*n*Pi2CV%*%condvarRS%*%t(Omega2CV) + 2*Pi2CV%*%ER%*%t(Omega2CV%*%ES) + n*Omega2CV%*%condvarS%*%t(Omega2CV) + Omega2CV%*%ES%*%t(Omega2CV%*%ES)))
    
    Model_Params <- MyPar 
    # Model_Params$PhiR     <- PhiRCV
    # Model_Params$PhiS1    <- PhiS1CV
    # Model_Params$PhiS2    <- PhiS2CV
    # Model_Params$Phi11    <- Phi11CV
    # Model_Params$Phi12    <- Phi12CV
    # Model_Params$A1       <- A1CV
    # Model_Params$A2       <- A2CV
    # Model_Params$B1       <- B1CV
    # Model_Params$B2       <- B2CV
    
    Model_Params$Phi3      <- Phi3CV
    Model_Params$Phi41     <- Phi41CV
    Model_Params$Phi42     <- Phi42CV
    Model_Params$Phi21     <- Phi21CV
    Model_Params$Phi22     <- Phi22CV
    Model_Params$Theta     <- ThetaCV
    Model_Params$Psi       <- PsiCV
    Model_Params$Pi1       <- Pi1CV
    Model_Params$Pi2       <- Pi2CV
    Model_Params$Omega1    <- Omega1CV
    Model_Params$Omega2    <- Omega2CV
    Model_Params$Gamma1    <- Gamma1CV
    Model_Params$Gamma2    <- Gamma2CV
    Model_Params$Delta1    <- Delta1CV
    Model_Params$Delta2    <- Delta2CV  
    
    ## 
    return(Model_Params)
  }
## Example: MStepRes <- EMAlgYAdLasso_MStep(Data = Data, old_Model_Param = old_Model_Param, E_res = EStepRes, tuningpGammaA = 0, tuningpDeltaB = 0) 

#########################################################################################################
## to calculate log likelihood 
## may need to add penalty part later 
#########################################################################################################
logLik <- function(Data = Data, E_res = EStepRes){
  ## to caclulate log-likelihood given current parameters   
  p1 <- nrow(Data$Y1)
  p2 <- nrow(Data$Y2)
  n <-  ncol(Data$Y2)
  
  log.lik <- -(n/2)*(sum(diag(solve(E_res$sigma11)%*%(cov(cbind(t(Data$X1), t(Data$X2), t(Data$Y1),t(Data$Y2)))+(colMeans(cbind(t(Data$X1), t(Data$X2), t(Data$Y1), t(Data$Y2))) - 
                     colMeans(cbind(t(Data$X1), t(Data$X2), t(Data$Y1), t(Data$Y2))))%*%t(colMeans(cbind(t(Data$X1), t(Data$X2), t(Data$Y1), t(Data$Y2))) - 
                     colMeans(cbind(t(Data$X1), t(Data$X2), t(Data$Y1), t(Data$Y2))))))) + (p1+p2)*log((2*pi))+log(det(E_res$sigma11)))
  return(log.lik)  
}
## Example   log.lik <- logLik(Data = Data, E_res = EStepRes) 

#########################################################################################################
## To check change in parameters 
##  may add more criteria if needed 
#########################################################################################################
# 
paraDiff <- function(old_Model_Param, new_Model_Param){
##   
## To calculate change in parameters before and after updates   
##  
  myDiffL0 <- myDiffL1 <- NULL  
  for (ii in 1:length(old_Model_Param)){
    myDiffL0 <- c(myDiffL0,   max( abs( as.vector(old_Model_Param[[ii]] - new_Model_Param[[ii]]) ) ) ) ## max One can calculate norml
    myDiffL1 <- c(myDiffL1,   norm(old_Model_Param[[ii]] - new_Model_Param[[ii]], type = "F") ) ## 
  }

  return( list(maxAbs = myDiffL0, normDiff = myDiffL1) )  
}
## Example  myDiff <- paraDiff(old_Model_Param, MStepRes)


#########################################################################################################
## The MAIN function for EM-algorithm
##  To just do the Y-part EM iterations given tuning parameters 
#################################################################################
## The EM algorithm function with one set of (tuningpGamma, tuningpDelta)
#################################################################################
EMAlgYAdLasso <-         function(Data = MyData,  ## (X1/X2,Y1/Y2)
                                 initial_Model_Param = Model_Params, ## all matrix model parameters
                                 tuningpGamma = 0, 
                                 tuningpDelta = 0, 
                                 weightsGamma = list(), 
                                 weightsDelta = list()){

  ## for code development and test purpose 
  ## Data = MyData; initial_Model_Param = Model_Params; tuningpGamma = 0; tuningpDelta = 0; 
  ## weightsGamma <- list(Model_Params$Gamma1, Model_Params$Gamma2);  weightsDelta = list(Model_Params$Delta1, Model_Params$Delta2)
  alg.start <- Sys.time()

  logLikList <- NULL 
  myDiffList <- NULL 
  
  ## define new_Model_Param to hold updated model parameters 
  old_Model_Param <- initial_Model_Param   
  
  ## The following EM-update (one iteration) should be a function ... 
  niter <- 0
  Continue_status <- 1  ## decide if continue the loop 
  
  # while(((norm(Gamma1 - Gamma1new) > 0.001) | (norm(Gamma2 - Gamma2new) > 0.001) | (norm(Delta1 - Delta1new) > 0.001) | (norm(Delta2 - Delta2new) > 0.001) | 
  #        (norm(Theta - Thetanew) > 0.001) | (norm(Pi1 - Pi1new, type = "F") > 0.001) | (norm(Psi - Psinew) > 0.001) | (norm(Omega1 - Omega1new) > 0.001) | 
  #        (norm(Omega2 - Omega2new) > 0.001) | (norm(Phi21 - Phi21new) > 0.001) | (norm(Phi22 - Phi22new) > 0.001) | (norm(Phi3 - Phi3new) > 0.001) |  
  #        (norm(Phi41 - Phi41new) > 0.001) | (norm(Phi42 - Phi42new) > 0.001)) & (niter < 30000)){
    
  while (Continue_status ){
    niter <- niter + 1 
    ## E-step 
    EStepRes <- EMAlgYAdLasso_EStep(Data = Data, old_Model_Param = old_Model_Param) 
    ## M-step
    ## 
    ## Note: to fix weightsGamma at the true Gamma / Delta might not be a good idea. We can discuss this. 
    ## 
    MStepRes <- EMAlgYAdLasso_MStep(Data = Data, old_Model_Param = old_Model_Param, E_res = EStepRes, tuningpGammaA = 0, tuningpDeltaB = 0,
                                  weightsGamma = weightsGamma, 
                                  weightsDelta = weightsDelta)

    
    ## get log.likelihood, and difference in parameters (max of norm of difference)  
    log.lik <- logLik(Data = Data, E_res = EStepRes) 
    myDiff <- paraDiff(old_Model_Param, MStepRes)
    # print(log.lik)
    
    ##############################################################################
    ## I feel that we should use change in log-likelihood as stopping criteria 
    ## 
    diff_criteria <-  max(myDiff$maxAbs)
    if (niter > 300 | diff_criteria< 0.001 ){
      Continue_status <- 0
    }
    ##############################################################################
    ## 
    ## To update Model parameters
    old_Model_Param <-  MStepRes
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria) 
  }
  
  time.diff <- Sys.time() - alg.start
  
  ## BIC, tuning parameters 
  ## To add  
  
  ## to save results 
  return( list(est_Model_param = old_Model_Param,  log.Like = logLikList, diffList = myDiffList, time.diff = time.diff) )
  
}

## Example:  
# myEMRes <-       EMAlgYAdLasso(Data = MyData,  
#                               initial_Model_Param = Model_Params, 
#                               tuningpGamma = 0, 
#                               tuningpDelta = 0, 
#                               weightsGamma = list(Model_Params$Gamma1, Model_Params$Gamma2), 
#                               weightsDelta = list(Model_Params$Delta1, Model_Params$Delta2) )

## 
##
## 
### End of code 
