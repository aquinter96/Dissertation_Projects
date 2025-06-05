Y_inits <- function(lats, Data, Xest){
  
  inits <- list()
  
  X1 <- Data$X1
  X2 <- Data$X2
  Y1 <- Data$Y1
  Y2 <- Data$Y2
  
  A1 <- Xest$X_final_pars$A1
  B1 <- Xest$X_final_pars$B1
  A2 <- Xest$X_final_pars$A2
  B2 <- Xest$X_final_pars$B2
  
  PhiR <- Xest$X_final_pars$PhiR
  PhiS1 <- Xest$X_final_pars$PhiS1
  PhiS2 <- Xest$X_final_pars$PhiS2
  
  j <- lats[[1]]
  z1 <- lats[[2]]
  z2 <- lats[[3]]
  p1 <- ncol(Y1)
  p2 <- ncol(Y2)
  m <- ncol(A1)
  u1 <- ncol(B1)
  u2 <- ncol(B2)
  u <- u1 + u2
  n <- nrow(Y1)
  
  PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
  PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
  
  C <- rbind(cbind(A1%*%PhiR, B1%*%PhiS1S), cbind(A2%*%PhiR, B2%*%PhiS2S))
  
  Ga <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(X1,Y1), cov(X2,Y1))
  Gb <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(X1,Y2), cov(X2,Y2))
  G1 <- t(Ga[1:m,])
  G2 <- t(Gb[1:m,])
  G3 <- t(Ga[(m+1):(m+u),])
  G4 <- t(Gb[(m+1):(m+u),])
  Ycov <- cov(cbind(Y1, Y2))
  D <- Ycov - rbind(t(Ga), t(Gb))%*%t(rbind(t(Ga), t(Gb)))
  
  G1G2svd <- svd(D[(1:p1), (p1+(1:p2))])
  G1init <- G1G2svd$u[,(1:j)]
  G2init <- G1G2svd$v[,(1:j)]
  if(j == 1){
    Gammainit <- c(G1init, G2init)
  }else{
    Gammainit <- rbind(G1init, G2init)
  }
  Gamma3 <- Gammainit%*%(eigen(1/n*crossprod(Gammainit))$vectors)
  Gamma5 <- Gamma3%*%qr.Q(qr(t(Gamma3[(1:j),(1:j)])))
  if(j == 1){
    phi3inh <- Gamma5[1]
  }else{
    phi3inh <- diag(diag(Gamma5[(1:j),(1:j)]))
  }
  Gammainit <- Gamma5 %*% solve(phi3inh)
  phi3in <- phi3inh%*%phi3inh
  if(j > 1){
    Gammainit[(1:j),(1:j)][upper.tri(Gammainit[(1:j),(1:j)], diag = F)] <- 0
  }
  if(j == 1){
    Gammainit[1] <- 1
  }else{
    diag(Gammainit[(1:j),(1:j)]) <- 1
  }
  Gamma1init <- Gammainit[(1:p1),]
  Gamma2init <- Gammainit[(p1+(1:p2)),]
  
  Delta1init <- svd(D[(1:p1), (1:p1)] - Gamma1init%*%phi3in%*%t(Gamma1init) - diag(p1))$u[,1:z1]
  Delta13 <- Delta1init%*%(eigen(1/n*crossprod(Delta1init))$vectors)
  Delta15 <- Delta13%*%qr.Q(qr(t(Delta13[(1:z1),(1:z1)])))
  if(z1 == 1){
    phi41inh <- Delta15[1]
  }else{
    phi41inh <- diag(diag(Delta15[(1:z1),(1:z1)]))
  }
  Delta1init <- Delta15 %*% solve(phi41inh)
  phi41in <- phi41inh%*%phi41inh
  if(z1 > 1){
    Delta1init[(1:z1),(1:z1)][upper.tri(Delta1init[(1:z1),(1:z1)], diag = F)] <- 0
  }
  if(z1 == 1){
    Delta1init[1] <- 1
  }else{
    diag(Delta1init[(1:z1),(1:z1)]) <- 1
  }
  
  Delta2init <- svd(D[(p1+(1:p2)), (p1+(1:p2))] - Gamma2init%*%phi3in%*%t(Gamma2init) - diag(p2))$u[,1:z2]
  Delta23 <- Delta2init%*%(eigen(1/n*crossprod(Delta2init))$vectors)
  Delta25 <- Delta23%*%qr.Q(qr(t(Delta23[(1:z2),(1:z2)])))
  if(z2 == 1){
    phi42inh <- Delta25[1]
  }else{
    phi42inh <- diag(diag(Delta25[(1:z2),(1:z2)]))
  }
  Delta2init <- Delta25 %*% solve(phi42inh)
  phi42in <- phi42inh%*%phi42inh
  if(z2 > 1){
    Delta2init[(1:z2),(1:z2)][upper.tri(Delta2init[(1:z2),(1:z2)], diag = F)] <- 0
  }
  if(z2 == 1){
    Delta2init[1] <- 1
  }else{
    diag(Delta2init[(1:z2),(1:z2)]) <- 1
  }
  
  ### Start of ADAM max for Y part 
  Lambdaold <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  Lambdanew <- as.matrix(rbind(cbind(Gamma1init%*%phi3in^(1/2), Delta1init%*%phi41in^(1/2), matrix(0, nrow = p1, ncol = z2)), cbind(Gamma2init%*%phi3in^(1/2),
                                                                                                                                    matrix(0, nrow = p2, ncol = z1), Delta2init%*%phi42in^(1/2))))
  k <- 1
  Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  alpha <- 0.001
  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 10^(-8)
  mold <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  uold <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  mnew <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  unew <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
  Phi2init <- diag(p1+p2)
  
  while((norm(Lambdanew - Lambdaold, type = "F") > 0.0001) & (k < 30000)){
    
    Lambdaold <- Lambdanew
    mold <- mnew
    uold <- unew
    
    gradMatrix <- 4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold)
    
    for(a in 1:j){
      for(b in 1:j){
        if(a >= b){
          grad.lam <- gradMatrix[a, b]
          mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
          unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
          Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        }
      }
    }
    for(a in (j+1):(p1+p2)){
      for(b in 1:j){
        grad.lam <- gradMatrix[a, b]
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      }
    }
    for(a in 1:z1){
      for(b in (j+1):(j+z1)){
        if(a >= (b - j)){
          grad.lam <- gradMatrix[a, b]
          mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
          unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
          Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        }
      }
    }
    for(a in (z1+1):p1){
      for(b in (j+1):(j+z1)){
        grad.lam <- gradMatrix[a, b]
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      }
    }
    for(a in (p1+1):(p1+z2)){
      for(b in (j+z1+1):(j+z1+z2)){
        if((a-p1) >= (b-j-z1)){
          grad.lam <- gradMatrix[a, b]
          mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
          unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
          Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        }
      }
    }
    for(a in (p1+z2+1):(p1+p2)){
      for(b in (j+z1+1):(j+z1+z2)){
        grad.lam <- gradMatrix[a, b]
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      }
    }
    k <- k+1
  }
  
  k2 <- k
  
  Smat <- matrix(0, nrow = p1+p2, ncol = p1+p2)
  mold <- rep(0, p1+p2)
  uold <- rep(0, p1+p2)
  mnew <- rep(0, p1+p2)
  unew <- rep(0, p1+p2)
  Phiold <- matrix(0, nrow = p1+p2, ncol = p1+p2)
  Phinew <-Phi2init
  grad.phi <- 0
  
  k <- 1
  
  while((norm(Phinew - Phiold, type = "F") > 0.0001) & (k < 30000)){
    
    Phiold <- Phinew
    mold <- mnew
    uold <- unew
    
    gradMatrix <- 2*(Lambdanew%*%t(Lambdanew) + Phiold - D)
    
    for(a in 1:(p1+p2)){
      grad.phi <- gradMatrix[a, a]
      mnew[a] <- beta1*mold[a] + (1-beta1)*grad.phi
      unew[a] <- max(beta2*uold[a], abs(grad.phi))
      Phinew[a,a] <- Phiold[a,a] - alpha/(1 - beta1^k)*mnew[a]/unew[a]
    }
    
    k <- k+1
    
  }
  
  Phinew[Phinew < 0] <- 1
  
  Gammaest <- Lambdanew[,1:j]
  if(j == 1){
    phi3fih <- Gammaest[1]
  }else{
    phi3fih <- diag(diag(Gammaest[(1:j),(1:j)]))
  }
  phi3fi <- phi3fih^2
  Gammaest <- Gammaest%*%solve(phi3fih)
  Gamma1est <- Gammaest[1:p1,]
  Gamma2est <- Gammaest[(p1+1):(p1+p2),]
  
  Delta1est <- Lambdanew[(1:p1),j+(1:z1)]
  if(z1 == 1){
    phi41fih <- Delta1est[1]
  }else{
    phi41fih <- diag(diag(Delta1est[(1:z1), (1:z1)]))
  }
  phi41fi <- phi41fih^2
  Delta1est <- Delta1est%*%solve(phi41fih)
  
  Delta2est <- Lambdanew[p1+(1:p2),j+z1+(1:z2)]
  if(z2 == 1){
    phi42fih <- Delta2est[1]
  }else{
    phi42fih <- diag(diag(Delta2est[(1:z2), (1:z2)]))
  }
  phi42fi <- phi42fih^2
  Delta2est <- Delta2est%*%solve(phi42fih)
  
  Phi21fi <- diag(diag(Phinew)[1:p1])
  Phi22fi <- diag(diag(Phinew)[p1+(1:p2)])
  
  ## Getting the loading matrix for the latent model connecting X and Y 
  Thetaest <- solve(t(Gamma1est)%*%Gamma1est)%*%t(Gamma1est)%*%G1
  Pi1est <- matrix(0, nrow = z1, ncol = m)
  Pi2est <- solve(t(Delta2est)%*%Delta2est)%*%t(Delta2est)%*%(G2 - Gamma2est%*%Thetaest)
  
  ## Getting the loading matrix for the latent model connecting X and Y 
  Psiest <- solve(t(Gamma1est)%*%Gamma1est)%*%t(Gamma1est)%*%G3
  Omega1est <- matrix(0, nrow = z1, ncol = u1 + u2)
  Omega2est <- solve(t(Delta2est)%*%Delta2est)%*%t(Delta2est)%*%(G4 - Gamma2est%*%Psiest)
  
  inits$Gamma1 <- as.matrix(Gamma1est)
  inits$Delta1 <- as.matrix(Delta1est)
  inits$Gamma2 <- as.matrix(Gamma2est)
  inits$Delta2 <- as.matrix(Delta2est)
  
  inits$Theta <- as.matrix(Thetaest)
  inits$Psi <- as.matrix(Psiest)
  inits$Pi1 <- as.matrix(Pi1est)
  inits$Omega1 <- as.matrix(Omega1est)
  inits$Pi2 <- as.matrix(Pi2est)
  inits$Omega2 <-as.matrix(Omega2est)
  
  inits$Phi21 <- as.matrix(Phi21fi)
  inits$Phi22 <- as.matrix(Phi22fi)
  
  inits$Phi3 <- as.matrix(phi3fi)
  inits$Phi41 <- as.matrix(phi41fi)
  inits$Phi42 <- as.matrix(phi42fi)
  
  return(inits)
  
}