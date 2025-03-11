## Initial_ValsX.R

Random_ValsX <- list()
Perturb_ValsX <- list()
Jive_ValsX <- list()
AdaMax_ValsX <- list()

X1 <- MyData$X1
X2 <- MyData$X2
q1 <- ncol(X1)
q2 <- ncol(X2)
n <- nrow(X1)
m <- ncol(Model_Params$A1)
u1 <- ncol(Model_Params$B1)
u2 <- ncol(Model_Params$B2)

Random_ValsX$A1 <- Model_Params$A1
Random_ValsX$B1 <- Model_Params$B1
Random_ValsX$A2 <- Model_Params$A2
Random_ValsX$B2 <- Model_Params$B2
Random_ValsX$PhiR <- diag(runif(m, -1, 1))
Random_ValsX$PhiS1 <- diag(runif(u1, -1, 1))
Random_ValsX$PhiS2 <- diag(runif(u2, -1, 1))
Random_ValsX$Phi11 <- diag(runif(q1, -1, 1))
Random_ValsX$Phi12 <- diag(runif(q2, -1, 1))

Perturb_ValsX$A1 <- Model_Params$A1
Perturb_ValsX$B1 <- Model_Params$B1
Perturb_ValsX$A2 <- Model_Params$A2
Perturb_ValsX$B2 <- Model_Params$B2
Perturb_ValsX$PhiR <- diag(diag(Model_Params$PhiR + rnorm(m, 0, 0.01)))
Perturb_ValsX$PhiS1 <- diag(diag(Model_Params$PhiS1 + rnorm(u1, 0, 0.01)))
Perturb_ValsX$PhiS2 <- diag(diag(Model_Params$PhiS2 + rnorm(u2, 0, 0.01)))
Perturb_ValsX$Phi11 <- diag(diag(Model_Params$Phi11 + rnorm(q1, 0, 0.01)))
Perturb_ValsX$Phi12 <- diag(diag(Model_Params$Phi12 + rnorm(q2, 0, 0.01)))

for(i in 1:q1){
  for(j in 1:m){
    if(i > j){
      Random_ValsX$A1[i,j] <- runif(1, -1, 1)
      Perturb_ValsX$A1[i,j] <- Perturb_ValsX$A1[i,j] + rnorm(1, 0, 0.01)
    }
  }
  for(j in 1:u1){
    if(i > j){
      Random_ValsX$B1[i,j] <- runif(1, -1, 1)
      Perturb_ValsX$B1[i,j] <- Perturb_ValsX$B1[i,j] + rnorm(1, 0, 0.01)
    }
  }
}
for(i in 1:q2){
  for(j in 1:m){
    if(i > j){
      Random_ValsX$A2[i,j] <- runif(1, -1, 1)
      Perturb_ValsX$A2[i,j] <- Perturb_ValsX$A2[i,j] + rnorm(1, 0, 0.01)
    }
  }
  for(j in 1:u2){
    if(i > j){
      Random_ValsX$B2[i,j] <- runif(1, -1, 1)
      Perturb_ValsX$B2[i,j] <- Perturb_ValsX$B2[i,j] + rnorm(1, 0, 0.01)
    }
  }
}




jivemod <- jive(list(t(X1), t(X2)), m, c(u1, u2), method = "given")
jiveloads <- jive.predict(list(t(X1), t(X2)), jivemod)
Ajive <- jiveloads$joint.load
A3 <- Ajive%*%(eigen(1/nrow(X1)*crossprod(Ajive))$vectors)
A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))

## Initial value for covariance matrix for R (joint latent factor for X)
phirinh <- diag(diag(A5[(1:m),(1:m)]))
Ajive <- A5 %*% solve(phirinh)
PhiRjive <- phirinh%*%phirinh
Ajive[(1:m),(1:m)][upper.tri(Ajive[(1:m),(1:m)], diag = F)] <- 0
diag(Ajive[(1:m),(1:m)]) <- 1
A1jive <- Ajive[1:q1,]
A2jive <- Ajive[(q1+1):(q1+q2),]
Jive_ValsX$A1 <- A1jive
Jive_ValsX$A2 <- A2jive
Jive_ValsX$PhiR <- PhiRjive

B1jive <- jiveloads$indiv.load[[1]]
B13 <- B1jive%*%(eigen(1/n*crossprod(B1jive))$vectors)
B15 <- B13%*%qr.Q(qr(t(B13[(1:u1),(1:u1)])))
phis1inh <- diag(diag(B15[(1:u1),(1:u1)]))
B1jive <- B15 %*% solve(phis1inh)
PhiS1jive <- phis1inh%*%phis1inh
B1jive[(1:u1),(1:u1)][upper.tri(B1jive[(1:u1),(1:u1)], diag = F)] <- 0
diag(B1jive[(1:u1),(1:u1)]) <- 1
Jive_ValsX$B1 <- B1jive
Jive_ValsX$PhiS1 <- PhiS1jive

B2jive <- jiveloads$indiv.load[[2]]
B23 <- B2jive%*%(eigen(1/n*crossprod(B2jive))$vectors)
B25 <- B23%*%qr.Q(qr(t(B23[(1:u2),(1:u2)])))
phis2inh <- diag(diag(B25[(1:u2),(1:u2)]))
B2jiveest <- B25 %*% solve(phis2inh)
PhiS2jive <- phis2inh%*%phis2inh
B2jive <- B25 %*% solve(diag(diag(B25[(1:u2),(1:u2)])))
B2jive[(1:u2),(1:u2)][upper.tri(B2jive[(1:u2),(1:u2)], diag = F)] <- 0
diag(B2jive[(1:u2),(1:u2)]) <- 1
Jive_ValsX$B2 <- B2jive
Jive_ValsX$PhiS2 <- PhiS2jive

## initial value for covariance of residuals, for X part  
Phi11jive <-diag(diag(cov(X1) - A1jive%*%t(A1jive) - B1jive%*%t(B1jive)))
Phi12jive <-diag(diag(cov(X2) - A2jive%*%t(A2jive) - B2jive%*%t(B2jive)))

Phi11jive[Phi11jive < 0] <- 1
Phi12jive[Phi12jive < 0] <- 1
Jive_ValsX$Phi11 <- Phi11jive
Jive_ValsX$Phi12 <- Phi12jive

Jive_ValsX <- Jive_ValsX[c(1, 4, 2, 6, 3, 5, 7:9)]

Lambdaold <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)

Lambdanew <- as.matrix(rbind(cbind(A1jive%*%PhiRjive^(1/2), B1jive%*%PhiS1jive^(1/2), matrix(0, nrow = q1, ncol = u2)), cbind(A2jive%*%PhiRjive^(1/2),
             matrix(0, nrow = q2, ncol = u1), B2jive%*%PhiS2jive^(1/2))))

k <- 1 ## a counter 
Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
alpha <- 0.001
beta1 <- 0.9
beta2 <- 0.999
epsilon <- 10^(-8)
mold <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
uold <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
mnew <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
unew <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
Phi1init <- diag(c(diag(Phi11jive), diag(Phi12jive)))
D <- cov(cbind(X1, X2))

## Adam for the X part: the loading matrices 
while((norm(Lambdanew - Lambdaold, type = "F") > 0.0001) & (k < 20000)){
  
  Lambdaold <- Lambdanew
  mold <- mnew
  uold <- unew
  
  for(a in 1:m){
    for(b in 1:m){
      if(a >= b){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
      }
    }
  }
  for(a in (m+1):(q1+q2)){
    for(b in 1:m){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
    }
  }
  
  for(a in 1:u1){
    for(b in (m+1):(m+u1)){
      if(a >= (b - m)){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
      }
    }
  }
  for(a in (u1+1):q1){
    for(b in (m+1):(m+u1)){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
    }
  }
  
  for(a in (q1+1):(q1+u2)){
    for(b in (m+u1+1):(m+u1+u2)){
      if((a-q1) >= (b-m-u1)){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
      }
    }
  }
  for(a in (q1+u2+1):(q1+q2)){
    for(b in (m+u1+1):(m+u1+u2)){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi1init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)
    }
  }
  k <- k+1
}

### ADAM max for error term covariance the X part 
Smat <- matrix(0, nrow = q1+q2, ncol = q1+q2)
mold <- rep(0, q1+q2)
uold <- rep(0, q1+q2)
mnew <- rep(0, q1+q2)
unew <- rep(0, q1+q2)
Phiold <- matrix(0, nrow = q1+q2, ncol = q1+q2)
Phinew <-Phi1init
grad.phi <- 0

k <- 1

while((norm(Phinew - Phiold, type = "F") > 0.0001) & (k < 10000)){
  
  Phiold <- Phinew
  mold <- mnew
  uold <- unew
  
  for(a in 2:(q1+q2)){
    Smat[a,a] <- 1
    grad.phi <- sum(diag(2*(Lambdanew%*%t(Lambdanew) + Phiold - D)%*%Smat))
    mnew[a] <- beta1*mold[a] + (1-beta1)*grad.phi
    unew[a] <- max(beta2*uold[a], abs(grad.phi))
    Phinew[a,a] <- Phiold[a,a] - alpha/(1 - beta1^k)*mnew[a]/unew[a]
    Smat <- matrix(0, nrow = q1+q2, ncol = q1+q2)
  }
  
  k <- k+1
}

Phinew[Phinew < 0] <- 1

Aadam <- Lambdanew[,1:m]
phirfih <- diag(diag(Aadam[(1:m), (1:m)]))
PhiRadam <- phirfih^2
Aadam <- Aadam%*%phirfih
A1adam <- Aadam[1:q1,]
A2adam <- Aadam[(q1+1):(q1+q2),]
AdaMax_ValsX$A1 <- A1adam
AdaMax_ValsX$A2 <- A2adam
AdaMax_ValsX$PhiR <- PhiRadam

B1adam <- Lambdanew[(1:q1),m+(1:u1)]
phis1fih <- diag(diag(B1adam[(1:u1), (1:u1)]))
PhiS1adam <- phis1fih^2
B1adam <- B1adam%*%phis1fih
AdaMax_ValsX$B1 <- B1adam
AdaMax_ValsX$PhiS1 <- PhiS1adam

B2adam <- Lambdanew[q1+(1:q2),m+u1+(1:u2)]
phis2fih <- diag(diag(B2adam[(1:u2), (1:u2)]))
PhiS2adam <- phis2fih^2
B2adam <- B2adam%*%phis2fih
AdaMax_ValsX$B2 <- B2adam
AdaMax_ValsX$PhiS2 <- PhiS2adam

AdaMax_ValsX$Phi11 <- diag(diag(Phinew)[1:q1])
AdaMax_ValsX$Phi12 <- diag(diag(Phinew)[q1+(1:q2)])

AdaMax_ValsX <- AdaMax_ValsX[c(1, 4, 2, 6, 3, 5, 7:9)]

