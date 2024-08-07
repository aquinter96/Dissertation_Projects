
X_inits <- function(Data, m, u){

inits <- list()

X1 <- Data$X1
X2 <- Data$X2

u1 <- u[1]
u2 <- u[2]
q1 <- ncol(X1)
q2 <- ncol(X2)
n <- nrow(X1)

### To get initial values on A/B matrices for the X part ... 
jivemod <- jive(list(t(X1), t(X2)), m, c(u1, u2), method = "given")

jiveloads <- jive.predict(list(t(X1), t(X2)), jivemod)
Ajive <- jiveloads$joint.load
A3 <- Ajive%*%(eigen(1/n*crossprod(Ajive))$vectors)
A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))

## Initial value for covariance matrix for R (joint latent factor for X)
phirinh <- diag(diag(A5[(1:m),(1:m)]))
Ajiveest <- A5 %*% solve(phirinh)
phirin <- phirinh%*%phirinh
Ajiveest[(1:m),(1:m)][upper.tri(Ajiveest[(1:m),(1:m)], diag = F)] <- 0
diag(Ajiveest[(1:m),(1:m)]) <- 1
A1jiveest <- Ajiveest[1:q1,]
A2jiveest <- Ajiveest[(q1+1):(q1+q2),]

B1jive <- jiveloads$indiv.load[[1]]
B13 <- B1jive%*%(eigen(1/n*crossprod(B1jive))$vectors)
B15 <- B13%*%qr.Q(qr(t(B13[(1:u1),(1:u1)])))
phis1inh <- diag(diag(B15[(1:u1),(1:u1)]))
B1jiveest <- B15 %*% solve(phis1inh)
phis1in <- phis1inh%*%phis1inh
B1jiveest[(1:u1),(1:u1)][upper.tri(B1jiveest[(1:u1),(1:u1)], diag = F)] <- 0
diag(B1jiveest[(1:u1),(1:u1)]) <- 1

B2jive <- jiveloads$indiv.load[[2]]
B23 <- B2jive%*%(eigen(1/n*crossprod(B2jive))$vectors)
B25 <- B23%*%qr.Q(qr(t(B23[(1:u2),(1:u2)])))
phis2inh <- diag(diag(B25[(1:u2),(1:u2)]))
B2jiveest <- B25 %*% solve(phis2inh)
phis2in <- phis2inh%*%phis2inh
B2jiveest <- B25 %*% solve(diag(diag(B25[(1:u2),(1:u2)])))
B2jiveest[(1:u2),(1:u2)][upper.tri(B2jiveest[(1:u2),(1:u2)], diag = F)] <- 0
diag(B2jiveest[(1:u2),(1:u2)]) <- 1

## initial value for covariance of residuals, for X part  
phi11init <-diag(diag(cov(X1) - A1jiveest%*%t(A1jiveest) - B1jiveest%*%t(B1jiveest)))
phi12init <-diag(diag(cov(X2) - A2jiveest%*%t(A2jiveest) - B2jiveest%*%t(B2jiveest)))

phi11init[phi11init < 0] <- 1
phi12init[phi12init < 0] <- 1


## The ADAM algorithm part 
Lambdaold <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)

Lambdanew <- as.matrix(rbind(cbind(A1jiveest%*%phirin^(1/2), B1jiveest%*%phis1in^(1/2), matrix(0, nrow = q1, ncol = u2)), cbind(A2jiveest%*%phirin^(1/2),
                                                                                                                                matrix(0, nrow = q2, ncol = u1), B2jiveest%*%phis2in^(1/2))))

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
Phi1init <- diag(c(diag(phi11init), diag(phi12init)))
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

k11 <- k


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

k12 <- k

Aest <- Lambdanew[,1:m]
phirfih <- diag(diag(Aest[(1:m), (1:m)]))
phirfi <- phirfih^2
Aest <- Aest%*%solve(phirfih)
A1est <- Aest[1:q1,]
A2est <- Aest[(q1+1):(q1+q2),]

B1est <- Lambdanew[(1:q1),m+(1:u1)]
phis1fih <- diag(diag(B1est[(1:u1), (1:u1)]))
phis1fi <- phis1fih^2
B1est <- B1est%*%solve(phis1fih)

B2est <- Lambdanew[q1+(1:q2),m+u1+(1:u2)]

phis2fih <- diag(diag(B2est[(1:u2), (1:u2)]))
phis2fi <- phis2fih^2
B2est <- B2est%*%solve(phis2fih)

Phi11est <- diag(diag(Phinew)[1:q1])
Phi12est <- diag(diag(Phinew)[q1+(1:q2)])

inits$A1 <- A1est
inits$B1 <- B1est
inits$A2 <- A2est
inits$B2 <- B2est

inits$PhiR <- phirfi
inits$PhiS1 <- phis1fi
inits$PhiS2 <- phis2fi
inits$Phi11 <- Phi11est
inits$Phi12 <- Phi12est

return(inits)

}
## EM algorithm for the X part parameters 

