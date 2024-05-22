library(MASS)
library(pracma)
library(psych)
library(r.jive)

args <- commandArgs(TRUE)

datset <- as.numeric(args[1])

set.seed(12543)

#initialize parameters
j = 5
m = 2
q1 = 20
u1 = 2
q2 = 30
u2 = 3
p1 = 30
z1 = 3
p2 = 20
z2 = 3
u = u1 + u2
A1 = matrix(rnorm(q1*m),nrow = q1)
A1[(1:m),(1:m)][upper.tri(A1[(1:m),(1:m)], diag = F)] <- 0
diag(A1[(1:m),(1:m)]) <- 1
B1 = matrix(rnorm(q1*u1),nrow = q1)
B1[(1:u1),(1:u1)][upper.tri(B1[(1:u1),(1:u1)], diag = F)] <- 0
diag(B1[(1:u1),(1:u1)]) <- 1
A2 = matrix(rnorm(q2*m),nrow = q2)
B2 = matrix(rnorm(q2*u2),nrow = q2)
B2[(1:u2),(1:u2)][upper.tri(B2[(1:u2),(1:u2)], diag = F)] <- 0
diag(B2[(1:u2),(1:u2)]) <- 1
Gamma1 = matrix(rnorm(p1*j),nrow = p1)
Gamma1[(1:j),(1:j)][upper.tri(Gamma1[(1:j),(1:j)], diag = F)] <- 0
diag(Gamma1[(1:j),(1:j)]) <- 1
Delta1 = nullspace(ginv(Gamma1))%*%matrix(rnorm(ncol(nullspace(ginv(Gamma1)))*z1),nrow=ncol(nullspace(ginv(Gamma1))))
Delta15 <- Delta1%*%qr.Q(qr(t(Delta1[(1:z1),(1:z1)])))
Delta1 <- Delta15 %*% solve(diag(diag(Delta15[(1:z1),(1:z1)])))
Delta1[(1:z1),(1:z1)][upper.tri(Delta1[(1:z1),(1:z1)], diag = F)] <- 0
diag(Delta1[(1:z1),(1:z1)]) <- 1
Gamma2 = matrix(rnorm(p2*j),nrow = p2)
Delta2 = nullspace(ginv(Gamma2))%*%matrix(rnorm(ncol(nullspace(ginv(Gamma2)))*z2),nrow=ncol(nullspace(ginv(Gamma2))))
Delta25 <- Delta2%*%qr.Q(qr(t(Delta2[(1:z2),(1:z2)])))
Delta2 <- Delta25 %*% solve(diag(diag(Delta25[(1:z2),(1:z2)])))
Delta2[(1:z2),(1:z2)][upper.tri(Delta2[(1:z2),(1:z2)], diag = F)] <- 0
diag(Delta2[(1:z2),(1:z2)]) <- 1
Theta <- matrix(rnorm(j*m),nrow = j)
Psi <- matrix(rnorm(j*u),nrow = j)
Pi1 <- matrix(rnorm(z1*m),nrow = z1)
Pi2 <- matrix(rnorm(z2*m),nrow = z2)
Omega1 <- matrix(rnorm(z1*u),nrow = z1)
Omega2 <- matrix(rnorm(z2*u),nrow = z2)

PhiR <- diag(abs(rnorm(m)))
PhiS1 <- diag(abs(rnorm(u1)))
PhiS2 <- diag(abs(rnorm(u2)))
Phi11 <- diag(abs(rnorm(q1)))
Phi12 <- diag(abs(rnorm(q2)))
Phi21 <- diag(abs(rnorm(p1)))
Phi22 <- diag(abs(rnorm(p2)))
Phi3 <- diag(abs(rnorm(j)))
Phi41 <- diag(abs(rnorm(z1)))
Phi42 <- diag(abs(rnorm(z2)))

PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))


MSEa1in <- NA
MSEa2in <- NA
MSEa1fi <- NA
MSEa2fi <- NA
MSEb1in <- NA
MSEb2in <- NA
MSEb1fi <- NA
MSEb2fi <- NA

MSEd1in <- NA
MSEd2in <- NA
MSEd1fi <- NA
MSEd2fi <- NA
MSEg1in <- NA
MSEg2in <- NA
MSEg1fi <- NA
MSEg2fi <- NA

set.seed(datset)
n = 1000

## to simulte latent factor and then X/Y
R <- t(MASS::mvrnorm(n, rep(0, m), PhiR))
S1 <- t(MASS::mvrnorm(n, rep(0, u1), PhiS1))
S2 <- t(MASS::mvrnorm(n, rep(0, u2), PhiS2))

X1 <- A1%*%R + B1%*%S1 + t(MASS::mvrnorm(n, rep(0, q1), Phi11))
X2 <- A2%*%R + B2%*%S2 + t(MASS::mvrnorm(n, rep(0, q2), Phi12))

V <- Theta%*%R + Psi%*%rbind(S1,S2) + t(MASS::mvrnorm(n, rep(0, j), Phi3))
W1 <- Pi1%*%R + Omega1%*%rbind(S1,S2) + t(MASS::mvrnorm(n, rep(0, z1), Phi41))
W2 <- Pi2%*%R + Omega2%*%rbind(S1,S2) + t(MASS::mvrnorm(n, rep(0, z2), Phi42))

Y1 <- Gamma1%*%V + Delta1%*%W1 + t(MASS::mvrnorm(n, rep(0, p1), Phi21))
Y2 <- Gamma2%*%V + Delta2%*%W2 + t(MASS::mvrnorm(n, rep(0, p2), Phi22))

#True values of G1-G4, for intermediate checking ... 
G1t <- Gamma1%*%Theta + Delta1%*%Pi1
G2t <- Gamma2%*%Theta + Delta2%*%Pi2
G3t <- Gamma1%*%Psi + Delta1%*%Omega1
G4t <- Gamma2%*%Psi + Delta2%*%Omega2
Gt <- rbind(cbind(G1t, G3t), cbind(G2t, G4t))

D <- cov(cbind(t(X1), t(X2)))


### To get initial values on A/B matrices for the X part ... 
# jivemod <- jive(list(X1, X2), m, c(2, 3), method = "given")
# jiveloads <- jive.predict(list(X1, X2), jivemod)
# Ajive <- jiveloads$joint.load
# A3 <- Ajive%*%(eigen(1/n*crossprod(Ajive))$vectors)
# A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))
# 
# ## Initial value for covariance matrix for R (joint latent factor for X)
# phirinh <- diag(diag(A5[(1:m),(1:m)]))
# Ajiveest <- A5 %*% solve(phirinh)
# phirin <- phirinh%*%phirinh
# Ajiveest[(1:m),(1:m)][upper.tri(Ajiveest[(1:m),(1:m)], diag = F)] <- 0
# diag(Ajiveest[(1:m),(1:m)]) <- 1
# A1jiveest <- Ajiveest[1:q1,]
# A2jiveest <- Ajiveest[(q1+1):(q1+q2),]
# 
# MSEa1in <- 1/((q1*m)-m*(m-1)/2)*sum((A1jiveest - A1)^2)
# MSEa2in <- 1/(q2*m)*sum((A2jiveest - A2)^2)
# 
# B1jive <- jiveloads$indiv.load[[1]]
# B13 <- B1jive%*%(eigen(1/n*crossprod(B1jive))$vectors)
# B15 <- B13%*%qr.Q(qr(t(B13[(1:u1),(1:u1)])))
# phis1inh <- diag(diag(B15[(1:u1),(1:u1)]))
# B1jiveest <- B15 %*% solve(phis1inh)
# phis1in <- phis1inh%*%phis1inh
# B1jiveest[(1:u1),(1:u1)][upper.tri(B1jiveest[(1:u1),(1:u1)], diag = F)] <- 0
# diag(B1jiveest[(1:u1),(1:u1)]) <- 1
# 
# MSEb1in <- 1/((q1*u1)-u1*(u1-1)/2)*sum((B1jiveest - B1)^2)
# 
# B2jive <- jiveloads$indiv.load[[2]]
# B23 <- B2jive%*%(eigen(1/n*crossprod(B2jive))$vectors)
# B25 <- B23%*%qr.Q(qr(t(B23[(1:u2),(1:u2)])))
# phis2inh <- diag(diag(B25[(1:u2),(1:u2)]))
# B2jiveest <- B25 %*% solve(phis2inh)
# phis2in <- phis2inh%*%phis2inh
# B2jiveest <- B25 %*% solve(diag(diag(B25[(1:u2),(1:u2)])))
# B2jiveest[(1:u2),(1:u2)][upper.tri(B2jiveest[(1:u2),(1:u2)], diag = F)] <- 0
# diag(B2jiveest[(1:u2),(1:u2)]) <- 1
# 
# MSEb2in <- 1/((q2*u2)-u2*(u2-1)/2)*sum((B2jiveest - B2)^2)
# 
# MSEphirinit <-1/m*sum((phirin - PhiR)^2)
# MSEphis1init <- 1/u1*sum((phis1in - PhiS1)^2)
# MSEphis2init <- 1/u2*sum((phis2in - PhiS2)^2)
# 
# 
# ## initial value for covariance of residuals, for X part  
# phi11init <-diag(diag(cov(t(X1)) - A1jiveest%*%t(A1jiveest) - B1jiveest%*%t(B1jiveest)))
# phi12init <-diag(diag(cov(t(X2)) - A2jiveest%*%t(A2jiveest) - B2jiveest%*%t(B2jiveest)))
# 
# phi11init[phi11init < 0] <- 1
# phi12init[phi12init < 0] <- 1
# 
# 
# MSEphi11in <- 1/q1*sum((phi11init - Phi11)^2)
# MSEphi12in <- 1/q2*sum((phi12init - Phi12)^2)


## The ADAM algorithm part 
Lambdaold <- matrix(0, nrow = q1+q2, ncol = m+u1+u2)

Lambdanew <- as.matrix(rbind(cbind(A1%*%PhiR^(1/2), B1%*%PhiS1^(1/2), matrix(0, nrow = q1, ncol = u2)), cbind(A2%*%PhiR^(1/2),
                                             matrix(0, nrow = q2, ncol = u1), B2%*%PhiS2^(1/2))))

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
Phi1init <- diag(c(diag(Phi11), diag(Phi12)))

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
Aest <- Aest%*%phirfih
A1est <- Aest[1:q1,]
A2est <- Aest[(q1+1):(q1+q2),]

MSEa1fi <- 1/((q1*m)-m*(m-1)/2)*sum((A1est - A1)^2)
MSEa2fi <- 1/(q2*m)*sum((A2est - A2)^2)

B1est <- Lambdanew[(1:q1),m+(1:u1)]
phis1fih <- diag(diag(B1est[(1:u1), (1:u1)]))
phis1fi <- phis1fih^2
B1est <- B1est%*%phis1fih

MSEb1fi <- 1/((q1*u1)-u1*(u1-1)/2)*sum((B1est - B1)^2)

B2est <- Lambdanew[q1+(1:q2),m+u1+(1:u2)]

phis2fih <- diag(diag(B2est[(1:u2), (1:u2)]))
phis2fi <- phis2fih^2
B2est <- B2est%*%phis2fih

MSEb2fi <- 1/((q2*u2)-u2*(u2-1)/2)*sum((B2est - B2)^2)

Phi11fi <- diag(diag(Phinew)[1:q1])
Phi12fi <- diag(diag(Phinew)[q1+(1:q2)])

MSEphi11fi <- 1/q1*sum((Phi11fi - Phi11)^2)
MSEphi12fi <- 1/q2*sum((Phi12fi - Phi12)^2)

MSEphirfi <- 1/m*sum((phirfi - PhiR)^2)
MSEphis1fi <- 1/u1*sum((phis1fi - PhiS1)^2)
MSEphis2fi <- 1/u2*sum((phis2fi - PhiS2)^2)

## EM algorithm for the X part parameters 
EMAlgXAdLassoCV <- function(Xpanel = list(), x_joint, x_indiv = c(), Ainits = list(), Binits = list(), vcovinits = list(), tuningpA = seq(5.1,5.5,0.1), tuningpB = seq(5.1,5.5,0.1), weightsA = list(), weightsB = list()){
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B
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
  
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  
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
    
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
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
            
            while((norm(A1update - A1CV, type = "F") > 0.005) | (norm(A2update - A2CV, type = "F") > 0.005) | (norm(B1update - B1CV) > 0.005) | (norm(B2update - B2CV) > 0.005)){
              
              A1update <- A1CV
              A2update <- A2CV
              B1update <- B1CV
              B2update <- B2CV
              
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
              #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
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
              #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
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
          
          #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
          log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]],Xpanel[[2]]))+(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]],Xpanel[[2]])) - colMeans(cbind(Xpanel[[1]],Xpanel[[2]])))))) + (q1+q2)*log((2*pi))+log(det(sigma11)))
          
          BICmat[a,b] <- log(n)*(sum(A1CV != 0) + sum(A2CV != 0) + sum(B1CV != 0) + sum(B2CV != 0) - (m+u1+u2)) - 2*log.lik
        }
      }
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

test <- EMAlgXAdLassoCV(Xpanel = list(t(X1), t(X2)), x_joint = m, x_indiv = c(u1, u2), Ainits = list(A1est, A2est), Binits = list(B1est, B2est), vcovinits = list(phirfi, phis1fi, phis2fi, Phi11fi, Phi12fi), tuningpA = 0, tuningpB = 0, weightsA = list(A1est, A2est), weightsB = list(B1est, B2est))

A1em <- test$A1
A2em <- test$A2

MSEa1em <- 1/((q1*m)-m*(m-1)/2)*sum((A1em - A1)^2)
MSEa2em <- 1/(q2*m)*sum((A2em - A2)^2)

B1em <- test$B1

MSEb1em <- 1/((q1*u1)-u1*(u1-1)/2)*sum((B1em - B1)^2)

B2em <- test$B2

MSEb2em <- 1/((q2*u2)-u2*(u2-1)/2)*sum((B2em - B2)^2)

PhiRem <- test$PhiR
PhiS1em <- test$PhiS1
PhiS2em <- test$PhiS2
Phi11em <- test$Phi11
Phi12em <- test$Phi12

MSEphirem <- 1/m*sum((PhiRem - PhiR)^2)
MSEphis1em <- 1/m*sum((PhiS1em - PhiS1)^2)
MSEphis2em <- 1/m*sum((PhiS2em - PhiS2)^2)
MSEphi11em <- 1/q1*sum((Phi11em - Phi11)^2)
MSEphi12em <- 1/q1*sum((Phi12em - Phi12)^2)

PhiS1Sem <- cbind(PhiS1em, matrix(0, nrow = u1, ncol = u2))
PhiS2Sem <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2em)

##### The Y part 
## to prepare initials (with estimates from X part fixed) for Y part parameters 
##   To serve as initial values for the ADAM Max 
# C <- rbind(cbind(A1em%*%test$PhiR, B1em%*%PhiS1Sem), cbind(A2em, B2em%*%PhiS2Sem))
# 
# Ga <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(t(X1),t(Y1)), cov(t(X2),t(Y1)))
# Gb <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(t(X1),t(Y2)), cov(t(X2),t(Y2)))
# G1 <- t(Ga[1:m,])
# G2 <- t(Gb[1:m,])
# G3 <- t(Ga[(m+1):(m+u),])
# G4 <- t(Gb[(m+1):(m+u),])
# Ycov <- cov(cbind(t(Y1), t(Y2)))
# D <- Ycov - rbind(t(Ga), t(Gb))%*%t(rbind(t(Ga), t(Gb)))
# 
# G1G2svd <- svd(D[(1:p1), (p1+(1:p2))])
# G1init <- G1G2svd$u[,(1:j)]
# G2init <- G1G2svd$v[,(1:j)]
# Gammainit <- rbind(G1init, G2init)
# Gamma3 <- Gammainit%*%(eigen(1/n*crossprod(Gammainit))$vectors)
# Gamma5 <- Gamma3%*%qr.Q(qr(t(Gamma3[(1:j),(1:j)])))
# phi3inh <- diag(diag(Gamma5[(1:j),(1:j)]))
# Gammainit <- Gamma5 %*% solve(phi3inh)
# phi3in <- phi3inh%*%phi3inh
# Gammainit[(1:j),(1:j)][upper.tri(Gammainit[(1:j),(1:j)], diag = F)] <- 0
# diag(Gammainit[(1:j),(1:j)]) <- 1
# Gamma1init <- Gammainit[(1:p1),]
# Gamma2init <- Gammainit[(p1+(1:p2)),]
# 
# MSEg1in <- 1/((p1*j)-j*(j-1)/2)*sum((Gamma1init - Gamma1)^2)
# MSEg2in <- 1/(p2*j)*sum((Gamma2init - Gamma2)^2)
# MSEphi3in <- 1/j*sum((phi3in - Phi3)^2)
# 
# # Delta1obj <-function(a){
# #   amat <- matrix(0, nrow = p1, ncol = z1)
# #   diag(amat[(1:z1),(1:z1)]) <- 1
# #   amat[(1:z1),(1:z1)][lower.tri(amat[(1:z1),(1:(z1))], diag = F)] <- a[1:(z1*(z1-1)/2)]
# #   amat[(z1+1):p1,] <- a[(z1*(z1-1)/2+1):(p1*z1-z1*(z1+1)/2)]
# #   norm((D[(1:p1), (1:p1)] - Gamma1init%*%phi3in%*%t(Gamma1init) - amat%*%t(amat) - diag(p1)), type = "F")
# # }
# 
# Delta1init <- svd(D[(1:p1), (1:p1)] - Gamma1init%*%phi3in%*%t(Gamma1init) - diag(p1))$u[,1:z1]
# Delta13 <- Delta1init%*%(eigen(1/n*crossprod(Delta1init))$vectors)
# Delta15 <- Delta13%*%qr.Q(qr(t(Delta13[(1:z1),(1:z1)])))
# phi41inh <- diag(diag(Delta15[(1:z1),(1:z1)]))
# Delta1init <- Delta15 %*% solve(phi41inh)
# phi41in <- phi41inh%*%phi41inh
# Delta1init[(1:z1),(1:z1)][upper.tri(Delta1init[(1:z1),(1:z1)], diag = F)] <- 0
# diag(Delta1init[(1:z1),(1:z1)]) <- 1
# 
# # Delta1initvec <- c(Delta1init[(1:z1),(1:z1)][lower.tri(Delta1init[(1:z1),(1:z1)], diag = F)], Delta1init[(z1+1):p1,])
# # Delta1newvec <- nlminb(Delta1initvec, Delta1obj)
# # Delta1init <- matrix(0, nrow = p1, ncol = z1)
# # diag(Delta1init[(1:z1),(1:z1)]) <- 1
# # Delta1init[(1:z1),(1:z1)][lower.tri(Delta1init[(1:z1),(1:(z1))], diag = F)] <- Delta1newvec$par[1:(z1*(z1-1)/2)]
# # Delta1init[(z1+1):p1,] <- Delta1newvec$par[(z1*(z1-1)/2+1):(p1*z1-z1*(z1+1)/2)]
# 
# MSEd1in <- 1/((p1*z1)-z1*(z1-1)/2)*sum((Delta1init - Delta1)^2)
# MSEphi41in <- 1/z1*sum((phi41in - Phi41)^2)
# 
# Phi21obj <-function(a){
#   amat <- matrix(0, nrow = p1, ncol = p1)
#   diag(amat) <- a
#   norm((D[(1:p1), (1:p1)] - Gamma1init%*%phi3in%*%t(Gamma1init) - Delta1init%*%t(Delta1init) - amat), type = "F")
# }
# 
# Phi21init <- diag(diag(D[(1:p1), (1:p1)] - Gamma1init%*%t(Gamma1init) - Delta1init%*%t(Delta1init)))
# Phi21init <-diag(p1)
# Phi21initvec <- rep(1,p1)
# Phi21newvec <- nlminb(Phi21initvec, Phi21obj)
# diag(Phi21init) <- Phi21newvec$par
# 
# MSEphi21in <- 1/p1*sum((Phi21init - Phi21)^2)
# 
# Delta2obj <-function(a){
#   amat <- matrix(0, nrow = p2, ncol = z2)
#   diag(amat[(1:z2),(1:(z2))]) <- 1
#   amat[(1:z2),(1:(z2))][lower.tri(amat[(1:z2),(1:(z2))], diag = F)] <- a[1:z2]
#   amat[(z2+1):p2,] <- a[(z2+1):(p2*z2-z2*(z2+1)/2)]
#   norm((D[(p1+(1:p2)), (p1+(1:p2))] - Gamma2init%*%t(Gamma2init) - amat%*%t(amat) - diag(p2)), type = "F")
# }
# 
# Delta2init <- svd(D[(p1+(1:p2)), (p1+(1:p2))] - Gamma2init%*%phi3in%*%t(Gamma2init) - diag(p2))$u[,1:z2]
# Delta23 <- Delta2init%*%(eigen(1/n*crossprod(Delta2init))$vectors)
# Delta25 <- Delta23%*%qr.Q(qr(t(Delta23[(1:z2),(1:z2)])))
# phi42inh <- diag(diag(Delta25[(1:z2),(1:z2)]))
# Delta2init <- Delta25 %*% solve(phi42inh)
# phi42in <- phi42inh%*%phi42inh
# Delta2init[(1:z2),(1:z2)][upper.tri(Delta2init[(1:z2),(1:z2)], diag = F)] <- 0
# diag(Delta2init[(1:z2),(1:z2)]) <- 1

# Delta2init <- matrix(0, nrow = p2, ncol = z2)
# Delta2initvec <- c(Delta2init[(1:z2),(1:(z2))][lower.tri(Delta2init[(1:z2),(1:(z2))], diag = F)], Delta2init[(z2+1):p2,])
# Delta2newvec <- nlminb(Delta2initvec, Delta2obj)
# Delta2init <- matrix(0, nrow = p2, ncol = z2)
# diag(Delta2init[(1:z2),(1:(z2))]) <- 1
# Delta2init[(1:z2),(1:(z2))][lower.tri(Delta2init[(1:z2),(1:(z2))], diag = F)] <- Delta2newvec$par[1:z2]
# Delta2init[(z2+1):p2,] <- Delta2newvec$par[(z2+1):(p2*z2-z2*(z2+1)/2)]

# MSEd2in <- 1/((p2*z2)-z2*(z2-1)/2)*sum((Delta2init - Delta2)^2)
# MSEphi42in <- 1/z2*sum((phi42in - Phi42)^2)
# 
# Phi22obj <-function(a){
#   amat <- matrix(0, nrow = p2, ncol = p2)
#   diag(amat) <- a
#   norm((D[(p1+(1:p2)), (p1+(1:p2))] - Gamma2init%*%t(Gamma2init) - Delta2init%*%t(Delta2init) - amat), type = "F")
# }
# 
# Phi22init <-diag(p2)
# Phi22initvec <- rep(1,p2)
# Phi22newvec <- nlminb(Phi22initvec, Phi22obj)
# diag(Phi22init) <- Phi22newvec$par

### Start of ADAM max for Y part 
Lambdaold <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
Lambdanew <- as.matrix(rbind(cbind(Gamma1%*%Phi3^(1/2), Delta1%*%Phi41^(1/2), matrix(0, nrow = p1, ncol = z2)), cbind(Gamma2%*%Phi3^(1/2),
                                              matrix(0, nrow = p2, ncol = z1), Delta2%*%Phi42^(1/2))))


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
Phi2init <- diag(c(diag(Phi21), diag(Phi22)))

while((norm(Lambdanew - Lambdaold, type = "F") > 0.0001) & (k < 30000)){
  
  Lambdaold <- Lambdanew
  mold <- mnew
  uold <- unew
  
  for(a in 1:j){
    for(b in 1:j){
      if(a >= b){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
      }
    }
  }
  for(a in (j+1):(p1+p2)){
    for(b in 1:j){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
    }
  }
  for(a in 1:z1){
    for(b in (j+1):(j+z1)){
      if(a >= (b - j)){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
      }
    }
  }
  for(a in (z1+1):p1){
    for(b in (j+1):(j+z1)){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
    }
  }
  for(a in (p1+1):(p1+z2)){
    for(b in (j+z1+1):(j+z1+z2)){
      if((a-p1) >= (b-j-z1)){
        Smat[a,b] <- 1
        grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
        mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
        unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
        Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
        Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
      }
    }
  }
  for(a in (p1+z2+1):(p1+p2)){
    for(b in (j+z1+1):(j+z1+z2)){
      Smat[a,b] <- 1
      grad.lam <- sum(diag(t(4*(Lambdaold%*%t(Lambdaold)%*%Lambdaold - D%*%Lambdaold + Phi2init%*%Lambdaold))%*%Smat))
      mnew[a,b] <- beta1*mold[a,b] + (1-beta1)*grad.lam
      unew[a,b] <- max(beta2*uold[a,b], abs(grad.lam))
      Lambdanew[a,b] <- Lambdaold[a,b] - alpha/(1 - beta1^k)*mnew[a,b]/unew[a,b]
      Smat <- matrix(0, nrow = p1+p2, ncol = j+z1+z2)
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
  
  for(a in 2:(p1+p2)){
    Smat[a,a] <- 1
    grad.phi <- sum(diag(2*(Lambdanew%*%t(Lambdanew) + Phiold - D)%*%Smat))
    mnew[a] <- beta1*mold[a] + (1-beta1)*grad.phi
    unew[a] <- max(beta2*uold[a], abs(grad.phi))
    Phinew[a,a] <- Phiold[a,a] - alpha/(1 - beta1^k)*mnew[a]/unew[a]
    Smat <- matrix(0, nrow = p1+p2, ncol = p1+p2)
  }
  
  k <- k+1
}

Phinew[Phinew < 0] <- 1

Gammaest <- Lambdanew[,1:j]
phi3fih <- diag(diag(Gammaest[(1:j), (1:j)]))
phi3fi <- phi3fih^2
Gammaest <- Gammaest%*%phi3fih
Gamma1est <- Gammaest[1:p1,]
Gamma2est <- Gammaest[(p1+1):(p1+p2),]

MSEg1fi <- 1/((p1*j)-j*(j-1)/2)*sum((Gamma1est - Gamma1)^2)
MSEg2fi <- 1/(p2*j)*sum((Gamma2est - Gamma2)^2)

Delta1est <- Lambdanew[(1:p1),j+(1:z1)]
phi41fih <- diag(diag(Delta1est[(1:z1), (1:z1)]))
phi41fi <- phi41fih^2
Delta1est <- Delta1est%*%phi41fih

MSEd1fi <- 1/((p1*z1)-z1*(z1-1)/2)*sum((Delta1est - Delta1)^2)

Delta2est <- Lambdanew[p1+(1:p2),j+z1+(1:z2)]
phi42fih <- diag(diag(Delta2est[(1:z2), (1:z2)]))
phi42fi <- phi42fih^2
Delta2est <- Delta2est%*%phi42fih

MSEd2fi <- 1/((p2*z2)-z2*(z2-1)/2)*sum((Delta2est - Delta2)^2)

Phi21fi <- diag(diag(Phinew)[1:p1])
Phi22fi <- diag(diag(Phinew)[p1+(1:p2)])

MSEphi21fi <- 1/p1*sum((Phi21fi - Phi21)^2)
MSEphi22fi <- 1/p2*sum((Phi22fi - Phi22)^2)

MSEphi3fi <- 1/m*sum((phi3fi - Phi3)^2)
MSEphi41fi <- 1/u1*sum((phi41fi - Phi41)^2)
MSEphi42fi <- 1/u2*sum((phi42fi - Phi42)^2)

PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)

C <- rbind(cbind(A1%*%PhiR, B1%*%PhiS1S), cbind(A2, B2%*%PhiS2S))

Ga <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(t(X1),t(Y1)), cov(t(X2),t(Y1)))
Gb <- solve(t(C)%*%C)%*%t(C)%*%rbind(cov(t(X1),t(Y2)), cov(t(X2),t(Y2)))
G1 <- t(Ga[1:m,])
G2 <- t(Gb[1:m,])
G3 <- t(Ga[(m+1):(m+u),])
G4 <- t(Gb[(m+1):(m+u),])

## Getting the loading matrix for the latent model connecting X and Y 
Thetaest <- 0.5*(solve(t(Gamma1est)%*%Gamma1est)%*%t(Gamma1est)%*%G1 + solve(t(Gamma2est)%*%Gamma2est)%*%t(Gamma2est)%*%G2)
Pi1est <- solve(t(Delta1est)%*%Delta1est)%*%t(Delta1est)%*%(G1 - Gamma1est%*%Thetaest)
Pi2est <- solve(t(Delta2est)%*%Delta2est)%*%t(Delta2est)%*%(G2 - Gamma2est%*%Thetaest)

## Getting the loading matrix for the latent model connecting X and Y 
Psiest <- 0.5*(solve(t(Gamma1est)%*%Gamma1est)%*%t(Gamma1est)%*%G3 + solve(t(Gamma2est)%*%Gamma2est)%*%t(Gamma2est)%*%G4)
Omega1est <- solve(t(Delta1est)%*%Delta1est)%*%t(Delta1est)%*%(G3 - Gamma1est%*%Psiest)
Omega2est <- solve(t(Delta2est)%*%Delta2est)%*%t(Delta2est)%*%(G4 - Gamma2est%*%Psiest)

MSEthest <- 1/(m*j)*sum((Theta-Thetaest)^2)
MSEpi1est <- 1/(z1*m)*sum((Pi1-Pi1est)^2)
MSEpi2est <- 1/(z2*m)*sum((Pi2-Pi2est)^2)
MSEpsiest <- 1/(j*(u1+u2))*sum((Psi-Psiest)^2)
MSEom1est <- 1/(z1*(u1+u2))*sum((Omega1-Omega1est)^2)
MSEom2est <- 1/(z2*(u1+u2))*sum((Omega2-Omega2est)^2)


### The EM step for the Y part, given initial values from the ADAM max 
###  However, the initials from the ADAM max sometimes not good, and hence EM did not improve much ...  

EMAlgYAdLassoCV <- function(Xpanel = list(), Ypanel = list(), Xests = list(), Gammainits = list(), Deltainits = list(), latentinits = list(), Xvcov = list(), vcovinits = list(), tuningpGamma = seq(5.1,5.5,0.1), tuningpDelta = seq(5.1,5.5,0.1), weightsGamma = list(), weightsDelta = list()){
  alg.start <- Sys.time()
  
  #define initial values for matrices based on input and extract matrix dimensions
  A1 <- Xests[[1]]
  A2 <- Xests[[2]]
  B1 <- Xests[[3]]
  B2 <- Xests[[4]]
  niter <- 0
  n <- nrow(Ypanel[[1]])
  q1 <- nrow(A1)
  q2 <- nrow(A2)
  m <- ncol(A1)
  u1 <- ncol(B1)
  u2 <- ncol(B2)
  p1 <- nrow(Gammainits[[1]])
  p2 <- nrow(Gammainits[[2]])
  z1 <- ncol(Deltainits[[1]])
  z2 <- ncol(Deltainits[[2]])
  y_joint <- ncol(Gammainits[[1]])
  Gamma1 <- matrix(0, nrow = p1, ncol = y_joint)
  Gamma2 <- matrix(0, nrow = p2, ncol = y_joint)
  Delta1 <- matrix(0, nrow = p1, ncol = z1)
  Delta2 <- matrix(0, nrow = p2, ncol = z2)
  Gamma1new <- Gammainits[[1]]
  Gamma2new <- Gammainits[[2]]
  Delta1new <- Deltainits[[1]]
  Delta2new <- Deltainits[[2]]
  Theta <- matrix(0, nrow = y_joint, ncol = m)
  Pi1 <- matrix(0, nrow = z1, ncol = m)
  Pi2 <- matrix(0, nrow = z2, ncol = m)
  Psi <- matrix(0, nrow = y_joint, ncol = u1 + u2)
  Omega1 <- matrix(0, nrow = z1, ncol = u1 + u2)
  Omega2 <- matrix(0, nrow = z2, ncol = u1 + u2)
  Thetanew <- latentinits[[1]]
  Pi1new <- latentinits[[2]]
  Pi2new <- latentinits[[3]]
  Psinew <- latentinits[[4]]
  Omega1new <- latentinits[[5]]
  Omega2new <- latentinits[[6]]
  PhiR <- Xvcov[[1]]
  PhiS1 <- Xvcov[[2]]
  PhiS2 <- Xvcov[[3]]
  Phi11 <- Xvcov[[4]]
  Phi12 <- Xvcov[[5]]
  Phi21 <- matrix(0, nrow = p1, ncol = p1)
  Phi22 <- matrix(0, nrow = p2, ncol = p2)
  Phi3 <- matrix(0, nrow = y_joint, ncol = y_joint)
  Phi41 <- matrix(0, nrow = z1, ncol = z1)
  Phi42 <- matrix(0, nrow = z2, ncol = z2)
  Phi21new <- vcovinits[[1]]
  Phi22new <- vcovinits[[2]]
  Phi3new <- vcovinits[[3]]
  Phi41new <- vcovinits[[4]]
  Phi42new <- vcovinits[[5]]
  
  #Used for the covariance between Si1 and Si, where Si = [Si1, Si2]. See derivations for details
  PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
  PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
  PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))
  
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  
  while(((norm(Gamma1 - Gamma1new) > 0.001) | (norm(Gamma2 - Gamma2new) > 0.001) | (norm(Delta1 - Delta1new) > 0.001) | (norm(Delta2 - Delta2new) > 0.001) | (norm(Theta - Thetanew) > 0.001) | (norm(Pi1 - Pi1new, type = "F") > 0.001) | (norm(Psi - Psinew) > 0.001) | (norm(Omega1 - Omega1new) > 0.001) | (norm(Omega2 - Omega2new) > 0.001) | (norm(Phi21 - Phi21new) > 0.001) | (norm(Phi22 - Phi22new) > 0.001) | (norm(Phi3 - Phi3new) > 0.001) |  (norm(Phi41 - Phi41new) > 0.001) | (norm(Phi42 - Phi42new) > 0.001)) & (niter < 30000)){
    
    Gamma1 <- Gamma1new
    Gamma2 <- Gamma2new
    Delta1 <- Delta1new
    Delta2 <- Delta2new
    Theta <- Thetanew
    Pi1 <- Pi1new
    Pi2 <- Pi2new
    Psi <- Psinew
    Omega1 <- Omega1new
    Omega2 <- Omega2new
    Phi21 <- Phi21new
    Phi22 <- Phi22new
    Phi3 <- Phi3new
    Phi41 <- Phi41new
    Phi42 <- Phi42new
    
    #for the first iteration of the algorithm ONLY, use CV to calculate the optimal tuning parameter for B
    if(niter == 0){
      #initialize norms for convergence, initial values for cross-validation, and
      #the coordinate descent update parameters for each matrix (epsilonijG1, thetaijG1, etc.)
      
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
      Gamma1CV <- matrix(0, nrow = p1, ncol = y_joint)
      Gamma2CV <- matrix(0, nrow = p2, ncol = y_joint)
      Delta1CV <- matrix(0, nrow = p1, ncol = z1)
      Delta2CV <- matrix(0, nrow = p2, ncol = z2)
      ThetaCV <- matrix(0, nrow = y_joint, ncol = m)
      Pi1CV <- matrix(0, nrow = z1, ncol = m)
      Pi2CV <- matrix(0, nrow = z2, ncol = m)
      PsiCV <- matrix(0, nrow = y_joint, ncol = u1 + u2)
      Omega1CV <- matrix(0, nrow = z1, ncol = u1 + u2)
      Omega2CV <- matrix(0, nrow = z2, ncol = u1 + u2)
      Phi21CV <- matrix(0, nrow = p1, ncol = p1)
      Phi22CV <- matrix(0, nrow = p2, ncol = p2)
      Phi3CV <- matrix(0, nrow = y_joint, ncol = y_joint)
      Phi41CV <- matrix(0, nrow = z1, ncol = z1)
      Phi42CV <- matrix(0, nrow = z2, ncol = z2)
      Gamma1start <- Gamma1
      Gamma2start <- Gamma2
      Delta1start <- Delta1
      Delta2start <- Delta2
      Thetastart <- Theta
      Pi1start <- Pi1
      Pi2start <- Pi2
      Psistart <- Psi
      Omega1start <- Omega1
      Omega2start <- Omega2
      Phi21start <- Phi21
      Phi22start <- Phi22
      Phi3start <- Phi3
      Phi41start <- Phi41
      Phi42start <- Phi42
      maxit <- 0
      
      #create matrix to store BIC for all combinations of tuning parameters
      BICmat <- matrix(0, nrow = length(tuningpGamma), ncol = length(tuningpDelta))
      
      cv.start <- Sys.time()
      
      for(a in 1:length(tuningpGamma)){
        for(b in length(tuningpDelta)){
          
          #reset norms for convergence and initial values for cross-validation for each new
          #combination of tuning parameters
          
          Gamma1CV <- Gamma1start
          Gamma2CV <- Gamma2start
          Delta1CV <- Delta1start
          Delta2CV <- Delta2start
          ThetaCV <- Thetastart
          Pi1CV <- Pi1start
          Pi2CV <- Pi2start
          PsiCV <- Psistart
          Omega1CV <- Omega1start
          Omega2CV <- Omega2start
          Phi21CV <- Phi21start
          Phi22CV <- Phi22start
          Phi3CV <- Phi3start
          Phi41CV <- Phi41start
          Phi42CV <- Phi42start
          
          
          Gamma1 <- matrix(0, nrow = p1, ncol = y_joint)
          Gamma2 <- matrix(0, nrow = p2, ncol = y_joint)
          Delta1 <- matrix(0, nrow = p1, ncol = z1)
          Delta2 <- matrix(0, nrow = p2, ncol = z2)
          Theta <- matrix(0, nrow = y_joint, ncol = m)
          Pi1 <- matrix(0, nrow = z1, ncol = m)
          Pi2 <- matrix(0, nrow = z2, ncol = m)
          Psi <- matrix(0, nrow = y_joint, ncol = u1 + u2)
          Omega1 <- matrix(0, nrow = z1, ncol = u1 + u2)
          Omega2 <- matrix(0, nrow = z2, ncol = u1 + u2)
          Phi21 <- matrix(0, nrow = p1, ncol = p1)
          Phi22 <- matrix(0, nrow = p2, ncol = p2)
          Phi3 <- matrix(0, nrow = y_joint, ncol = y_joint)
          Phi41 <- matrix(0, nrow = z1, ncol = z1)
          Phi42 <- matrix(0, nrow = z2, ncol = z2)
          
          maxit <- 0
          
          while(((norm(Gamma1 - Gamma1CV) > 0.001) | (norm(Gamma2 - Gamma2CV) > 0.001) | (norm(Delta1 - Delta1CV) > 0.001) | (norm(Delta2 - Delta2CV) > 0.001) | (norm(Theta - ThetaCV) > 0.001) | (norm(Pi1 - Pi1CV, type = "F") > 0.001) | (norm(Psi - PsiCV) > 0.001) | (norm(Omega1 - Omega1CV) > 0.001) | (norm(Omega2 - Omega2CV) > 0.001) | (norm(Phi21 - Phi21CV) > 0.001) | (norm(Phi22 - Phi22CV) > 0.001) | (norm(Phi3 - Phi3CV) > 0.001) |  (norm(Phi41 - Phi41CV) > 0.001) | (norm(Phi42 - Phi42CV) > 0.001)) & (maxit < 30000)){
            
            Gamma1 <- Gamma1CV
            Gamma2 <- Gamma2CV
            Delta1 <- Delta1CV
            Delta2 <- Delta2CV
            Theta <- ThetaCV
            Pi1 <- Pi1CV
            Pi2 <- Pi2CV
            Psi <- PsiCV
            Omega1 <- Omega1CV
            Omega2 <- Omega2CV
            Phi21 <- Phi21CV
            Phi22 <- Phi22CV
            Phi3 <- Phi3CV
            Phi41 <- Phi41CV
            Phi42 <- Phi42CV
            
            #calculate covariance of latent factors with observed data, i.e. Cov(Ei, Fi), where
            #E = [Ri, Si1, Si2, Vi, Wi1, Wi2] and F = [Xi1, Xi2, Yi1, Yi2]
            
            sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2), PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1), PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2)),
                             cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS1S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS1S))),
                             cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS2S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS2S))),
                             cbind(t(A1%*%PhiR%*%t(Theta) + B1%*%PhiS1S%*%t(Psi)), t(A2%*%PhiR%*%t(Theta) + B2%*%PhiS2S%*%t(Psi)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Theta) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Psi) + Gamma1%*%Phi3), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Theta) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Psi) + Gamma2%*%Phi3)),
                             cbind(t(A1%*%PhiR%*%t(Pi1) + B1%*%PhiS1S%*%t(Omega1)), t(A2%*%PhiR%*%t(Pi1) + B2%*%PhiS2S%*%t(Omega1)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega1) + Delta1%*%Phi41), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi1) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega1))),
                             cbind(t(A1%*%PhiR%*%t(Pi2) + B1%*%PhiS1S%*%t(Omega2)), t(A2%*%PhiR%*%t(Pi2) + B2%*%PhiS2S%*%t(Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega2)), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega2) + Delta2%*%Phi42)))
            
            #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
            
            sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                             cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                             cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)),
                             cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
            
            #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
            
            sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2), PhiR%*%t(Theta), PhiR%*%t(Pi1), PhiR%*%t(Pi2)),
                             cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2), PhiS1S%*%t(Psi), PhiS1S%*%t(Omega1), PhiS1S%*%t(Omega2)),
                             cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2, PhiS2S%*%t(Psi), PhiS2S%*%t(Omega1), PhiS2S%*%t(Omega2)),
                             cbind(Theta%*%PhiR, Psi%*%t(PhiS1S), Psi%*%t(PhiS2S), Theta%*%PhiR%*%t(Theta) + Psi%*%PhiS%*%t(Psi) + Phi3, Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1), Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)),
                             cbind(Pi1%*%PhiR, Omega1%*%t(PhiS1S), Omega1%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1)), Pi1%*%PhiR%*%t(Pi1) + Omega1%*%PhiS%*%t(Omega1) + Phi41, Pi1%*%PhiR%*%t(Pi2) + Omega1%*%PhiS%*%t(Omega2)),
                             cbind(Pi2%*%PhiR, Omega2%*%t(PhiS1S), Omega2%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)), Pi2%*%PhiR%*%t(Pi1) + Omega2%*%PhiS%*%t(Omega1), Pi2%*%PhiR%*%t(Pi2) + Omega2%*%PhiS%*%t(Omega2) + Phi42))
            
            #calculate conditional covariance of latent factors given observed data, i.e. Cov(Ei|Fi)
            #and then extract the conditional covariances of each latent factor (Cov(Ri|Fi), Cov(Si1|Fi), etc)
            
            condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
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
            
            #calculate conditional expectations of latent factors given observed data
            
            ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
            ER <- ERSVW[1:m,]
            ES1 <- ERSVW[m+(1:u1),]
            ES2 <- ERSVW[m+u1+(1:u2),]
            ES <- ERSVW[m+(1:(u1+u2)),]
            EV <- ERSVW[m+u1+u2+(1:y_joint),]
            EW1 <- ERSVW[m+u1+u2+y_joint+(1:z1),]
            EW2 <- ERSVW[m+u1+u2+y_joint+z1+(1:z2),]
            
            #coordinate descent for Gamma1 (note that y_joint = j, i.e. the number of joint
            #latent factors for Y. However, as j is being used to index for the loop, we instead
            #use the notation y_joint)
            
            Gamma1update <- matrix(0, nrow = p1, ncol = y_joint)
            Gamma2update <- matrix(0, nrow = p2, ncol = y_joint)
            Delta1update <- matrix(0, nrow = p1, ncol = z1)
            Delta2update <- matrix(0, nrow = p2, ncol = z2)
            
            while((norm(Gamma1update - Gamma1CV, type = "F") > 0.005) | (norm(Gamma2update - Gamma2CV, type = "F") > 0.005) | (norm(Delta1update - Delta1CV) > 0.005) | (norm(Delta2update - Delta2CV) > 0.005)){
              
              Gamma1update <- Gamma1CV
              Gamma2update <- Gamma2CV
              Delta1update <- Delta1CV
              Delta2update <- Delta2CV
              
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
                    thetaijG1 <- crossprod(EV[j,],Ypanel[[1]][,i])
                    omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1CV[i,-j]
                    for(k in 1:z1){
                      tauijG1 <- tauijG1 + Delta1CV[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
                    }
                    g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
                    if(weightsGamma[[1]][i,j] == 0){
                      lambdaG <- (tuningpGamma[a]*Phi21CV[i,i])/epsilonijG1*(1/1e-5)
                    }
                    else{
                      lambdaG <- (tuningpGamma[a]*Phi21CV[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
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
                  thetaijG1 <- crossprod(EV[j,],Ypanel[[1]][,i])
                  omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1CV[i,-j]
                  for(k in 1:z1){
                    tauijG1 <- tauijG1 + Delta1CV[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
                  }
                  g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
                  if(weightsGamma[[1]][i,j] == 0){
                    lambdaG <- (tuningpGamma[a]*Phi21CV[i,i])/epsilonijG1*(1/1e-5)
                  }
                  else{
                    lambdaG <- (tuningpGamma[a]*Phi21CV[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
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
                  thetaijG2 <- crossprod(EV[j,],Ypanel[[2]][,i])
                  omegaijG2 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma2CV[i,-j]
                  for(k in 1:z2){
                    tauijG2 <- tauijG2 + Delta2CV[i,k]*t(EV%*%t(EW2) + n*condvarVW2)[k,j]
                  }
                  g2bar <- ((thetaijG2-omegaijG2-tauijG2))/(epsilonijG2)
                  if(weightsGamma[[2]][i,j] == 0){
                    lambdaG <- (tuningpGamma[a]*Phi22CV[i,i])/epsilonijG2*(1/1e-5)
                  }
                  else{
                    lambdaG <-(tuningpGamma[a]*Phi22CV[i,i])/epsilonijG2*(1/abs(weightsGamma[[2]][i,j]))
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
                    thetaijD1 <- crossprod(EW1[j,],Ypanel[[1]][,i])
                    omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1CV[i,-j]
                    for(k in 1:y_joint){
                      tauijD1 <- tauijD1 + Gamma1CV[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
                    }
                    d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
                    if(weightsDelta[[1]][i,j] == 0){
                      lambdaD <- (tuningpDelta[b]*Phi21CV[i,i])/epsilonijD1*(1/1e-5)
                    }
                    else{
                      lambdaD <- (tuningpDelta[b]*Phi21CV[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
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
                  thetaijD1 <- crossprod(EW1[j,],Ypanel[[1]][,i])
                  omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1CV[i,-j]
                  for(k in 1:y_joint){
                    tauijD1 <- tauijD1 + Gamma1CV[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
                  }
                  d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
                  if(weightsDelta[[1]][i,j] == 0){
                    lambdaD <- (tuningpDelta[b]*Phi21CV[i,i])/epsilonijD1*(1/1e-5)
                  }
                  else{
                    lambdaD <- (tuningpDelta[b]*Phi21CV[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
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
                    thetaijD2 <- crossprod(EW2[j,],Ypanel[[2]][,i])
                    omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2CV[i,-j]
                    for(k in 1:y_joint){
                      tauijD2 <- tauijD2 + Gamma2CV[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
                    }
                    d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
                    if(weightsDelta[[2]][i,j] == 0){
                      lambdaD <- (tuningpDelta[b]*Phi22CV[i,i])/epsilonijD2*(1/1e-5)
                    }
                    else{
                      lambdaD <- (tuningpDelta[b]*Phi22CV[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
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
                  thetaijD2 <- crossprod(EW2[j,],Ypanel[[2]][,i])
                  omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2CV[i,-j]
                  for(k in 1:y_joint){
                    tauijD2 <- tauijD2 + Gamma2CV[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
                  }
                  d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
                  if(weightsDelta[[2]][i,j] == 0){
                    lambdaD <- (tuningpDelta[b]*Phi22CV[i,i])/epsilonijD2*(1/1e-5)
                  }
                  else{
                    lambdaD <- (tuningpDelta[b]*Phi22CV[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
                  }
                  Delta2CV[i,j] <- lasso(d2bar, lambdaD)
                }
              }
            }
            
            ThetaCV <- (n*condvarVR + EV%*%t(ER) - n*Psi%*%t(condvarRS) - Psi%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
            Pi1CV <- (n*t(condvarRW1) + EW1%*%t(ER) - n*Omega1%*%t(condvarRS) - Omega1%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
            Pi2CV <- (n*t(condvarRW2) + EW2%*%t(ER) - n*Omega2%*%t(condvarRS) - Omega2%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
            PsiCV <- (n*condvarVS + EV%*%t(ES) - n*ThetaCV%*%condvarRS - ThetaCV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
            Omega1CV <- (n*t(condvarSW1) + EW1%*%t(ES) - n*Pi1CV%*%condvarRS - Pi1CV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
            Omega2CV <- (n*t(condvarSW2) + EW2%*%t(ES) - n*Pi2CV%*%condvarRS - Pi2CV%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
            
            Phi21CV <- 1/n*diag(diag(t(Ypanel[[1]])%*%Ypanel[[1]] - 2*t(Ypanel[[1]])%*%t(Gamma1CV%*%EV) - 2*t(Ypanel[[1]])%*%t(Delta1CV%*%EW1) + 2*n*Gamma1CV%*%condvarVW1%*%t(Delta1CV)
                                     + 2*Gamma1CV%*%EV%*%t(Delta1CV%*%EW1) + n*Gamma1CV%*%condvarV%*%t(Gamma1CV) + Gamma1CV%*%EV%*%t(Gamma1CV%*%EV) + n*Delta1CV%*%condvarW1%*%t(Delta1CV) + Delta1CV%*%EW1%*%t(Delta1CV%*%EW1)))
            Phi22CV <- 1/n*diag(diag(t(Ypanel[[2]])%*%Ypanel[[2]] - 2*t(Ypanel[[2]])%*%t(Gamma2CV%*%EV) - 2*t(Ypanel[[2]])%*%t(Delta2CV%*%EW2) + 2*n*Gamma2CV%*%condvarVW2%*%t(Delta2CV)
                                     + 2*Gamma2CV%*%EV%*%t(Delta2CV%*%EW2) + n*Gamma2CV%*%condvarV%*%t(Gamma2CV) + Gamma2CV%*%EV%*%t(Gamma2CV%*%EV) + n*Delta2CV%*%condvarW2%*%t(Delta2CV) + Delta2CV%*%EW2%*%t(Delta2CV%*%EW2)))
            Phi3CV <- 1/n*diag(diag(n*condvarV + EV%*%t(EV) - 2*n*condvarVR%*%t(ThetaCV) - 2*EV%*%t(ThetaCV%*%ER) - 2*n*condvarVS%*%t(PsiCV) - 2*EV%*%t(PsiCV%*%ES) + n*ThetaCV%*%condvarR%*%t(ThetaCV) + ThetaCV%*%ER%*%t(ThetaCV%*%ER)
                                    + 2*n*ThetaCV%*%condvarRS%*%t(PsiCV) + 2*ThetaCV%*%ER%*%t(PsiCV%*%ES) + n*PsiCV%*%condvarS%*%t(PsiCV) + PsiCV%*%ES%*%t(PsiCV%*%ES)))
            Phi41CV <- 1/n*diag(diag(n*condvarW1 + EW1%*%t(EW1) - 2*n*t(Pi1CV%*%condvarRW1) -2*EW1%*%t(Pi1CV%*%ER) -2*n*t(Omega1CV%*%condvarSW1) - 2*EW1%*%t(Omega1CV%*%ES) + n*Pi1CV%*%condvarR%*%t(Pi1CV) + Pi1CV%*%ER%*%t(Pi1CV%*%ER)
                                     + 2*n*Pi1CV%*%condvarRS%*%t(Omega1CV) + 2*Pi1CV%*%ER%*%t(Omega1CV%*%ES) + n*Omega1CV%*%condvarS%*%t(Omega1CV) + Omega1CV%*%ES%*%t(Omega1CV%*%ES)))
            Phi42CV <- 1/n*diag(diag(n*condvarW2 + EW2%*%t(EW2) - 2*n*t(Pi2CV%*%condvarRW2) -2*EW2%*%t(Pi2CV%*%ER) -2*n*t(Omega2CV%*%condvarSW2) - 2*EW2%*%t(Omega2CV%*%ES) + n*Pi2CV%*%condvarR%*%t(Pi2CV) + Pi2CV%*%ER%*%t(Pi2CV%*%ER)
                                     + 2*n*Pi2CV%*%condvarRS%*%t(Omega2CV) + 2*Pi2CV%*%ER%*%t(Omega2CV%*%ES) + n*Omega2CV%*%condvarS%*%t(Omega2CV) + Omega2CV%*%ES%*%t(Omega2CV%*%ES)))
            
            #once new updates for every matrix are calculated, calculate Frobenius norm of new matrices for comparison to norm of previous matrices
            #to check for convergence
            maxit <- maxit + 1
            log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]],Ypanel[[2]]))+(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))))) + (p1+p2)*log((2*pi))+log(det(sigma11)))
            print(log.lik)
          }
          
          #recalculate model covariance of observed data for calculation of likelihood for BIC
          
          sigma11special <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                                  cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                                  cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)),
                                  cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
          
          #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
          log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]],Ypanel[[2]]))+(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))))) + (p1+p2)*log((2*pi))+log(det(sigma11)))
          
          BICmat[a,b] <- log(n)*(sum(Gamma1CV != 0) + sum(Gamma2CV != 0) + sum(Delta1CV != 0) + sum(Delta2CV != 0) - (y_joint+z1+z2)) - 2*log.lik
          print("CV done")
        }
      }
      
      #select optimal tuning parameters for Gamma and Delta
      tuningpGammahat <- tuningpGamma[which.min(apply(BICmat, MARGIN = 1, min))]
      tuningpDeltahat <- tuningpDelta[which.min(apply(BICmat, MARGIN = 2, min))]
      
      Gamma1 <- Gamma1start
      Gamma2 <- Gamma2start
      Delta1 <- Delta1start
      Delta2 <- Delta2start
      Theta <- Thetastart
      Pi1 <- Pi1start
      Pi2 <- Pi2start
      Psi <- Psistart
      Omega1 <- Omega1start
      Omega2 <- Omega2start
      Phi21 <- Phi21start
      Phi22 <- Phi22start
      Phi3 <- Phi3start
      Phi41 <- Phi41start
      Phi42 <- Phi42start
      
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    #now that the optimal tuning parameters have been selected, simply perform coordinate descent
    #using the optimal tuning parameters
    sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2), PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1), PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2)),
                     cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS1S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS1S))),
                     cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS2S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS2S))),
                     cbind(t(A1%*%PhiR%*%t(Theta) + B1%*%PhiS1S%*%t(Psi)), t(A2%*%PhiR%*%t(Theta) + B2%*%PhiS2S%*%t(Psi)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Theta) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Psi) + Gamma1%*%Phi3), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Theta) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Psi) + Gamma2%*%Phi3)),
                     cbind(t(A1%*%PhiR%*%t(Pi1) + B1%*%PhiS1S%*%t(Omega1)), t(A2%*%PhiR%*%t(Pi1) + B2%*%PhiS2S%*%t(Omega1)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega1) + Delta1%*%Phi41), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi1) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega1))),
                     cbind(t(A1%*%PhiR%*%t(Pi2) + B1%*%PhiS1S%*%t(Omega2)), t(A2%*%PhiR%*%t(Pi2) + B2%*%PhiS2S%*%t(Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega2)), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega2) + Delta2%*%Phi42)))
    
    #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
    
    sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                     cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                     cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)),
                     cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
    
    #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
    
    sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2), PhiR%*%t(Theta), PhiR%*%t(Pi1), PhiR%*%t(Pi2)),
                     cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2), PhiS1S%*%t(Psi), PhiS1S%*%t(Omega1), PhiS1S%*%t(Omega2)),
                     cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2, PhiS2S%*%t(Psi), PhiS2S%*%t(Omega1), PhiS2S%*%t(Omega2)),
                     cbind(Theta%*%PhiR, Psi%*%t(PhiS1S), Psi%*%t(PhiS2S), Theta%*%PhiR%*%t(Theta) + Psi%*%PhiS%*%t(Psi) + Phi3, Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1), Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)),
                     cbind(Pi1%*%PhiR, Omega1%*%t(PhiS1S), Omega1%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1)), Pi1%*%PhiR%*%t(Pi1) + Omega1%*%PhiS%*%t(Omega1) + Phi41, Pi1%*%PhiR%*%t(Pi2) + Omega1%*%PhiS%*%t(Omega2)),
                     cbind(Pi2%*%PhiR, Omega2%*%t(PhiS1S), Omega2%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)), Pi2%*%PhiR%*%t(Pi1) + Omega2%*%PhiS%*%t(Omega1), Pi2%*%PhiR%*%t(Pi2) + Omega2%*%PhiS%*%t(Omega2) + Phi42))
    
    condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
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
    
    ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
    ER <- ERSVW[1:m,]
    ES1 <- ERSVW[m+(1:u1),]
    ES2 <- ERSVW[m+u1+(1:u2),]
    ES <- ERSVW[m+(1:(u1+u2)),]
    EV <- ERSVW[m+u1+u2+(1:y_joint),]
    EW1 <- ERSVW[m+u1+u2+y_joint+(1:z1),]
    EW2 <- ERSVW[m+u1+u2+y_joint+z1+(1:z2),]
    
    #now that the optimal tuning parameter has been chosen, run the EM algorithm as normal, using the optimal tuning parameter for all subsequent iterations
    
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
    
    Gamma1update <- matrix(0, nrow = p1, ncol = y_joint)
    Gamma2update <- matrix(0, nrow = p2, ncol = y_joint)
    Delta1update <- matrix(0, nrow = p1, ncol = z1)
    Delta2update <- matrix(0, nrow = p2, ncol = z2)
    
    while((norm(Gamma1update - Gamma1new, type = "F") > 0.005) | (norm(Gamma2update - Gamma2new, type = "F") > 0.005) | (norm(Delta1update - Delta1new) > 0.005) | (norm(Delta2update - Delta2new) > 0.005)){
      
      Gamma1update <- Gamma1new
      Gamma2update <- Gamma2new
      Delta1update <- Delta1new
      Delta2update <- Delta2new
      
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
            thetaijG1 <- crossprod(EV[j,],Ypanel[[1]][,i])
            omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1new[i,-j]
            for(k in 1:z1){
              tauijG1 <- tauijG1 + Delta1new[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
            }
            g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
            if(weightsGamma[[1]][i,j] == 0){
              lambdaG <- (tuningpGammahat*Phi21new[i,i])/epsilonijG1*(1/1e-5)
            }
            else{
              lambdaG <- (tuningpGammahat*Phi21new[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
            }
            Gamma1new[i,j] <- lasso(g1bar, lambdaG)
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
          thetaijG1 <- crossprod(EV[j,],Ypanel[[1]][,i])
          omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma1new[i,-j]
          for(k in 1:z1){
            tauijG1 <- tauijG1 + Delta1new[i,k]*t(EV%*%t(EW1) + n*condvarVW1)[k,j]
          }
          g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
          if(weightsGamma[[1]][i,j] == 0){
            lambdaG <- (tuningpGammahat*Phi21new[i,i])/epsilonijG1*(1/1e-5)
          }
          else{
            lambdaG <- (tuningpGammahat*Phi21new[i,i])/epsilonijG1*(1/abs(weightsGamma[[1]][i,j]))
          }
          Gamma1new[i,j] <- lasso(g1bar, lambdaG)
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
          thetaijG2 <- crossprod(EV[j,],Ypanel[[2]][,i])
          omegaijG2 <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%Gamma2new[i,-j]
          for(k in 1:z2){
            tauijG2 <- tauijG2 + Delta2new[i,k]*t(EV%*%t(EW2) + n*condvarVW2)[k,j]
          }
          g2bar <- ((thetaijG2-omegaijG2-tauijG2))/(epsilonijG2)
          if(weightsGamma[[2]][i,j] == 0){
            lambdaG <- (tuningpGammahat*Phi22new[i,i])/epsilonijG2*(1/1e-5)
          }
          else{
            lambdaG <-(tuningpGammahat*Phi22new[i,i])/epsilonijG2*(1/abs(weightsGamma[[2]][i,j]))
          }
          Gamma2new[i,j] <- lasso(g2bar, lambdaG)
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
            thetaijD1 <- crossprod(EW1[j,],Ypanel[[1]][,i])
            omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1new[i,-j]
            for(k in 1:y_joint){
              tauijD1 <- tauijD1 + Gamma1new[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
            }
            d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
            if(weightsDelta[[1]][i,j] == 0){
              lambdaD <- (tuningpDeltahat*Phi21new[i,i])/epsilonijD1*(1/1e-5)
            }
            else{
              lambdaD <- (tuningpDeltahat*Phi21new[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
            }
            Delta1new[i,j] <- lasso(d1bar, lambdaD)
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
          thetaijD1 <- crossprod(EW1[j,],Ypanel[[1]][,i])
          omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-j,j]%*%Delta1new[i,-j]
          for(k in 1:y_joint){
            tauijD1 <- tauijD1 + Gamma1new[i,k]*t(EW1%*%t(EV) + n*t(condvarVW1))[k,j]
          }
          d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
          if(weightsDelta[[1]][i,j] == 0){
            lambdaD <- (tuningpDeltahat*Phi21new[i,i])/epsilonijD1*(1/1e-5)
          }
          else{
            lambdaD <- (tuningpDeltahat*Phi21new[i,i])/epsilonijD1*(1/abs(weightsDelta[[1]][i,j]))
          }
          Delta1new[i,j] <- lasso(d1bar, lambdaD)
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
            thetaijD2 <- crossprod(EW2[j,],Ypanel[[2]][,i])
            omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2new[i,-j]
            for(k in 1:y_joint){
              tauijD2 <- tauijD2 + Gamma2new[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
            }
            d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
            if(weightsDelta[[2]][i,j] == 0){
              lambdaD <- (tuningpDeltahat*Phi22new[i,i])/epsilonijD2*(1/1e-5)
            }
            else{
              lambdaD <- (tuningpDeltahat*Phi22new[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
            }
            Delta2new[i,j] <- lasso(d2bar, lambdaD)
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
          thetaijD2 <- crossprod(EW2[j,],Ypanel[[2]][,i])
          omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-j,j]%*%Delta2new[i,-j]
          for(k in 1:y_joint){
            tauijD2 <- tauijD2 + Gamma2new[i,k]*t(EW2%*%t(EV) + n*t(condvarVW2))[k,j]
          }
          d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
          if(weightsDelta[[2]][i,j] == 0){
            lambdaD <- (tuningpDeltahat*Phi22new[i,i])/epsilonijD2*(1/1e-5)
          }
          else{
            lambdaD <- (tuningpDeltahat*Phi22new[i,i])/epsilonijD2*(1/abs(weightsDelta[[2]][i,j]))
          }
          Delta2new[i,j] <- lasso(d2bar, lambdaD)
        }
      }
    }
    
    Thetanew <- (n*condvarVR + EV%*%t(ER) - n*Psi%*%t(condvarRS) - Psi%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    Pi1new <- (n*t(condvarRW1) + EW1%*%t(ER) - n*Omega1%*%t(condvarRS) - Omega1%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    Pi2new <- (n*t(condvarRW2) + EW2%*%t(ER) - n*Omega2%*%t(condvarRS) - Omega2%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
    Psinew <- (n*condvarVS + EV%*%t(ES) - n*Theta%*%condvarRS - Theta%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    Omega1new <- (n*t(condvarSW1) + EW1%*%t(ES) - n*Pi1%*%condvarRS - Pi1%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    Omega2new <- (n*t(condvarSW2) + EW2%*%t(ES) - n*Pi2%*%condvarRS - Pi2%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
    
    Phi21new <- 1/n*diag(diag(t(Ypanel[[1]])%*%Ypanel[[1]] - 2*t(Ypanel[[1]])%*%t(Gamma1new%*%EV) - 2*t(Ypanel[[1]])%*%t(Delta1new%*%EW1) + 2*n*Gamma1new%*%condvarVW1%*%t(Delta1new)
                              + 2*Gamma1new%*%EV%*%t(Delta1new%*%EW1) + n*Gamma1new%*%condvarV%*%t(Gamma1new) + Gamma1new%*%EV%*%t(Gamma1new%*%EV) + n*Delta1new%*%condvarW1%*%t(Delta1new) + Delta1new%*%EW1%*%t(Delta1new%*%EW1)))
    Phi22new <- 1/n*diag(diag(t(Ypanel[[2]])%*%Ypanel[[2]] - 2*t(Ypanel[[2]])%*%t(Gamma2new%*%EV) - 2*t(Ypanel[[2]])%*%t(Delta2new%*%EW2) + 2*n*Gamma2new%*%condvarVW2%*%t(Delta2new)
                              + 2*Gamma2new%*%EV%*%t(Delta2new%*%EW2) + n*Gamma2new%*%condvarV%*%t(Gamma2new) + Gamma2new%*%EV%*%t(Gamma2new%*%EV) + n*Delta2new%*%condvarW2%*%t(Delta2new) + Delta2new%*%EW2%*%t(Delta2new%*%EW2)))
    Phi3new <- 1/n*diag(diag(n*condvarV + EV%*%t(EV) -2*n*(condvarVR%*%t(Thetanew)) - 2*(EV%*%t(Thetanew%*%ER)) - 2*n*(condvarVS%*%t(Psinew)) - 2*(EV%*%t(Psinew%*%ES)) + n*Thetanew%*%condvarR%*%t(Thetanew) + Thetanew%*%ER%*%t(Thetanew%*%ER)
                             + 2*n*Thetanew%*%condvarRS%*%t(Psinew) + 2*Thetanew%*%ER%*%t(Psinew%*%ES) + n*Psinew%*%condvarS%*%t(Psinew) + Psinew%*%ES%*%t(Psinew%*%ES)))
    Phi41new <- 1/n*diag(diag(n*condvarW1 + EW1%*%t(EW1) - 2*n*t(Pi1new%*%condvarRW1) -2*EW1%*%t(Pi1new%*%ER) -2*n*t(Omega1new%*%condvarSW1) - 2*EW1%*%t(Omega1new%*%ES) + n*Pi1new%*%condvarR%*%t(Pi1new) + Pi1new%*%ER%*%t(Pi1new%*%ER)
                              + 2*n*Pi1new%*%condvarRS%*%t(Omega1new) + 2*Pi1new%*%ER%*%t(Omega1new%*%ES) + n*Omega1new%*%condvarS%*%t(Omega1new) + Omega1new%*%ES%*%t(Omega1new%*%ES)))
    Phi42new <- 1/n*diag(diag(n*condvarW2 + EW2%*%t(EW2) - 2*n*t(Pi2new%*%condvarRW2) -2*EW2%*%t(Pi2new%*%ER) -2*n*t(Omega2new%*%condvarSW2) - 2*EW2%*%t(Omega2new%*%ES) + n*Pi2new%*%condvarR%*%t(Pi2new) + Pi2new%*%ER%*%t(Pi2new%*%ER)
                              + 2*n*Pi2new%*%condvarRS%*%t(Omega2new) + 2*Pi2new%*%ER%*%t(Omega2new%*%ES) + n*Omega2new%*%condvarS%*%t(Omega2new) + Omega2new%*%ES%*%t(Omega2new%*%ES)))
    
    niter <- niter + 1
    print(-(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]],Ypanel[[2]]))+(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))))) + (p1+p2)*log((2*pi))+log(det(sigma11))))
  }
  
  
  #extract final matrix estmates, conditional expectations, etc
  Gamma1 <- Gamma1new
  Gamma2 <- Gamma2new
  Delta1 <- Delta1new
  Delta2 <- Delta2new
  Theta <- Thetanew
  Pi1 <- Pi1new
  Pi2 <- Pi2new
  Psi <- Psinew
  Omega1 <- Omega1new
  Omega2 <- Omega2new
  Phi21 <- Phi21new
  Phi22 <- Phi22new
  Phi3 <- Phi3new
  Phi41 <- Phi41new
  Phi42 <- Phi42new
  
  sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2), PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1), PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2)),
                   cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS1S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS1S))),
                   cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS2S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS2S))),
                   cbind(t(A1%*%PhiR%*%t(Theta) + B1%*%PhiS1S%*%t(Psi)), t(A2%*%PhiR%*%t(Theta) + B2%*%PhiS2S%*%t(Psi)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Theta) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Psi) + Gamma1%*%Phi3), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Theta) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Psi) + Gamma2%*%Phi3)),
                   cbind(t(A1%*%PhiR%*%t(Pi1) + B1%*%PhiS1S%*%t(Omega1)), t(A2%*%PhiR%*%t(Pi1) + B2%*%PhiS2S%*%t(Omega1)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega1) + Delta1%*%Phi41), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi1) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega1))),
                   cbind(t(A1%*%PhiR%*%t(Pi2) + B1%*%PhiS1S%*%t(Omega2)), t(A2%*%PhiR%*%t(Pi2) + B2%*%PhiS2S%*%t(Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega2)), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega2) + Delta2%*%Phi42)))
  
  #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
  
  sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                   cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                   cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)),
                   cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
  
  #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
  
  sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2), PhiR%*%t(Theta), PhiR%*%t(Pi1), PhiR%*%t(Pi2)),
                   cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2), PhiS1S%*%t(Psi), PhiS1S%*%t(Omega1), PhiS1S%*%t(Omega2)),
                   cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2, PhiS2S%*%t(Psi), PhiS2S%*%t(Omega1), PhiS2S%*%t(Omega2)),
                   cbind(Theta%*%PhiR, Psi%*%t(PhiS1S), Psi%*%t(PhiS2S), Theta%*%PhiR%*%t(Theta) + Psi%*%PhiS%*%t(Psi) + Phi3, Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1), Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)),
                   cbind(Pi1%*%PhiR, Omega1%*%t(PhiS1S), Omega1%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1)), Pi1%*%PhiR%*%t(Pi1) + Omega1%*%PhiS%*%t(Omega1) + Phi41, Pi1%*%PhiR%*%t(Pi2) + Omega1%*%PhiS%*%t(Omega2)),
                   cbind(Pi2%*%PhiR, Omega2%*%t(PhiS1S), Omega2%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)), Pi2%*%PhiR%*%t(Pi1) + Omega2%*%PhiS%*%t(Omega1), Pi2%*%PhiR%*%t(Pi2) + Omega2%*%PhiS%*%t(Omega2) + Phi42))
  
  condvar <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)
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
  
  ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
  ER <- ERSVW[1:m,]
  ES1 <- ERSVW[m+(1:u1),]
  ES2 <- ERSVW[m+u1+(1:u2),]
  ES <- ERSVW[m+(1:(u1+u2)),]
  EV <- ERSVW[m+u1+u2+(1:y_joint),]
  EW1 <- ERSVW[m+u1+u2+y_joint+(1:z1),]
  EW2 <- ERSVW[m+u1+u2+y_joint+z1+(1:z2),]
  
  log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]],Ypanel[[2]]))+(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))%*%t(colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])) - colMeans(cbind(Xpanel[[1]], Xpanel[[2]], Ypanel[[1]], Ypanel[[2]])))))) + (p1+p2)*log((2*pi))+log(det(sigma11)))
  
  BICopt <- log(n)*(sum(Gamma1 != 0) + sum(Gamma2 != 0) + sum(Delta1 != 0) + sum(Delta2 != 0) - (y_joint+z1+z2)) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("Gamma1" = Gamma1, "Gamma2" = Gamma2, "Delta1" = Delta1, "Delta2" = Delta2, "Theta" = Theta, "Pi1" = Pi1, "Pi2" = Pi2, "Psi" = Psi, "Omega1" = Omega1, "Omega2" = Omega2, "Phi21" = Phi21, "Phi22" = Phi22, "Phi3" = Phi3, "Phi41" = Phi41, "Phi42" = Phi42, "iterations" = niter, "optimal lambdaGamma" = tuningpGammahat, "optimal lambdaDelta" = tuningpDeltahat, "ERSVW" = ERSVW, "BICopt" = BICopt)
  return(output)
}

test2 <- EMAlgYAdLassoCV(Xpanel = list(t(X1), t(X2)), Ypanel = list(t(Y1), t(Y2)), Xests = list(A1em, A2em, B1em, B2em), Gammainits = list(Gamma1est, Gamma2est), Deltainits = list(Delta1est, Delta2est), latentinits = list(Thetaest, Pi1est, Pi2est, Psiest, Omega1est, Omega2est), Xvcov = list(PhiRem, PhiS1em, PhiS2em, Phi11em, Phi12em), vcovinits = list(Phi21fi, Phi22fi, phi3fi, phi41fi, phi42fi), tuningpGamma = 0, tuningpDelta = 0, weightsGamma = list(Gamma1est, Gamma2est), weightsDelta = list(Delta1est, Delta2est))

#test3 <- EMAlgYAdLassoCV(Xpanel = list(t(X1), t(X2)), Ypanel = list(t(Y1), t(Y2)), Xests = list(A1, A2, B1, B2), Gammainits = list(Gamma1, Gamma2), Deltainits = list(Delta1, Delta2), latentinits = list(Theta, Pi1, Pi2, Psi, Omega1, Omega2), Xvcov = list(PhiR, PhiS1, PhiS2, Phi11, Phi12), vcovinits = list(Phi21, Phi22, Phi3, Phi41, Phi42), tuningpGamma = 0, tuningpDelta = 0, weightsGamma = list(Gamma1, Gamma2), weightsDelta = list(Delta1, Delta2))

MSEg1em <- 1/(p1*j-j*(j+1)/2)*sum((Gamma1 - test2$Gamma1)^2)
MSEg2em <- 1/(p2*j-j*(j+1)/2)*sum((Gamma2 - test2$Gamma2)^2)
MSEd1em <- 1/(p1*z1-z1*(z1+1)/2)*sum((Delta1 - test2$Delta1)^2)
MSEd2em <- 1/(p2*z2-z2*(z2+1)/2)*sum((Delta2 - test2$Delta2)^2)
MSEthem <- 1/(m*j)*sum((Theta-test2$Theta)^2)
MSEpi1em <- 1/(z1*m)*sum((Pi1-test2$Pi1)^2)
MSEpi2em <- 1/(z2*m)*sum((Pi2-test2$Pi2)^2)
MSEpsem <- 1/(j*(u1+u2))*sum((Psi-test2$Psi)^2)
MSEom1em <- 1/(z1*(u1+u2))*sum((Omega1-test2$Omega1)^2)
MSEom2em <- 1/(z2*(u1+u2))*sum((Omega2-test2$Omega2)^2)
MSEphi21em <- 1/p1*sum((Phi21-test2$Phi21)^2)
MSEphi22em <- 1/p2*sum((Phi22-test2$Phi22)^2)
MSEphi3em <- 1/j*sum((Phi3-test2$Phi3)^2)
MSEphi41em <- 1/z1*sum((Phi41-test2$Phi41)^2)
MSEphi42em <- 1/z2*sum((Phi42-test2$Phi42)^2)

Results <- list(test, test2)

saveRDS(Results, file=paste0("AdamAltEstrep", datset, ".rds"))
