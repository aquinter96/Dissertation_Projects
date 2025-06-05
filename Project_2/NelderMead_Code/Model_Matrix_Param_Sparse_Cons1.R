## Model_Matrix_Param_Sparse_cons1.R

##### parameters on dimensions of X, Y, and latent factors  

Model_Params <- list()

############################################  
###########  X part parameters A1/A2/B1/B2 ###########  
############################################  

## A1: q1*m: the factor loading matrix corresponds to the joint latent factor for X1
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis), for A = (A1': A2')'
Model_Params$A1 = matrix(0, nrow = meta_param$q1, ncol = meta_param$m)
diag(Model_Params$A1[(1:meta_param$m),(1:meta_param$m)]) = 1
Model_Params$A1[(1:meta_param$m),(1:meta_param$m)][lower.tri(Model_Params$A1[(1:meta_param$m),(1:meta_param$m)], diag = F)] = rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)*(2*rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)-1)*sample(c(runif(meta_param$m*(meta_param$m-1)/2,0.38,0.42),runif(meta_param$m*(meta_param$m-1)/2,0.48,0.52),runif(meta_param$m*(meta_param$m-1)/2,0.98,1.02)),meta_param$m*(meta_param$m-1)/2)
Model_Params$A1[(meta_param$m+1):meta_param$q1,] = rbinom(meta_param$q1*meta_param$m-meta_param$m^2,1,0.5)*(2*rbinom(meta_param$q1*meta_param$m-meta_param$m^2,1,0.5)-1)*sample(c(runif(meta_param$q1*meta_param$m-meta_param$m^2,0.38,0.42),runif(meta_param$q1*meta_param$m-meta_param$m^2,0.48,0.52),runif(meta_param$q1*meta_param$m-meta_param$m^2,0.98,1.02)),meta_param$q1*meta_param$m-meta_param$m^2)

## B1: q1*u1: the factor loading matrix corresponds to the separate latent factor for X1
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$B1 = matrix(0, nrow = meta_param$q1, ncol = meta_param$u1)
diag(Model_Params$B1[(1:meta_param$u1),(1:meta_param$u1)]) = 1
Model_Params$B1[(1:meta_param$u1),(1:meta_param$u1)][lower.tri(Model_Params$B1[(1:meta_param$u1),(1:meta_param$u1)], diag = F)] = rbinom(meta_param$u1*(meta_param$u1-1)/2,1,0.5)*(2*rbinom(meta_param$u1*(meta_param$u1-1)/2,1,0.5)-1)*sample(c(runif(meta_param$u1*(meta_param$u1-1)/2,0.38,0.42),runif(meta_param$u1*(meta_param$u1-1)/2,0.48,0.52),runif(meta_param$u1*(meta_param$u1-1)/2,0.98,1.02)),meta_param$u1*(meta_param$u1-1)/2)
Model_Params$B1[(meta_param$u1+1):meta_param$q1,] = rbinom(meta_param$q1*meta_param$u1-meta_param$u1^2,1,0.5)*(2*rbinom(meta_param$q1*meta_param$u1-meta_param$u1^2,1,0.5)-1)*sample(c(runif(meta_param$q1*meta_param$u1-meta_param$u1^2,0.38,0.42),runif(meta_param$q1*meta_param$u1-meta_param$u1^2,0.48,0.52),runif(meta_param$q1*meta_param$u1-meta_param$u1^2,0.98,1.02)),meta_param$q1*meta_param$u1-meta_param$u1^2)

## A2: q2*m: the factor loading matrix corresponds to the joint latent factor for X2
Model_Params$A2 = matrix(rbinom(meta_param$q2*meta_param$m,1,0.5)*(2*rbinom(meta_param$q2*meta_param$m,1,0.5)-1)*sample(c(runif(meta_param$q2*meta_param$m,0.38,0.42),runif(meta_param$q2*meta_param$m,0.48,0.52),runif(meta_param$q2*meta_param$m,0.98,1.02)),meta_param$q2*meta_param$m), nrow = meta_param$q2)

## B2: q2*u2: the factor loading matrix corresponds to the separate latent factor for X2
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$B2 = matrix(0, nrow = meta_param$q2, ncol = meta_param$u2)
diag(Model_Params$B2[(1:meta_param$u2),(1:meta_param$u2)]) = 1
Model_Params$B2[(1:meta_param$u2),(1:meta_param$u2)][lower.tri(Model_Params$B2[(1:meta_param$u2),(1:meta_param$u2)], diag = F)] = rbinom(meta_param$u2*(meta_param$u2-1)/2,1,0.5)*(2*rbinom(meta_param$u2*(meta_param$u2-1)/2,1,0.5)-1)*sample(c(runif(meta_param$u2*(meta_param$u2-1)/2,0.38,0.42),runif(meta_param$u2*(meta_param$u2-1)/2,0.48,0.52),runif(meta_param$u2*(meta_param$u2-1)/2,0.98,1.02)),meta_param$u2*(meta_param$u2-1)/2)
Model_Params$B2[(meta_param$u2+1):meta_param$q2,] = rbinom(meta_param$q2*meta_param$u2-meta_param$u2^2,1,0.5)*(2*rbinom(meta_param$q2*meta_param$u2-meta_param$u2^2,1,0.5)-1)*sample(c(runif(meta_param$q2*meta_param$u2-meta_param$u2^2,0.38,0.42),runif(meta_param$q2*meta_param$u2-meta_param$u2^2,0.48,0.52),runif(meta_param$q2*meta_param$u2-meta_param$u2^2,0.98,1.02)),meta_param$q2*meta_param$u2-meta_param$u2^2)

############################################  
###########  Y part parameters Gamma/Delta###########  
############################################  

## Gamma1: p1*j: the factor loading matrix corresponds to the joint latent factor for X1
Model_Params$Gamma1 = matrix(0, nrow = meta_param$p1, ncol = meta_param$j)
diag(Model_Params$Gamma1[(1:meta_param$j),(1:meta_param$j)]) = 1
Model_Params$Gamma1[(1:meta_param$j),(1:meta_param$j)][lower.tri(Model_Params$Gamma1[(1:meta_param$j),(1:meta_param$j)], diag = F)] = rbinom(meta_param$j*(meta_param$j-1)/2,1,0.5)*(2*rbinom(meta_param$j*(meta_param$j-1)/2,1,0.5)-1)*sample(c(runif(meta_param$j*(meta_param$j-1)/2,0.38,0.42),runif(meta_param$j*(meta_param$j-1)/2,0.48,0.52),runif(meta_param$j*(meta_param$j-1)/2,0.98,1.02)),meta_param$j*(meta_param$j-1)/2)
Model_Params$Gamma1[(meta_param$j+1):meta_param$p1,] = rbinom(meta_param$p1*meta_param$j-meta_param$j^2,1,0.5)*(2*rbinom(meta_param$p1*meta_param$j-meta_param$j^2,1,0.5)-1)*sample(c(runif(meta_param$p1*meta_param$j-meta_param$j^2,0.38,0.42),runif(meta_param$p1*meta_param$j-meta_param$j^2,0.48,0.52),runif(meta_param$p1*meta_param$j-meta_param$j^2,0.98,1.02)),meta_param$p1*meta_param$j-meta_param$j^2)

## Delta1: p1*z1: the factor loading matrix corresponds to the separate latent factor for Y1
##   Delta1 is in the null space of Gamma1   (need to find a note for this??)
Model_Params$Delta1 = matrix(0, nrow = meta_param$p1, ncol = meta_param$z1)
diag(Model_Params$Delta1[(1:meta_param$z1),(1:meta_param$z1)]) = 1
Model_Params$Delta1[(1:meta_param$z1),(1:meta_param$z1)][lower.tri(Model_Params$Delta1[(1:meta_param$z1),(1:meta_param$z1)], diag = F)] = rbinom(meta_param$z1*(meta_param$z1-1)/2,1,0.5)*(2*rbinom(meta_param$z1*(meta_param$z1-1)/2,1,0.5)-1)*sample(c(runif(meta_param$z1*(meta_param$z1-1)/2,0.38,0.42),runif(meta_param$z1*(meta_param$z1-1)/2,0.48,0.52),runif(meta_param$z1*(meta_param$z1-1)/2,0.98,1.02)),meta_param$z1*(meta_param$z1-1)/2)
Model_Params$Delta1[(meta_param$z1+1):meta_param$p1,] = rbinom(meta_param$p1*meta_param$z1-meta_param$z1^2,1,0.5)*(2*rbinom(meta_param$p1*meta_param$z1-meta_param$z1^2,1,0.5)-1)*sample(c(runif(meta_param$p1*meta_param$z1-meta_param$z1^2,0.38,0.42),runif(meta_param$p1*meta_param$z1-meta_param$z1^2,0.48,0.52),runif(meta_param$p1*meta_param$z1-meta_param$z1^2,0.98,1.02)),meta_param$p1*meta_param$z1-meta_param$z1^2)

## Gamma1: p2*j: the factor loading matrix corresponds to the joint latent factor for Y1
Model_Params$Gamma2 = matrix(rbinom(meta_param$p2*meta_param$j,1,0.5)*(2*rbinom(meta_param$p2*meta_param$j,1,0.5)-1)*sample(c(runif(meta_param$p2*meta_param$j,0.38,0.42),runif(meta_param$p2*meta_param$j,0.48,0.52),runif(meta_param$p2*meta_param$j,0.98,1.02)),meta_param$p2*meta_param$j), nrow = meta_param$p2)

## Delta2: p2*z2: the factor loading matrix corresponds to the separate latent factor for Y1
##   Delta2 is in the null space of Gamma2  (need to find a note for this??)
Model_Params$Delta2 = matrix(0, nrow = meta_param$p2, ncol = meta_param$z2)
diag(Model_Params$Delta2[(1:meta_param$z2),(1:meta_param$z2)]) = 1
Model_Params$Delta2[(1:meta_param$z2),(1:meta_param$z2)][lower.tri(Model_Params$Delta2[(1:meta_param$z2),(1:meta_param$z2)], diag = F)] = rbinom(meta_param$z2*(meta_param$z2-1)/2,1,0.5)*(2*rbinom(meta_param$z2*(meta_param$z2-1)/2,1,0.5)-1)*sample(c(runif(meta_param$z2*(meta_param$z2-1)/2,0.38,0.42),runif(meta_param$z2*(meta_param$z2-1)/2,0.48,0.52),runif(meta_param$z2*(meta_param$z2-1)/2,0.98,1.02)),meta_param$z2*(meta_param$z2-1)/2)
Model_Params$Delta2[(meta_param$z2+1):meta_param$p2,] = rbinom(meta_param$p2*meta_param$z2-meta_param$z2^2,1,0.5)*(2*rbinom(meta_param$p2*meta_param$z2-meta_param$z2^2,1,0.5)-1)*sample(c(runif(meta_param$p2*meta_param$z2-meta_param$z2^2,0.38,0.42),runif(meta_param$p2*meta_param$z2-meta_param$z2^2,0.48,0.52),runif(meta_param$p2*meta_param$z2-meta_param$z2^2,0.98,1.02)),meta_param$p2*meta_param$z2-meta_param$z2^2)

#############################################################################  
###########  Y~X part (associations between latent factors) /Pi/Theta ###########  
#############################################################################  

## Theta: j*m:  the association between the joint latent factors (m in total) for X and the joint latent factors (j in total) for Y
Model_Params$Theta <- matrix(runif(meta_param$j*meta_param$m, -1, 1),nrow = meta_param$j)
## Psi: j*u:  the association between the separate latent factors (u = u1 + u2 in total) for X and the joint latent factors (j in total) for Y
## Psi is Lambda in the model (overleaf) 
Model_Params$Psi <- matrix(runif(meta_param$j*meta_param$u, -1, 1),nrow = meta_param$j)

## z1*m the association between the joint latent factors (m in total) for X and the joint latent factors (z1 in total) for Y1
Model_Params$Pi1 <- matrix(0, nrow = meta_param$z1, ncol = meta_param$m)
## z2*m the association between the joint latent factors (m in total) for X and the joint latent factors (z2 in total) for Y2
Model_Params$Pi2 <- matrix(runif(meta_param$z2*meta_param$m, -1, 1),nrow = meta_param$z2)
## Omega1: z1*u:  the association between the separate latent factors (u = u1+u2 in total) for X (X1+X2) and the joint latent factors (z1 in total) for Y1
Model_Params$Omega1 <- matrix(0, nrow = meta_param$z1, ncol = meta_param$u)
## Omega2: z2*u:  the association between the separate latent factors (u = u1+u2 in total) for X (X1+X2) and the joint latent factors (z2 in total) for Y2
Model_Params$Omega2 <- matrix(runif(meta_param$z2*meta_param$u, -1, 1),nrow = meta_param$z2)


#############################################################################  
## residual part and covariance matrix of latent factors
##   Phi.../
#############################################################################  

## covariance of shared latent factors (m in toal) for X part 
Model_Params$PhiR <- diag(runif(meta_param$m, 0.8, 1.2))

## covariance of separate latent factors (u1 in toal) for X1 part 
Model_Params$PhiS1 <- diag(runif(meta_param$u1, 0.8, 1.2))
## covariance of separate latent factors (u2 in toal) for X2 part 
Model_Params$PhiS2 <- diag(runif(meta_param$u2, 0.8, 1.2))

## covariance os S1 vs S 
Model_Params$PhiS1S <- cbind(Model_Params$PhiS1, matrix(0, nrow = meta_param$u1, ncol = meta_param$u2)) ## u1*(u1+u2)
## covariance os S2 vs S 
Model_Params$PhiS2S <- cbind(matrix(0, nrow = meta_param$u2, ncol = meta_param$u1), Model_Params$PhiS2) ## u2*(u1+u2)
## What is this: PhiS: the diagonal matrix of covariance of S1 and S2
Model_Params$PhiS <- diag(c(diag(Model_Params$PhiS1), diag(Model_Params$PhiS2))) ##  (u1+u2)*(u1+u2)

## covariance of the residuals of X1 part  
Model_Params$Phi11 <- diag(runif(meta_param$q1, 0.8, 1.2))
## covariance of the residuals of X2 part  
Model_Params$Phi12 <- diag(runif(meta_param$q2, 0.8, 1.2))

## covariance of the residuals of Y1 part  
Model_Params$Phi21 <- diag(runif(meta_param$p1, 0.8, 1.2))
## covariance of the residuals of Y2 part  
Model_Params$Phi22 <- diag(runif(meta_param$p2, 0.8, 1.2))

## covariance of the residuals of V (joint latent factor of Y) ~ R  (joint latent factor of X) +S  (separate latent factor of X)     
Model_Params$Phi3 <- diag(runif(meta_param$j, 0.8, 1.2))
## covariance of the residuals of W1 (separate latent factor of Y1) ~ R  (joint latent factor of X) +S  (separate latent factor of X)     
Model_Params$Phi41 <- diag(runif(meta_param$z1, 0.8, 1.2))
## covariance of the residuals of W2 (separate latent factor of Y2) ~ R  (joint latent factor of X) +S  (separate latent factor of X)     
Model_Params$Phi42 <- diag(runif(meta_param$z2, 0.8, 1.2))

## put all in a list
# Model_Params <- list()
# Model_Params$PhiR     <- PhiR
# Model_Params$PhiS1    <- PhiS1  
# Model_Params$PhiS2    <- PhiS2
# Model_Params$Phi11    <- Phi11  
# Model_Params$Phi12    <- Phi12
# Model_Params$A1       <- A1     
# Model_Params$A2       <- A2     
# Model_Params$B1       <- B1     
# Model_Params$B2       <- B2     
# 
# Model_Params$Phi3      <- Phi3
# Model_Params$Phi41     <- Phi41  
# Model_Params$Phi42     <- Phi42  
# Model_Params$Phi21     <- Phi21  
# Model_Params$Phi22     <- Phi22
# Model_Params$Theta     <- Theta  
# Model_Params$Psi       <- Psi    
# Model_Params$Pi1       <- Pi1    
# Model_Params$Pi2       <- Pi2
# Model_Params$Omega1    <- Omega1
# Model_Params$Omega2    <- Omega2 
# Model_Params$Gamma1    <- Gamma1 
# Model_Params$Gamma2    <- Gamma2 
# Model_Params$Delta1    <- Delta1 
# Model_Params$Delta2    <- Delta2 


## end of code 
