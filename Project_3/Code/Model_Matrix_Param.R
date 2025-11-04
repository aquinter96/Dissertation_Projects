## Model_Matrix_Param.R

##### parameters on dimensions of X, Y, and latent factors  

Model_Params <- list()

############################################  
###########  X part parameters A1/A2/B1/B2 ###########  
############################################  

Model_Params$B0 = runif(meta_param$q, -1, 1)

## B2: q2*u2: the factor loading matrix corresponds to the separate latent factor for X2
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$B <- matrix(0, nrow = meta_param$q, ncol = meta_param$s)
diag(Model_Params$B[(1:meta_param$s),(1:meta_param$s)]) <- 1
Model_Params$B[(1:meta_param$s),(1:meta_param$s)][lower.tri(Model_Params$B[(1:meta_param$s),(1:meta_param$s)], diag = F)] <- rbinom(meta_param$s*(meta_param$s-1)/2,1,0.5)*(2*rbinom(meta_param$s*(meta_param$s-1)/2,1,0.5)-1)*sample(c(runif(meta_param$s*(meta_param$s-1)/2,0.38,0.42),runif(meta_param$s*(meta_param$s-1)/2,0.48,0.52),runif(meta_param$s*(meta_param$s-1)/2,0.98,1.02)),meta_param$s*(meta_param$s-1)/2)
Model_Params$B[(meta_param$s+1):meta_param$q,] <- rbinom(meta_param$q*meta_param$s-meta_param$s^2,1,0.5)*(2*rbinom(meta_param$q*meta_param$s-meta_param$s^2,1,0.5)-1)*sample(c(runif(meta_param$q*meta_param$s-meta_param$s^2,0.38,0.42),runif(meta_param$q*meta_param$s-meta_param$s^2,0.48,0.52),runif(meta_param$q*meta_param$s-meta_param$s^2,0.98,1.02)),meta_param$q*meta_param$s-meta_param$s^2)

Model_Params$A0 = runif(meta_param$p, -1, 1)

## A1: q1*m: the factor loading matrix corresponds to the joint latent factor for X1
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis), for A = (A1': A2')'
Model_Params$A <- matrix(0, nrow = meta_param$p, ncol = meta_param$m)
diag(Model_Params$A[(1:meta_param$m),(1:meta_param$m)]) <- 1
Model_Params$A[(1:meta_param$m),(1:meta_param$m)][lower.tri(Model_Params$A[(1:meta_param$m),(1:meta_param$m)], diag = F)] <- rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)*(2*rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)-1)*sample(c(runif(meta_param$m*(meta_param$m-1)/2,0.38,0.42),runif(meta_param$m*(meta_param$m-1)/2,0.48,0.52),runif(meta_param$m*(meta_param$m-1)/2,0.98,1.02)),meta_param$m*(meta_param$m-1)/2)
Model_Params$A[(meta_param$m+1):meta_param$p,] <- rbinom(meta_param$p*meta_param$m-meta_param$m^2,1,0.5)*(2*rbinom(meta_param$p*meta_param$m-meta_param$m^2,1,0.5)-1)*sample(c(runif(meta_param$p*meta_param$m-meta_param$m^2,0.38,0.42),runif(meta_param$p*meta_param$m-meta_param$m^2,0.48,0.52),runif(meta_param$p*meta_param$m-meta_param$m^2,0.98,1.02)),meta_param$p*meta_param$m-meta_param$m^2)

Model_Params$C0 = runif(meta_param$r, -1, 1)

## A1: q1*m: the factor loading matrix corresponds to the joint latent factor for X1
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis), for A = (A1': A2')'
Model_Params$C <- matrix(0, nrow = meta_param$r, ncol = meta_param$t)
diag(Model_Params$C[(1:meta_param$t),(1:meta_param$t)]) <- 1
Model_Params$C[(1:meta_param$t),(1:meta_param$t)][lower.tri(Model_Params$C[(1:meta_param$t),(1:meta_param$t)], diag = F)] <- rbinom(meta_param$t*(meta_param$t-1)/2,1,0.5)*(2*rbinom(meta_param$t*(meta_param$t-1)/2,1,0.5)-1)*sample(c(runif(meta_param$t*(meta_param$t-1)/2,0.38,0.42),runif(meta_param$t*(meta_param$t-1)/2,0.48,0.52),runif(meta_param$t*(meta_param$t-1)/2,0.98,1.02)),meta_param$t*(meta_param$t-1)/2)
Model_Params$C[(meta_param$t+1):meta_param$r,] <- rbinom(meta_param$r*meta_param$t-meta_param$t^2,1,0.5)*(2*rbinom(meta_param$r*meta_param$t-meta_param$t^2,1,0.5)-1)*sample(c(runif(meta_param$r*meta_param$t-meta_param$t^2,0.38,0.42),runif(meta_param$r*meta_param$t-meta_param$t^2,0.48,0.52),runif(meta_param$r*meta_param$t-meta_param$t^2,0.98,1.02)),meta_param$r*meta_param$t-meta_param$t^2)

## B2: q2*u2: the factor loading matrix corresponds to the separate latent factor for X2
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$D <- matrix(0, nrow = meta_param$r, ncol = meta_param$s)
Model_Params$D[1:meta_param$r,] <- rbinom(meta_param$r*meta_param$s,1,0.5)*(2*rbinom(meta_param$r*meta_param$s,1,0.5)-1)*sample(c(runif(meta_param$r*meta_param$s,0.38,0.42),runif(meta_param$r*meta_param$s,0.48,0.52),runif(meta_param$r*meta_param$s,0.98,1.02)),meta_param$r*meta_param$s)


############################################  
###########  Y part parameters Gamma/Delta###########  
############################################  

## Gamma1: p2*j: the factor loading matrix corresponds to the joint latent factor for Y1
Model_Params$Gamma = matrix(runif(meta_param$m*meta_param$s, -1, 1),nrow = meta_param$m)

## Gamma1: p2*j: the factor loading matrix corresponds to the joint latent factor for Y1
Model_Params$Delta = matrix(runif(meta_param$t*meta_param$m, -1, 1),nrow = meta_param$t)

#############################################################################  
## residual part and covariance matrix of latent factors
##   Phi.../
#############################################################################  

## covariance of the residuals of X1 part  
Model_Params$Psi <- diag(runif(meta_param$s, 0.8, 1.2))

## covariance of the residuals of X1 part  
Model_Params$Phi2 <- diag(runif(meta_param$m, 0.8, 1.2))

## covariance of the residuals of X1 part  
Model_Params$Phi3 <- diag(runif(meta_param$r, 0.8, 1.2))

## covariance of the residuals of X1 part  
Model_Params$Phi4 <- diag(runif(meta_param$t, 0.8, 1.2))

## covariance of the residuals of X1 part  
Model_Params$Phi5 <- diag(runif(meta_param$p, 0.8, 1.2))

## covariance of the residuals of X1 part  
Model_Params$Phi1 <- diag(runif(meta_param$q, 0.8, 1.2))

## end of code
