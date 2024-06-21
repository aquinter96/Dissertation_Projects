## Model_Matrix_Param.R

##### To generate data (X, Y) given model parameters #####

Model_Params <- list()

############################################  
###########  X part parameters B0/B ###########  
############################################  

## B0: q*1: the intercept vector corresponding to the X data
Model_Params$B0 = rnorm(meta_param$q)

## B: q*s: the factor loading matrix corresponding to the latent factor for X
## randomly generate individual elements from uniform distribution, randomly assign them as positive or negative,
## and lastly randomly assign some elements to be 0 in order to impose sparsity
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$B = matrix(0, nrow = meta_param$q, ncol = meta_param$s)
diag(Model_Params$B[(1:meta_param$s),(1:meta_param$s)]) = 1
Model_Params$B[(1:meta_param$s),(1:meta_param$s)][lower.tri(Model_Params$B[(1:meta_param$s),(1:meta_param$s)], diag = F)] = rbinom(meta_param$s,1,0.5)*(2*rbinom(meta_param$s,1,0.5)-1)*runif(meta_param$s, 0.35, 1)
Model_Params$B[(meta_param$s+1):meta_param$q,] = rbinom((meta_param$q*meta_param$s - meta_param$s^2),1,0.5)*(2*rbinom((meta_param$q*meta_param$s - meta_param$s^2),1,0.5)-1)*runif((meta_param$q*meta_param$s - meta_param$s^2), 0.35, 1)

############################################  
###########  Y part parameters A0/A ###########  
############################################  

## A0: p*1: the intercept vector corresponding to the Y data
Model_Params$A0 = rnorm(meta_param$p)

## A: p*m: the factor loading matrix corresponding to the latent factor for Y
## randomly generate individual elements from uniform distribution, randomly assign them as positive or negative,
## and lastly randomly assign some elements to be 0 in order to impose sparsity
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$A <- matrix(0, nrow = meta_param$p, ncol = meta_param$m)
diag(Model_Params$A[(1:meta_param$m),(1:meta_param$m)]) <- 1
Model_Params$A[(1:meta_param$m),(1:meta_param$m)][lower.tri(Model_Params$A[(1:meta_param$m),(1:meta_param$m)], diag = F)] <- rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)*(2*rbinom(meta_param$m*(meta_param$m-1)/2,1,0.5)-1)*runif(meta_param$m*(meta_param$m-1)/2, 0.35, 1)
Model_Params$A[(meta_param$m+1):meta_param$p,] <- rbinom((meta_param$p*meta_param$m - meta_param$m^2),1,0.5)*(2*rbinom((meta_param$p*meta_param$m - meta_param$m^2),1,0.5)-1)*runif((meta_param$p*meta_param$m - meta_param$m^2), 0.35, 1)

#############################################################################  
###########  Y~X part (associations between latent factors) Gamma ###########  
#############################################################################  

## Gamma: m*s:  the association between the latent factors (s in total) for X and the latent factors (m in total) for Y
Model_Params$Gamma <- matrix(rnorm(meta_param$s*meta_param$m),nrow = meta_param$m)

#############################################################################  
## residual part and covariance matrix of latent factors
##   Phi1/Phi2/Phi3
#############################################################################  

## covariance of the residuals of X part  
Model_Params$Phi1 <- diag(abs(rnorm(meta_param$q)))

## covariance of the residuals of Y part  
Model_Params$Phi2 <- diag(abs(rnorm(meta_param$p)))

## covariance of the residuals of Z (latent factor of Y) ~ W  (latent factor of X)    
Model_Params$Phi3 <- diag(abs(rnorm(meta_param$m)))

## end of code 
