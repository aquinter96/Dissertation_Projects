## Model1_Matrix_Param.R

##### To generate data (X, Y) given model parameters #####

Model_Params <- list()

############################################  
###########  X part parameters B0/B ###########  
############################################  

## B0: q*1: the intercept vector corresponding to the X data
Model_Params$B0 = runif(meta_param$q, -1, 1)

## B: q*s: the factor loading matrix corresponding to the latent factor for X
## randomly generate individual elements from uniform distribution, randomly assign them as positive or negative,
## and lastly randomly assign some elements to be 0 in order to impose sparsity
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$B = matrix(c(1, -0.65, 0.65, rep(0,7), rep(0.70, 10), 0, 1, 0.70, rep(-0.85, 7), rep(0, 10), 0, 0, 1, rep(0, 7), rep(-0.85, 10)), nrow = meta_param$q, ncol = meta_param$s)

############################################  
###########  Y part parameters A0/A ###########  
############################################  

## A0: p*1: the intercept vector corresponding to the Y data
Model_Params$A0 = runif(meta_param$p, -1, 1)

## A: p*m: the factor loading matrix corresponding to the latent factor for Y
## randomly generate individual elements from uniform distribution, randomly assign them as positive or negative,
## and lastly randomly assign some elements to be 0 in order to impose sparsity
##  lower-triangular matrix with diagonal elements being 1 (identifiability for factor analysis)
Model_Params$A <- matrix(c(1, 0.45, -0.80, rep(0.80, 7), rep(0, 10), 0, 1, 0.70, rep(0, 7), rep(-0.95, 10)), nrow = meta_param$p, ncol = meta_param$m)


#############################################################################  
###########  Y~X part (associations between latent factors) Gamma ###########  
#############################################################################  

## Gamma: m*s:  the association between the latent factors (s in total) for X and the latent factors (m in total) for Y
Model_Params$Gamma <- matrix(c(0.60, -0.75, -0.70, 0.65, 0.80, -0.55), nrow = meta_param$m, ncol = meta_param$s)

#############################################################################  
## residual part and covariance matrix of latent factors
##   Phi1/Phi2/Phi3
#############################################################################  

## covariance of the residuals of X part  
Model_Params$Phi1 <- diag(runif(meta_param$q, 0.9, 1.1))

## covariance of the residuals of Y part  
Model_Params$Phi2 <- diag(runif(meta_param$p, 0.9, 1.1))

## covariance of the residuals of Z (latent factor of Y) ~ W  (latent factor of X)    
Model_Params$Phi3 <- diag(runif(meta_param$m, 0.9, 1.1))

## end of code 
