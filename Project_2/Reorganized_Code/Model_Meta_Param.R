## Model_Meta_Param.R

##### parameters on dimensions of X, Y, and latent factors  

meta_param <- list()

meta_param$m = 2  ## number of joint/shared latent factors for X part 
meta_param$q1 = 20  ## number of variables in X1 (block 1 of X) 
meta_param$u1 = 2  ## number of separate latent factors for X1 part
meta_param$q2 = 30  ## number of variables in X2 (block 2 of X)
meta_param$u2 = 3  ## number of separate latent factors for X2 part
meta_param$u = meta_param$u1 + meta_param$u2

meta_param$j = 5  ## number of joint/shared latent factors for Y part  
meta_param$p1 = 30   ## number of variables in Y1 (block 1 of Y)
meta_param$z1 = 3   ## number of variables in Y1 (block 1 of Y) 
meta_param$p2 = 20    ## number of variables in Y2 (block 2 of Y)
meta_param$z2 = 3   ## number of variables in Y2 (block 2 of Y) 

## end of code 
