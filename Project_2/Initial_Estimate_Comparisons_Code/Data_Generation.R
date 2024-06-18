## Data_Generation.R

##### to generate data (X,Y), given model parameters ... 

Data_Generation <- function(n = 1000, Model_Params){

  Data <- list()
  
  MyPar <- Model_Params
  
  m <- ncol(MyPar$PhiR)
  u1 <- ncol(MyPar$PhiS1)
  u2 <-  ncol(MyPar$PhiS2)   
  q1 <-  ncol(MyPar$Phi11)
  q2 <-  ncol(MyPar$Phi12)
  
  p1 <-  ncol(MyPar$Phi21)
  p2 <-   ncol(MyPar$Phi22)
  j <- ncol(MyPar$Phi3)
  z1 <- ncol(MyPar$Phi41)
  z2 <- ncol(MyPar$Phi42)  
  
  ## to simulate latent factor and then X/Y
  ## latent factors for X 
  Data$R <- t(MASS::mvrnorm(n, rep(0, m), MyPar$PhiR))    ##  shared latent factors for X 
  Data$S1 <- t(MASS::mvrnorm(n, rep(0, u1), MyPar$PhiS1)) ##  Separate latent factors for X1
  Data$S2 <- t(MASS::mvrnorm(n, rep(0, u2), MyPar$PhiS2)) ##  Separate latent factors for X2

  Data$X1 <- MyPar$A1%*%Data$R + MyPar$B1%*%Data$S1 + t(MASS::mvrnorm(n, rep(0, q1), MyPar$Phi11)) ## X1
  Data$X2 <- MyPar$A2%*%Data$R + MyPar$B2%*%Data$S2 + t(MASS::mvrnorm(n, rep(0, q2), MyPar$Phi12)) ## X2 

  ## latent factors for Y 
  Data$V <- MyPar$Theta%*%Data$R + MyPar$Psi%*%rbind(Data$S1,Data$S2) + t(MASS::mvrnorm(n, rep(0, j), MyPar$Phi3)) ##  shared latent factors for Y 
  Data$W1 <- MyPar$Pi1%*%Data$R + MyPar$Omega1%*%rbind(Data$S1,Data$S2) + t(MASS::mvrnorm(n, rep(0, z1), MyPar$Phi41)) ##  Separate latent factors for Y1
  Data$W2 <- MyPar$Pi2%*%Data$R + MyPar$Omega2%*%rbind(Data$S1,Data$S2) + t(MASS::mvrnorm(n, rep(0, z2), MyPar$Phi42)) ##  Separate latent factors for Y2

  Data$Y1 <- MyPar$Gamma1%*%Data$V + MyPar$Delta1%*%Data$W1 + t(MASS::mvrnorm(n, rep(0, p1), MyPar$Phi21))
  Data$Y2 <- MyPar$Gamma2%*%Data$V + MyPar$Delta2%*%Data$W2 + t(MASS::mvrnorm(n, rep(0, p2), MyPar$Phi22))

  return(Data)  

}

## end of code 
