
library(MASS)
library(pracma)
library(ggplot2)

setwd("~/Paper2_Parallel")

set.seed(2435)

### To initialize parameters

code_dir <- paste0("~/Paper2_Parallel/")
source(file = paste0(code_dir, "Model_Meta_Param.R") )
#source(file = paste0(code_dir, "Model_Matrix_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param_Sparse.R"))

XMSE <- data.frame(MSE = double(), sampsize = double(), Matrix = character())
XeMSE <- data.frame(MSE = double(), sampsize = double(), Matrix = character())
Xtnr <- data.frame(tnr = double(), sampsize = double(), Matrix = character())
Xtpr <- data.frame(tpr = double(), sampsize = double(), Matrix = character())
YMSE <- data.frame(MSE = double(), sampsize = double(), Matrix = character())
YeMSE <- data.frame(MSE = double(), sampsize = double(), Matrix = character())
YfMSE <- data.frame(MSE = double(), sampsize = double(), Matrix = character())
Ytnr <- data.frame(tnr = double(), sampsize = double(), Matrix = character())
Ytpr <- data.frame(tpr = double(), sampsize = double(), Matrix = character())
Ytime <- data.frame(time = double(), sampsize = double())

Ests <- NULL

iter1 <- 1
iter2 <- 1

for(i in 1:50){
for(j in c(200, 400, 600, 800)){
  if(file.exists(paste0("NMflatrep", i, "n", j, ".rds"))){
  Ests <- readRDS(paste0("NMflatrep", i, "n", j, ".rds"))
  XMSE[iter1,] <- c(1/(meta_param$q1*meta_param$m-meta_param$m*(meta_param$m+1)/2)*sum((Ests$X_estimates$X_final_pars$A1 - Model_Params$A1)^2), j, "A1")
  Xtnr[iter1,] <- c(0, j, "A1")
  Xtpr[iter1,] <- c(0, j, "A1")
  for(a in 1:meta_param$q1){
    for(b in 1:meta_param$m){
      if(a > b){
        if((Ests$X_estimates$X_final_pars$A1[a,b] == 0) & (Model_Params$A1[a,b] == 0)){
          Xtnr[iter1,1] <- as.numeric(Xtnr[iter1,1]) + 1/(sum(Model_Params$A1 == 0) - meta_param$m*(meta_param$m-1)/2)
        }
        if((Ests$X_estimates$X_final_pars$A1[a,b] != 0) & (Model_Params$A1[a,b] != 0)){
          Xtpr[iter1,1] <- as.numeric(Xtpr[iter1,1]) + 1/(sum(Model_Params$A1 != 0) - meta_param$m)
        }
      }
    }
  }
  iter1 <- iter1 + 1
  XMSE[iter1,] <- c(1/(meta_param$q1*meta_param$u1-meta_param$u1*(meta_param$u1+1)/2)*sum((Ests$X_estimates$X_final_pars$B1 - Model_Params$B1)^2), j, "B1")
  Xtnr[iter1,] <- c(0, j, "B1")
  Xtpr[iter1,] <- c(0, j, "B1")
  for(a in 1:meta_param$q1){
    for(b in 1:meta_param$u1){
      if(a > b){
        if((Ests$X_estimates$X_final_pars$B1[a,b] == 0) & (Model_Params$B1[a,b] == 0)){
          Xtnr[iter1,1] <- as.numeric(Xtnr[iter1,1]) + 1/(sum(Model_Params$B1 == 0) - meta_param$u1*(meta_param$u1-1)/2)
        }
        if((Ests$X_estimates$X_final_pars$B1[a,b] != 0) & (Model_Params$B1[a,b] != 0)){
          Xtpr[iter1,1] <- as.numeric(Xtpr[iter1,1]) + 1/(sum(Model_Params$B1 != 0) - meta_param$u1)
        }
      }
    }
  }
  iter1 <- iter1 + 1
  XMSE[iter1,] <- c(1/(meta_param$q2*meta_param$m-meta_param$m*(meta_param$m+1)/2)*sum((Ests$X_estimates$X_final_pars$A2 - Model_Params$A2)^2), j, "A2")
  Xtnr[iter1,] <- c(0, j, "A2")
  Xtpr[iter1,] <- c(0, j, "A2")
  for(a in 1:meta_param$q2){
    for(b in 1:meta_param$m){
      if((Ests$X_estimates$X_final_pars$A2[a,b] == 0) & (Model_Params$A2[a,b] == 0)){
        Xtnr[iter1,1] <- as.numeric(Xtnr[iter1,1]) + 1/(sum(Model_Params$A2 == 0))
      }
      if((Ests$X_estimates$X_final_pars$A2[a,b] != 0) & (Model_Params$A2[a,b] != 0)){
        Xtpr[iter1,1] <- as.numeric(Xtpr[iter1,1]) + 1/(sum(Model_Params$A2 != 0))
      }
    }
  }
  iter1 <- iter1 + 1
  XMSE[iter1,] <- c(1/(meta_param$q2*meta_param$u2-meta_param$u2*(meta_param$u2+1)/2)*sum((Ests$X_estimates$X_final_pars$B2 - Model_Params$B2)^2), j, "B2")
  Xtnr[iter1,] <- c(0, j, "B2")
  Xtpr[iter1,] <- c(0, j, "B2")
  for(a in 1:meta_param$q2){
    for(b in 1:meta_param$u2){
      if(a > b){
        if((Ests$X_estimates$X_final_pars$B2[a,b] == 0) & (Model_Params$B2[a,b] == 0)){
          Xtnr[iter1,1] <- as.numeric(Xtnr[iter1,1]) + 1/(sum(Model_Params$B2 == 0) - meta_param$u2*(meta_param$u2-1)/2)
        }
        if((Ests$X_estimates$X_final_pars$B2[a,b] != 0) & (Model_Params$B2[a,b] != 0)){
          Xtpr[iter1,1] <- as.numeric(Xtpr[iter1,1]) + 1/(sum(Model_Params$B2 != 0) - meta_param$u2)
        }
      }
    }
  }
  iter1 <- iter1 + 1
  XeMSE[iter2,] <- c(1/meta_param$m*sum((Ests$X_estimates$X_final_pars$PhiR - Model_Params$PhiR)^2), j, "PhiR")
  iter2 <- iter2 + 1
  XeMSE[iter2,] <- c(1/meta_param$u1*sum((Ests$X_estimates$X_final_pars$PhiS1 - Model_Params$PhiS1)^2), j, "PhiS1")
  iter2 <- iter2 + 1
  XeMSE[iter2,] <- c(1/meta_param$u2*sum((Ests$X_estimates$X_final_pars$PhiS2 - Model_Params$PhiS2)^2), j, "PhiS2")
  iter2 <- iter2 + 1
  XeMSE[iter2,] <- c(1/meta_param$q1*sum((Ests$X_estimates$X_final_pars$Phi11 - Model_Params$Phi11)^2), j, "Phi11")
  iter2 <- iter2 + 1
  XeMSE[iter2,] <- c(1/meta_param$q2*sum((Ests$X_estimates$X_final_pars$Phi12 - Model_Params$Phi12)^2), j, "Phi12")
  iter2 <- iter2 + 1
}
}
}

tapply(as.numeric(Xtnr$tnr), list(Xtnr$Matrix, Xtnr$sampsize), mean)
tapply(as.numeric(Xtpr$tpr), list(Xtpr$Matrix, Xtnr$sampsize), mean)




iter1 <- 1
iter2 <- 1
iter3 <- 1
k <- 0

  for(i in 1:50){
  for(j in c(200, 400, 600, 800)){
    if(file.exists(paste0("NMflatrep", i, "n", j, ".rds"))){
    k <- k + 1
    Ests <- readRDS(paste0("NMflatrep", i, "n", j, ".rds"))
    YMSE[iter1,] <- c(1/(meta_param$p1*meta_param$j-meta_param$j*(meta_param$j+1)/2)*sum((Ests$Y_estimates$Y_final_pars$Gamma1 - Model_Params$Gamma1)^2), j, "Gamma1")
    Ytnr[iter1,] <- c(0, j, "Gamma1")
    Ytpr[iter1,] <- c(0, j, "Gamma1")
    for(a in 1:meta_param$p1){
      for(b in 1:meta_param$j){
        if(a > b){
          if((Ests$Y_estimates$Y_final_pars$Gamma1[a,b] == 0) & (Model_Params$Gamma1[a,b] == 0)){
            Ytnr[iter1,1] <- as.numeric(Ytnr[iter1,1]) + 1/(sum(Model_Params$Gamma1 == 0) - meta_param$j*(meta_param$j-1)/2)
          }
          if((Ests$Y_estimates$Y_final_pars$Gamma1[a,b] != 0) & (Model_Params$Gamma1[a,b] != 0)){
            Ytpr[iter1,1] <- as.numeric(Ytpr[iter1,1]) + 1/(sum(Model_Params$Gamma1 != 0) - meta_param$j)
          }
        }
      }
    }
    iter1 <- iter1 + 1
    YMSE[iter1,] <- c(1/(meta_param$p1*meta_param$z1-meta_param$z1*(meta_param$z1+1)/2)*sum((Ests$Y_estimates$Y_final_pars$Delta1 - Model_Params$Delta1)^2), j, "Delta1")
    Ytnr[iter1,] <- c(0, j, "Delta1")
    Ytpr[iter1,] <- c(0, j, "Delta1")
    for(a in 1:meta_param$p1){
      for(b in 1:meta_param$z1){
        if(a > b){
          if((Ests$Y_estimates$Y_final_pars$Delta1[a,b] == 0) & (Model_Params$Delta1[a,b] == 0)){
            Ytnr[iter1,1] <- as.numeric(Ytnr[iter1,1]) + 1/(sum(Model_Params$Delta1 == 0) - meta_param$z1*(meta_param$z1-1)/2)
          }
          if((Ests$Y_estimates$Y_final_pars$Delta1[a,b] != 0) & (Model_Params$Delta1[a,b] != 0)){
            Ytpr[iter1,1] <- as.numeric(Ytpr[iter1,1]) + 1/(sum(Model_Params$Delta1 != 0) - meta_param$z1)
          }
        }
      }
    }
    iter1 <- iter1 + 1
    YMSE[iter1,] <- c(1/(meta_param$p2*meta_param$j-meta_param$j*(meta_param$j+1)/2)*sum((Ests$Y_estimates$Y_final_pars$Gamma2 - Model_Params$Gamma2)^2), j, "Gamma2")
    Ytnr[iter1,] <- c(0, j, "Gamma2")
    Ytpr[iter1,] <- c(0, j, "Gamma2")
    for(a in 1:meta_param$p2){
      for(b in 1:meta_param$j){
        if((Ests$Y_estimates$Y_final_pars$Gamma2[a,b] == 0) & (Model_Params$Gamma2[a,b] == 0)){
          Ytnr[iter1,1] <- as.numeric(Ytnr[iter1,1]) + 1/(sum(Model_Params$Gamma2 == 0))
        }
        if((Ests$Y_estimates$Y_final_pars$Gamma2[a,b] != 0) & (Model_Params$Gamma2[a,b] != 0)){
          Ytpr[iter1,1] <- as.numeric(Ytpr[iter1,1]) + 1/(sum(Model_Params$Gamma2 != 0))
        }
      }
    }
    iter1 <- iter1 + 1
    YMSE[iter1,] <- c(1/(meta_param$p2*meta_param$z2-meta_param$z2*(meta_param$z2+1)/2)*sum((Ests$Y_estimates$Y_final_pars$Delta2 - Model_Params$Delta2)^2), j, "Delta2")
    Ytnr[iter1,] <- c(0, j, "Delta2")
    Ytpr[iter1,] <- c(0, j, "Delta2")
    for(a in 1:meta_param$p2){
      for(b in 1:meta_param$z2){
        if(a > b){
          if((Ests$Y_estimates$Y_final_pars$Delta2[a,b] == 0) & (Model_Params$Delta2[a,b] == 0)){
            Ytnr[iter1,1] <- as.numeric(Ytnr[iter1,1]) + 1/(sum(Model_Params$Delta2 == 0) - meta_param$z2*(meta_param$z2-1)/2)
          }
          if((Ests$Y_estimates$Y_final_pars$Delta2[a,b] != 0) & (Model_Params$Delta2[a,b] != 0)){
            Ytpr[iter1,1] <- as.numeric(Ytpr[iter1,1]) + 1/(sum(Model_Params$Delta2 != 0) - meta_param$z2)
          }
        }
      }
    }
    Ytime[k,] <- c(Ests$tot.time[[3]], j)
    iter1 <- iter1 + 1
    YeMSE[iter2,] <- c(1/meta_param$j*sum((Ests$Y_estimates$Y_final_pars$Phi3 - Model_Params$Phi3)^2), j, "Phi3")
    iter2 <- iter2 + 1
    YeMSE[iter2,] <- c(1/meta_param$z1*sum((Ests$Y_estimates$Y_final_pars$Phi41 - Model_Params$Phi41)^2), j, "Phi41")
    iter2 <- iter2 + 1
    YeMSE[iter2,] <- c(1/meta_param$z2*sum((Ests$Y_estimates$Y_final_pars$Phi42 - Model_Params$Phi42)^2), j, "Phi42")
    iter2 <- iter2 + 1
    YeMSE[iter2,] <- c(1/meta_param$p1*sum((Ests$Y_estimates$Y_final_pars$Phi21 - Model_Params$Phi21)^2), j, "Phi21")
    iter2 <- iter2 + 1
    YeMSE[iter2,] <- c(1/meta_param$p2*sum((Ests$Y_estimates$Y_final_pars$Phi22 - Model_Params$Phi22)^2), j, "Phi22")
    iter2 <- iter2 + 1
    YfMSE[iter3,] <- c(1/10*sum((Ests$Y_estimates$Y_final_pars$Theta - Model_Params$Theta)^2), j, "Theta")
    iter3 <- iter3 + 1
    YfMSE[iter3,] <- c(1/25*sum((Ests$Y_estimates$Y_final_pars$Psi - Model_Params$Psi)^2), j, "Psi")
    iter3 <- iter3 + 1
    # YfMSE[iter3,] <- c(1/6*sum((Ests$Y_estimates$Y_final_pars$Pi1 - Model_Params$Pi1)^2), j, "Pi1")
    # iter3 <- iter3 + 1
    YfMSE[iter3,] <- c(1/6*sum((Ests$Y_estimates$Y_final_pars$Pi2 - Model_Params$Pi2)^2), j, "Pi2")
    iter3 <- iter3 + 1
    # YfMSE[iter3,] <- c(1/15*sum((Ests$Y_estimates$Y_final_pars$Omega1 - Model_Params$Omega1)^2), j, "Omega1")
    # iter3 <- iter3 + 1
    YfMSE[iter3,] <- c(1/15*sum((Ests$Y_estimates$Y_final_pars$Omega2 - Model_Params$Omega2)^2), j, "Omega2")
    iter3 <- iter3 + 1
}
}
}

tapply(as.numeric(Ytnr$tnr), list(Ytnr$Matrix, Ytnr$sampsize), mean)
tapply(as.numeric(Ytpr$tpr), list(Ytpr$Matrix, Ytpr$sampsize), mean)
tapply(as.numeric(XMSE$MSE), list(XMSE$Matrix, XMSE$sampsize), mean)
tapply(as.numeric(YMSE$MSE), list(YMSE$Matrix, YMSE$sampsize), mean)
tapply(as.numeric(Ytime$time), Ytime$sampsize, mean)/60
