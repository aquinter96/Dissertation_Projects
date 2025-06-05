## AlexProject1_CodeTesting_v01

### Random seed 

library(pracma)
library(MASS)
library(r.jive)
library(parallel)
library(nloptr)

args <- commandArgs(TRUE)
datset <- as.numeric(args[1])
sampsize <- as.numeric(args[2])
set.seed(2435)


### To initialize parameters


code_dir <- paste0("~/Paper2_Parallel/")
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param_Sparse_Cons2.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(datset)
MyData <- Data_Generation(n = sampsize, Model_Params=Model_Params)


### To get estimates for the X part ...

source(file = paste0(code_dir, "Lasso.R") )
source(file = paste0(code_dir, "X_Inits2.R") )
source(file = paste0(code_dir, "Estep_X.R") )
source(file = paste0(code_dir, "Mstep_X.R") )
source(file = paste0(code_dir, "logLik_X.R") )
source(file = paste0(code_dir, "Convergence_check.R") )
source(file = paste0(code_dir, "BIC_X.R") )
source(file = paste0(code_dir, "EMAlg_XNM.R") )
source(file = paste0(code_dir, "OverallXAlgNM.R") )

# Singular_ErrorX <- function(tuningList, Data, initial_Model_Params, weights){
#   return(tryCatch(EMAlg_X(tuningList, Data, initial_Model_Params, weights), error = function(e) NULL))
# }
# 
# Singular_ErrorX2 <- function(initial_Model_Params, Data, tuningpA1, tuningpB1, tuningpA2, tuningpB2){
#   return(tryCatch(EMAlg_X(initial_Model_Params, Data, tuningpA1, tuningpB1, tuningpA2, tuningpB2), error = function(e) NULL))
# }

NMWrapperX <- function(Xinit, Data, tuningvals){
  return(neldermead(tuningvals, EMAlg_X, lower = rep(0, 4), Data = Data, initial_Model_Params = Xinit))
}

### UNCOMMENT THIS CHUNK TO ESTIMATE X ###


numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "jive", "jive.predict", "EMAlg_X", "neldermead"))

Xest <- OverallXAlg(MyData, m_seq = 3:5, u_seq = list(1:3, 2:4))

stopCluster(cl)


#saveRDS(Xest, paste0("test13", sampsize, "rep", datset, ".rds"))

# ### To get estimates for the Y part ...

source(file = paste0(code_dir, "Y_Inits_Cons2.R") )
source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y_Cons2.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlg_YNM.R") )
source(file = paste0(code_dir, "OverallYAlgNM.R") )

### To fix X part at truth when estimating Y part

### COMMENT OUT THIS CHUNK WHEN ESTIMATING X INSTEAD OF FIXING AT TRUTH ###

# Xest <- list()
# X_Param <- Model_Params[c(1:4, 15:17, 21:22)]
# Xest$X_final_pars <- X_Param
# Xest$m_opt <- 4
# Xest$u1_opt <- 2
# Xest$u2_opt <- 3
# Xest$X_log.Lik <- 2
# Xest$diffList <- NA
# Xest$X_BIC <- 432
# Xest$tuningpA1 <- 0
# Xest$tuningpB1 <- 0
# Xest$tuningpA2 <- 0
# Xest$tuningpB2 <- 0
#
# Xest$X_time_diff <- 0

# Singular_ErrorY2 <- function(initial_Model_Params, Data, Xest, tuningpG1, tuningpD1, tuningpG2, tuningpD2){
#   return(tryCatch(EMAlg_Y(initial_Model_Params, Data, Xest, tuningpG1, tuningpD1, tuningpG2, tuningpD2), error = function(e) NULL))
# }

NMWrapperY <- function(Yinit, Data, Xest, tuningvals){
  return(neldermead(tuningvals, EMAlg_Y, lower = rep(0, 4), Data = Data, Xest = Xest, initial_Model_Params = Yinit))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y", "EMAlg_Y", "neldermead"))

Yest <- OverallYAlg(MyData, Xest, j_seq = 4:6, z_seq = list(3:5, 2:4))

stopCluster(cl)


## end of code
saveRDS(Yest, paste0("Yest2rep", datset, "n", sampsize, ".rds"))