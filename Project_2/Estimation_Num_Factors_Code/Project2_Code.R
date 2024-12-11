## AlexProject1_CodeTesting_v01

### Random seed 

library(pracma)
library(MASS)
library(r.jive)
library(parallel)

args <- commandArgs(TRUE)
datset <- as.numeric(args[1])
sampsize <- as.numeric(args[2])
set.seed(2435)


### To initialize parameters

code_dir <- paste0("~/Paper2_Parallel/")
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


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
source(file = paste0(code_dir, "EMAlg_X2.R") )
source(file = paste0(code_dir, "OverallXAlg2.R") )

### UNCOMMENT THIS CHUNK TO ESTIMATE X ###

# numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
# cl <- makeCluster(as.numeric(numCores))
# clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "jive", "jive.predict"))
# 
# Xest <- OverallXAlg(MyData, 0, 0, 0, 0, 1:5, list(1:5, 1:5))
# 
# stopCluster(cl)

#saveRDS(Xest, paste0("X2_lat_estn", sampsize, "rep", datset, ".rds"))

### To get estimates for the Y part ...

source(file = paste0(code_dir, "Y_Inits2.R") )
source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlg_Y2.R") )
source(file = paste0(code_dir, "OverallYAlg2.R") )

### To fix X part at truth when estimating Y part

### COMMENT OUT THIS CHUNK WHEN ESTIMATING X INSTEAD OF FIXING AT TRUTH ###

Xest <- list()
X_Param <- Model_Params[c(1:4, 15:17, 21:22)]
Xest$X_final_pars <- X_Param
Xest$m_opt <- 4
Xest$u1_opt <- 2
Xest$u2_opt <- 3
Xest$X_log.Lik <- 2
Xest$diffList <- NA
Xest$X_BIC <- 432
Xest$tuningpA1 <- 0
Xest$tuningpB1 <- 0
Xest$tuningpA2 <- 0
Xest$tuningpB2 <- 0

Xest$X_time_diff <- 0

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c( "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y"))

Yest <- OverallYAlg(MyData, Xest, 0, 0, 0, 0, j_seq = 1:5, z_seq = list(1:5, 1:5))

stopCluster(cl)
## end of code

saveRDS(Yest, paste0("Y2_lat_estn", sampsize, "rep", datset, ".rds"))