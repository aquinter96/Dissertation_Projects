## AlexProject1_CodeTesting_v01

### Random seed 

library(pracma)
library(MASS)
library(r.jive)
library(parallel)

args <- commandArgs(TRUE)
datset <- as.numeric(args[1])

set.seed(2435)

### To initialize parameters

code_dir <- paste0("~/Paper2_Parallel/")
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(datset)
MyData <- Data_Generation(n = 1000, Model_Params=Model_Params)


### To get estimates for the X part ...

source(file = paste0(code_dir, "Lasso.R") )
source(file = paste0(code_dir, "Estep_X.R") )
source(file = paste0(code_dir, "Mstep_X.R") )
source(file = paste0(code_dir, "logLik_X.R") )
source(file = paste0(code_dir, "Convergence_check.R") )
source(file = paste0(code_dir, "BIC_X.R") )
source(file = paste0(code_dir, "EMAlg_X.R") )
source(file = paste0(code_dir, "OverallXAlg0.R") )

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X"))

Xest <- OverallXAlg(Model_Params[c(1:4, 15:17, 21:22)], MyData, 0:3, 0:3, 0:3, 0:3)

stopCluster(cl)


### To get estimates for the Y part ...

source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlg_Y.R") )
source(file = paste0(code_dir, "OverallYAlg0.R") )

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y"))

Yest <- OverallYAlg(Model_Params[c(5:11, 13, 12, 14, 23:27)], MyData, Xest, tuningpG1 = 0:3, tuningpD1 = 0:3, tuningpG2 = 0:3, tuningpD2 = 0:3)

stopCluster(cl)
## end of code

saveRDS(Yest, paste0("Sim", datset, ".rds"))