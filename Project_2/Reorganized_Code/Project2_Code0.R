## AlexProject1_CodeTesting_v01

### Random seed 

library(pracma)
library(MASS)
library(r.jive)
library(parallel)

## args <- commandArgs(TRUE)
# datset <- 1 ## as.numeric(args[1])
set.seed(2435)


### To initialize parameters

code_dir <- paste0("~/Paper2_Parallel/")
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(1)
MyData <- Data_Generation(n = 100, Model_Params=Model_Params)


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

Xest <- OverallXAlg0(Model_Params[c(1:4, 15:17, 21:22)], MyData, 1:4, 1:4, 1:4, 1:4)

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

Yest <- OverallYAlg0(Model_Params[c(5:11, 13, 12, 14, 23:27)], MyData, Xest, tuningpG1 = 1:4, tuningpD1 = 1:4, tuningpG2 = 1:4, tuningpD2 = 1:4)

## end of code

Ytest <- saveRDS(Yest, "XYtestout.rds")