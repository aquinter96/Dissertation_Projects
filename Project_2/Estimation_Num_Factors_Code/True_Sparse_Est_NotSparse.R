## True_Sparse_Est_NotSparse.R

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
source(file = paste0(code_dir, "Model_Matrix_Param_Sparse.R") )


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

Singular_ErrorX2 <- function(initial_Model_Params, Data, tuningpA1, tuningpB1, tuningpA2, tuningpB2){
  return(tryCatch(EMAlg_X(initial_Model_Params, Data, tuningpA1, tuningpB1, tuningpA2, tuningpB2), error = function(e) NULL))
}

### UNCOMMENT THIS CHUNK TO ESTIMATE X ###

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "jive", "jive.predict", "EMAlg_X"))

Xest <- OverallXAlg(MyData, 0, 0, 0, 0, 1:5, list(1:5, 1:5))

stopCluster(cl)

saveRDS(Xest, paste0("Xhs_lat_estn", sampsize, "rep", datset, ".rds"))

