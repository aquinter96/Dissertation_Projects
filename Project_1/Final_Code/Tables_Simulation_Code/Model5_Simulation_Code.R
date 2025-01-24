## Model5_Simulation_Code.R
args <- commandArgs(TRUE)
datset <- as.numeric(args[1])
sampsize <- as.numeric(args[2])

## Load required packages
library(parallel)

### Random seed 
set.seed(2435)


### To initialize parameters

source(file = "Model5_Meta_Param.R")
source(file = "Model45_Matrix_Param.R")


### To generate simulated data 

source(file = "Data_Generation.R")

set.seed(datset)
MyData <- Data_Generation(n = sampsize, Model_Params=Model_Params)


### Load "auxillary" functions

## Functions defined within other functions are not automatically exported to nodes
## during parallelization and instead must be defined separately and exported

source(file = "Lasso.R")

### Intialize cluster with 4 nodes and export auxillary functions to nodes
tottime <- system.time({
  
  ### To get estimates for the X part ...
  
  source(file = "B_inits.R")
  source(file = "Estep_X.R")
  source(file = "Mstep_X.R")
  source(file = "logLik_X.R")
  source(file = "Convergence_check.R")
  source(file = "BIC_X.R")
  source(file = "EMAlgBAdLassoCV.R")
  source(file = "OverallBAlg.R")
  
  numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
  cl <- makeCluster(as.numeric(numCores))
  clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X"))
  
  Best <- OverallBAlg(MyData$X)
  
  stopCluster(cl)
  
  source(file = "A_inits.R")
  source(file = "Estep_Y.R")
  source(file = "Mstep_Y.R")
  source(file = "logLik_Y.R")
  source(file = "BIC_Y.R")
  source(file = "EMAlgAGammaAdLassoCV.R")
  source(file = "OverallAGAlg.R")
  
  numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
  cl <- makeCluster(as.numeric(numCores))
  clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y"))
  
  Aest <- OverallAGAlg(MyData, Best)
  
  stopCluster(cl)
  
})

saveRDS(Aest, file=paste0("Model5rep", datset, "n", sampsize, ".rds"))