## Model1_Simulation_Code.R
args <- commandArgs(TRUE)
datset <- as.numeric(args[1])
sampsize <- as.numeric(args[2])

## Load required packages
library(parallel)

### Random seed 
set.seed(2435)

### To initialize parameters

source(file = "Model_Meta_Param.R")
source(file = "Model_Matrix_Param.R")

### To generate simulated data 

source(file = "Data_Generation.R")

set.seed(datset)
MyData <- Data_Generation(n = sampsize, Model_Params=Model_Params)


### Load "auxillary" functions

## Functions defined within other functions are not automatically exported to nodes
## during parallelization and instead must be defined separately and exported

source(file = "Lasso.R")

MyData$X <- MyData$M
MyData$M <- MyData$Y
### Intialize cluster with 4 nodes and export auxillary functions to nodes
# tottime <- system.time({
  
  ### To get estimates for the X part ...
  
  source(file = "B_inits.R")
  source(file = "Estep_X.R")
  source(file = "Mstep_X.R")
  source(file = "logLik_X.R")
  source(file = "Convergence_check.R")
  source(file = "BIC_X.R")
  source(file = "EMAlgX.R")
  source(file = "OverallXAlg.R")

  Singular_ErrorX <- function(tuningpB, Data, initial_Model_Params, weights){
    return(tryCatch(EMAlgX(tuningpB, Data, initial_Model_Params, weights), error = function(e) NULL))
  }

  numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
  cl <- makeCluster(as.numeric(numCores))
  clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "EMAlgX"))

  Xest <- OverallXAlg(MyData$X)

  stopCluster(cl)

  source(file = "A_inits.R")
  source(file = "Estep_M.R")
  source(file = "Mstep_M.R")
  source(file = "logLik_M.R")
  source(file = "BIC_M.R")
  source(file = "EMAlgM.R")
  source(file = "OverallMAlg.R")

  Singular_ErrorM <- function(tuningpA, Data, Xest, initial_Model_Params, weights){
    return(tryCatch(EMAlgM(tuningpA, Data, Xest, initial_Model_Params, weights), error = function(e) NULL))
  }

  numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
  cl <- makeCluster(as.numeric(numCores))
  clusterExport(cl, c("lasso", "Estep_M", "Mstep_M", "logLik_M", "Convergence_check", "BIC_M", "EMAlgM"))

  Mest <- OverallMAlg(MyData, Xest)

  stopCluster(cl)

  # source(file = "C_inits.R")
  # source(file = "Estep_Y.R")
  # source(file = "Mstep_Y.R")
  # source(file = "logLik_Y.R")
  # source(file = "BIC_Y.R")
  # source(file = "EMAlgY.R")
  # source(file = "OverallYAlg.R")
  # 
  # Singular_ErrorY <- function(tuningList, Data, Mest, initial_Model_Params){
  #   return(tryCatch(EMAlgY(tuningList, Data, Mest, initial_Model_Params), error = function(e) NULL))
  # }
  # 
  # numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
  # cl <- makeCluster(as.numeric(numCores))
  # clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y", "EMAlgY"))
  # 
  # #Yest <- OverallYAlg(MyData, Mest)
  # Yest <- OverallYAlg(MyData, Model_Params, tuningpC = 0, tuningpD = 0)
  
# })

#Yest$tot.time <- tottime

saveRDS(Mest, file=paste0("Model1MYrep", datset, "n", sampsize, ".rds"))