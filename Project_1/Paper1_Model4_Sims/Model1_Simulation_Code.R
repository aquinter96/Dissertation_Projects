## Model1_Simulation_Code.R
args <- commandArgs(TRUE)
datset <- as.numeric(args[1])
sampsize <- as.numeric(args[2])

## Load required packages
library(parallel)

### Random seed 
set.seed(2435)


### To initialize parameters

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")


### To generate simulated data 

source(file = "Data_Generation.R")

set.seed(datset)
MyData <- Data_Generation(n = sampsize, Model_Params=Model_Params)


### Load "auxillary" functions

## Functions defined within other functions are not automatically exported to nodes
## during parallelization and instead must be defined separately and exported

source(file = "Lasso.R")
source(file = "TuneSearchB.R")
source(file = "TuneSearchA.R")


### Intialize cluster with 4 nodes and export auxillary functions to nodes

cl <- makeCluster(detectCores())

clusterExport(cl, c("lasso", "tunesearchB", "tunesearchA"))


### To get estimates for the X part ...

source(file = "B_inits.R")
source(file = "Estep_X.R")
source(file = "Mstep_X.R")
source(file = "logLik_X.R")
source(file = "Convergence_check.R")
source(file = "BIC_X.R")
source(file = "EMAlgBAdLassoCV.R")
source(file = "OverallBAlg.R")

Best <- OverallBAlg(MyData$X)


### To get estimates for the Y part ...

source(file = "A_inits.R")
source(file = "Estep_Y.R")
source(file = "Mstep_Y.R")
source(file = "logLik_Y.R")
source(file = "BIC_Y.R")
source(file = "EMAlgAGammaAdLassoCV.R")
source(file = "OverallAGAlg.R")

Aest <- OverallAGAlg(MyData, Best)


### Stop cluster once the calculations are finished

stopCluster(cl)

saveRDS(Aest, file=paste("Model1rep", datset, "n", sampsize, ".rds", sep=""))

## end of code