# AlexProject1_CodeTesting_Parallel.R

### Load libraries

library(parallel)


### Random seed 


## args <- commandArgs(TRUE)
# datset <- 1 ## as.numeric(args[1])
set.seed(2435)


### To initialize parameters

code_dir <- "C:/Users/alex/OneDrive/Documents/GitHub/Dissertation_Projects/Project_1/Reorganized_Code/"
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(1)
MyData <- Data_Generation(n = 50, Model_Params=Model_Params)


### Load "auxillary" functions

## Functions defined within other functions are not automatically exported to nodes
## during parallelization and instead must be defined separately and exported

source(file = paste0(code_dir, "Lasso.R") )
source(file = paste0(code_dir, "TuneSearchB.R") )
source(file = paste0(code_dir, "TuneSearchA.R") )


### Intialize cluster with 4 nodes and export auxillary functions to nodes

cl <- makeCluster(4)

clusterExport(cl, c("lasso", "tunesearchB", "tunesearchA"))


### To get estimates for the X part ...

source(file = paste0(code_dir, "B_inits.R") )
source(file = paste0(code_dir, "EMAlgBAdLassoCV.R") )
source(file = paste0(code_dir, "OverallBAlg.R") )

Best <- OverallBAlg(MyData$X)


### To get estimates for the Y part ...

source(file = paste0(code_dir, "initvalcalc.R") )
source(file = paste0(code_dir, "EMAlgAGammaAdLassoCV.R") )
source(file = paste0(code_dir, "OverallAGAlg.R") )

Aest <- OverallAGAlg(MyData$X, MyData$Y, Best)


### Stop cluster once the calculations are finished

stopCluster(cl)

## end of code