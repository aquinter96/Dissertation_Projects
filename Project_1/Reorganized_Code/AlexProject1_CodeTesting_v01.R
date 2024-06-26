## AlexProject1_CodeTesting_v01

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


### To get estimates for the X part ...

source(file = paste0(code_dir, "Lasso.R") )
source(file = paste0(code_dir, "B_inits.R") )
source(file = paste0(code_dir, "Estep_X.R") )
source(file = paste0(code_dir, "Mstep_X.R") )
source(file = paste0(code_dir, "logLik_X.R") )
source(file = paste0(code_dir, "Convergence_check.R") )
source(file = paste0(code_dir, "BIC_X.R") )
source(file = paste0(code_dir, "EMAlgBAdLassoCV.R") )
source(file = paste0(code_dir, "OverallBAlg.R") )

Best <- OverallBAlg(MyData$X, 1:3, 2:3)


### To get estimates for the Y part ...

source(file = paste0(code_dir, "A_inits.R") )
source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlgAGammaAdLassoCV.R") )
source(file = paste0(code_dir, "OverallAGAlg.R") )

Aest <- OverallAGAlg(MyData, Best, 1, 2)

## end of code