## AlexProject1_CodeTesting_v01
args <- commandArgs(TRUE)
datset <- as.numeric(args[1])

### Random seed 

## args <- commandArgs(TRUE)
# datset <- 1 ## as.numeric(args[1])
set.seed(2435)


### To initialize parameters

code_dir <- "~/Reformatted_Code/"
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(datset)
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
source(file = paste0(code_dir, "A_inits.R") )
source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlgAGammaAdLassoCV.R") )
source(file = paste0(code_dir, "OverallAGAlg.R") )

runtime <- system.time({
  
  saveB <- OverallBAlg(MyData$X, seq(0, 10, 1), s_seq = 3)

  saveA <- OverallAGAlg(MyData, saveB, seq(0, 10, 1), m_seq = 2)

})

saveRDS(list(saveA, runtime), file=paste("norm100rep", datset, "n50", ".rds", sep=""))

## end of code