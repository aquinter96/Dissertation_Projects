## AlexProject1_CodeTesting_v01
args <- commandArgs(TRUE)
datset <- as.numeric(args[1])

### Random seed 
library(parallel)

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

MyFB <- function(tuningpB, MyData){

  MyRes <- OverallBAlg(MyData$X, tuningpB, s_seq = 3)
  
  return(MyRes)
  
}

MyFA <- function(tuningpA, MyData, Best){
  
  MyRes <- OverallAGAlg(MyData, Best, tuningpA, m_seq = 2)
  
  return(MyRes)
  
}
### To get estimates for the Y part ...

runtime <- system.time({
  cl <- makeCluster(4) ## keep one core for other things. 
  clusterExport(cl, varlist = c("OverallBAlg", "EMAlgBAdLassoCV", "B_inits", "MyFB", "Estep_X", "Mstep_X", "logLik_X", "BIC_X",
                                "OverallAGAlg", "EMAlgAGammaAdLassoCV", "A_inits", "MyFA", "Estep_Y", "Mstep_Y", "logLik_Y",
                                "BIC_Y", "Convergence_check", "lasso"), envir = environment())
  clusterEvalQ(cl, library(Matrix))
  saveB <- parSapply(cl, seq(0, 10, 1), MyFB, MyData)
  BICs <- rep(0, length(seq(0, 10, 1)))
  for(i in 1:length(seq(0, 10, 1))){
    BICs[i] <- saveB[,i]$B_BIC
  }
  saveA <- parSapply(cl, seq(0, 10, 1), MyFA, MyData, saveB[,which(BICs == min(BICs))])
  stopCluster(cl)
})

saveRDS(list(saveA, runtime), file=paste("pars100rep", datset, "n50", ".rds", sep=""))

## end of code