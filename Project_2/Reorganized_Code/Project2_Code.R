## AlexProject1_CodeTesting_v01

### Random seed 

library(pracma)
library(MASS)
library(r.jive)

## args <- commandArgs(TRUE)
# datset <- 1 ## as.numeric(args[1])
set.seed(2435)


### To initialize parameters

code_dir <- "C:/Users/alex/OneDrive/Documents/GitHub/Dissertation_Projects/Project_2/Reorganized_Code/"
source(file = paste0(code_dir, "Model_Meta_Param.R") )
source(file = paste0(code_dir, "Model_Matrix_Param.R") )


### To generate simulated data 

source(file = paste0(code_dir, "Data_Generation.R") )

set.seed(1)
MyData <- Data_Generation(n = 100, Model_Params=Model_Params)


### To get estimates for the X part ...

source(file = paste0(code_dir, "Lasso.R") )
source(file = paste0(code_dir, "X_inits.R") )
source(file = paste0(code_dir, "Estep_X.R") )
source(file = paste0(code_dir, "Mstep_X.R") )
source(file = paste0(code_dir, "logLik_X.R") )
source(file = paste0(code_dir, "Convergence_check.R") )
source(file = paste0(code_dir, "BIC_X.R") )
source(file = paste0(code_dir, "EMAlg_X.R") )
source(file = paste0(code_dir, "OverallXAlg.R") )

Xest <- OverallXAlg(MyData, tuningpA1 = c(0, 1), tuningpB1 = c(0, 1), tuningpA2 = c(0, 1), tuningpB2 = c(0, 1), m_seq = c(2, 3), u_seq = list(2, c(2, 3)))


### To get estimates for the Y part ...

source(file = paste0(code_dir, "Y_inits.R") )
source(file = paste0(code_dir, "Estep_Y.R") )
source(file = paste0(code_dir, "Mstep_Y.R") )
source(file = paste0(code_dir, "logLik_Y.R") )
source(file = paste0(code_dir, "BIC_Y.R") )
source(file = paste0(code_dir, "EMAlg_Y.R") )
source(file = paste0(code_dir, "OverallYAlg.R") )

Yest <- OverallYAlg(MyData, Xest, tuningpG1 = c(0, 1), tuningpD1 = c(0, 1), tuningpG2 = c(0, 1), tuningpD2 = c(0, 1), j_seq = 2, z_seq = list(c(2, 3), 3))

## end of code