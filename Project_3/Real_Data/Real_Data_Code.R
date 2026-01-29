## Load required packages
library(parallel)

### Random seed 
set.seed(2435)

setwd("~/Paper3_Real")

load(file = "~/Paper3_Real/phenotypic_map/TCGA/groupsBuenoTCGAred_LM.RData")
X <- groupsBuenoTCGAred[!is.na(groupsBuenoTCGAred$TCGA_iCluster),]
X <- X[,c(1, 5, 11:12, 20)]
rownames(X) <- X$Sample
X <- X[,-1]
X <- X[order(row.names(X)),]
X <- as.matrix(X)

load("~/Paper3_Real/phenotypic_map/TCGA/D_met.bodB_MOFA.RData")
load("~/Paper3_Real/phenotypic_map/TCGA/D_met.enhB_MOFA.RData")
load("~/Paper3_Real/phenotypic_map/TCGA/D_met.proB_MOFA.RData")

D_met.bodB_MOFA <- t(D_met.bodB_MOFA)
D_met.enhB_MOFA <- t(D_met.enhB_MOFA)
D_met.proB_MOFA <- t(D_met.proB_MOFA)

D_met.bodB_MOFA <- D_met.bodB_MOFA[,1:500]
D_met.enhB_MOFA <- D_met.enhB_MOFA[,1:500]
D_met.proB_MOFA <- D_met.proB_MOFA[,1:500]

D_met.be <- merge(D_met.bodB_MOFA, D_met.enhB_MOFA, by = 0, all = T)
rownames(D_met.be) <- D_met.be$Row.names
D_met.be <- D_met.be[,-1]
M <- merge(D_met.be, D_met.proB_MOFA, by = 0, all = T)
rownames(M) <- M$Row.names
M <- M[,-1]
M <- M[order(row.names(M)),]
M <- as.matrix(M)

load("~/Paper3_Real/phenotypic_map/TCGA/D_expr_MOFA.RData")
Y <- t(D_expr_MOFA)
Y <- Y[,1:500]
Y <- Y[order(row.names(Y)),]

rm(D_expr_MOFA)
rm(D_met.bodB_MOFA)
rm(D_met.enhB_MOFA)
rm(D_met.proB_MOFA)
rm(D_met.be)
rm(groupsBuenoTCGAred)

MyData <- list()
MyData$X <- X
MyData$M <- M
MyData$Y <- Y

### Load "auxillary" functions

## Functions defined within other functions are not automatically exported to nodes
## during parallelization and instead must be defined separately and exported

source(file = "~/Paper3_Code_AltY/Lasso.R")

### Intialize cluster with 4 nodes and export auxillary functions to nodes
# tottime <- system.time({

### To get estimates for the X part ...

source(file = "~/Paper3_Code_AltY/B_inits.R")
source(file = "~/Paper3_Code_AltY/Estep_X.R")
source(file = "~/Paper3_Code_AltY/Mstep_X.R")
source(file = "~/Paper3_Code_AltY/logLik_X.R")
source(file = "~/Paper3_Code_AltY/Convergence_check.R")
source(file = "~/Paper3_Code_AltY/BIC_X.R")
source(file = "~/Paper3_Code_AltY/EMAlgX.R")
source(file = "~/Paper3_Code_AltY/OverallXAlg.R")

Singular_ErrorX <- function(tuningpB, Data, initial_Model_Params, weights){
  return(tryCatch(EMAlgX(tuningpB, Data, initial_Model_Params, weights), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "EMAlgX"))

Xest <- OverallXAlg(MyData$X, tuningpB = seq(0, 15, 0.1), s_seq = 1:4)

stopCluster(cl)

source(file = "~/Paper3_Code_AltY/A_inits.R")
source(file = "~/Paper3_Code_AltY/Estep_M.R")
source(file = "~/Paper3_Code_AltY/Mstep_M.R")
source(file = "~/Paper3_Code_AltY/logLik_M.R")
source(file = "~/Paper3_Code_AltY/BIC_M.R")
source(file = "~/Paper3_Code_AltY/EMAlgM.R")
source(file = "~/Paper3_Code_AltY/OverallMAlg.R")

Singular_ErrorM <- function(tuningpA, Data, Xest, initial_Model_Params, weights){
  return(tryCatch(EMAlgM(tuningpA, Data, Xest, initial_Model_Params, weights), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_M", "Mstep_M", "logLik_M", "Convergence_check", "BIC_M", "EMAlgM"))

Mest <- OverallMAlg(MyData, Xest, tuningpA = seq(0, 15, 0.1), m_seq = 1:8)

stopCluster(cl)

source(file = "~/Paper3_Code_AltY/C_inits.R")
source(file = "~/Paper3_Code_AltY/Estep_Y.R")
source(file = "~/Paper3_Code_AltY/Mstep_Y.R")
source(file = "~/Paper3_Code_AltY/logLik_Y.R")
source(file = "~/Paper3_Code_AltY/BIC_Y.R")
source(file = "~/Paper3_Code_AltY/EMAlgY.R")
source(file = "~/Paper3_Code_AltY/OverallYAlg.R")

Singular_ErrorY <- function(tuningList, Data, Mest, initial_Model_Params){
  return(tryCatch(EMAlgY(tuningList, Data, Mest, initial_Model_Params), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y", "EMAlgY"))

Yest <- OverallYAlg(MyData, Mest, tuningpC = seq(0, 15, 0.1), tuningpD = seq(0, 15, 0.1), t_seq = 1:8)

stopCluster(cl)

saveRDS(Yest, file=paste0("Paper_Real_Results_500.rds"))