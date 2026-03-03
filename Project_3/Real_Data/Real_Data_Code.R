## Load required packages
library(parallel)
library(readxl)
library(nloptr)

args <- commandArgs(TRUE)
s <- as.numeric(args[1])
m <- as.numeric(args[2])
t <- as.numeric(args[3])

### Random seed 
set.seed(2435)

setwd("~/Paper3_Real")

clin.supp <- read_excel("~/Paper3_Real/phenotypic_map/205173_2_supp_5073239_pg14sl.xlsx", sheet = "1B_MPM_Master_Patient_Table")
load(file = "~/Paper3_Real/phenotypic_map/TCGA/groupsBuenoTCGAred_LM.RData")
X <- groupsBuenoTCGAred[!is.na(groupsBuenoTCGAred$TCGA_iCluster),]
X <- merge(X, clin.supp, by.x = "Sample", by.y = "TCGA_barcode", all.x = T)
X <- X[,c(1, 12, 41:43, 46)]
rownames(X) <- X$Sample
X <- X[,-1]
X <- X[order(row.names(X)),]

load("~/Paper3_Real/phenotypic_map/TCGA/D_met.bodB_MOFA.RData")
load("~/Paper3_Real/phenotypic_map/TCGA/D_met.enhB_MOFA.RData")
load("~/Paper3_Real/phenotypic_map/TCGA/D_met.proB_MOFA.RData")

D_met.bodB_MOFA <- t(D_met.bodB_MOFA)
D_met.enhB_MOFA <- t(D_met.enhB_MOFA)
D_met.proB_MOFA <- t(D_met.proB_MOFA)

D_met.bodB_MOFA <- D_met.bodB_MOFA[,1:100]
D_met.enhB_MOFA <- D_met.enhB_MOFA[,1:100]
D_met.proB_MOFA <- D_met.proB_MOFA[,1:100]

D_met.be <- merge(D_met.bodB_MOFA, D_met.enhB_MOFA, by = 0, all = T)
rownames(D_met.be) <- D_met.be$Row.names
D_met.be <- D_met.be[,-1]
M <- merge(D_met.be, D_met.proB_MOFA, by = 0, all = T)
rownames(M) <- M$Row.names
M <- M[,-1]
M <- M[order(row.names(M)),]

load("~/Paper3_Real/phenotypic_map/TCGA/D_expr_MOFA.RData")
Y <- t(D_expr_MOFA)
Y <- Y[,1:300]
Y <- Y[order(row.names(Y)),]

XM <- merge(X, M, by = "row.names")
rownames(XM) <- XM$Row.names
XM <- XM[,-1]

XMY <- merge(XM, Y, by = "row.names")
XMY <- XMY[,-1]
XMY <- matrix(as.numeric(unlist(XMY)), nrow = 73)
XMY <- XMY[complete.cases(XMY),]

X <- XMY[,1:5]
M <- XMY[,6:305]
Y <- XMY[,306:605]

rm(D_expr_MOFA)
rm(D_met.bodB_MOFA)
rm(D_met.enhB_MOFA)
rm(D_met.proB_MOFA)
rm(D_met.be)
rm(groupsBuenoTCGAred)
rm(clin.supp)

MyData <- list()
MyData$X <- X
MyData$M <- M
MyData$Y <- Y

rm(X)
rm(M)
rm(Y)
rm(XM)
rm(XMY)

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

NMWrapperXinit <- function(lat_grid, X){
  return(tryCatch(B_inits(lat_grid, X = X), error = function(e) NULL))
}

NMWrapperX <- function(Xinit, X, tuningvals){
  return(tryCatch(sbplx(tuningvals, NMWrapperXEM, lower = 0, X = X, initial_Model_Params = Xinit), error = function(e) NULL))
}

NMWrapperXEM <- function(tuningvals, X, initial_Model_Params){
  return(tryCatch(EMAlgX(tuningvals = tuningvals, Data = X, initial_Model_Params = initial_Model_Params, NM = T), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "EMAlgX", "sbplx", "NMWrapperXEM"))

Xest <- OverallXAlg(MyData$X, s_seq = s)

stopCluster(cl)

source(file = "~/Paper3_Code_AltY/A_inits.R")
source(file = "~/Paper3_Code_AltY/Estep_M.R")
source(file = "~/Paper3_Code_AltY/Mstep_M.R")
source(file = "~/Paper3_Code_AltY/logLik_M.R")
source(file = "~/Paper3_Code_AltY/BIC_M.R")
source(file = "~/Paper3_Code_AltY/EMAlgM.R")
source(file = "~/Paper3_Code_AltY/OverallMAlg.R")

NMWrapperMinit <- function(lat_grid, Data, Xest){
  return(tryCatch(A_inits(lat_grid, Data = Data, Xest = Xest), error = function(e) NULL))
}

NMWrapperM <- function(Minit, Data, Xest, tuningvals){
  return(tryCatch(sbplx(tuningvals, NMWrapperMEM, lower = 0, control = list(maxeval = 10000), Data = Data, Xest = Xest, initial_Model_Params = Minit), error = function(e) NULL))
}

NMWrapperMEM <- function(tuningvals, Data, Xest, initial_Model_Params){
  return(tryCatch(EMAlgM(tuningvals = tuningvals, Data = Data, Xest = Xest, initial_Model_Params = initial_Model_Params, NM = T), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_M", "Mstep_M", "logLik_M", "Convergence_check", "BIC_M", "EMAlgM"))

Mest <- OverallMAlg(MyData, Xest, m_seq = m)

stopCluster(cl)

source(file = "~/Paper3_Code_AltY/C_inits.R")
source(file = "~/Paper3_Code_AltY/Estep_Y.R")
source(file = "~/Paper3_Code_AltY/Mstep_Y.R")
source(file = "~/Paper3_Code_AltY/logLik_Y.R")
source(file = "~/Paper3_Code_AltY/BIC_Y.R")
source(file = "~/Paper3_Code_AltY/EMAlgY.R")
source(file = "~/Paper3_Code_AltY/OverallYAlg.R")

NMWrapperYinit <- function(lat_grid, Data, Mest){
  return(tryCatch(C_inits(lat_grid, Data = Data, Mest = Mest), error = function(e) NULL))
}

NMWrapperY <- function(Yinit, Data, Mest, tuningvals){
  return(tryCatch(sbplx(tuningvals, NMWrapperYEM, lower = rep(0, 2), control = list(maxeval = 10000), Data = Data, Mest = Mest, initial_Model_Params = Yinit), error = function(e) NULL))
}

NMWrapperYEM <- function(tuningvals, Data, Mest, initial_Model_Params){
  return(tryCatch(EMAlgY(tuningvals = tuningvals, Data = Data, Mest = Mest, initial_Model_Params = initial_Model_Params, NM = T), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y", "EMAlgY"))

Yest <- OverallYAlg(MyData, Mest, t_seq = t)

stopCluster(cl)

saveRDS(Yest, file=paste0("Paper3_Real_s", s, "m", m, "t", t, ".rds"))