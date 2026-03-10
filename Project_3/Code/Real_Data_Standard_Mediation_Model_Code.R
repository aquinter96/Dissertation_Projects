## Load required packages
library(readxl)
library(glmnet)

setwd("~/Paper3_Real")

datfiles <- list.files(pattern=paste0("Paper3_Real_s.*.rds"))

BICmat <- matrix(0, nrow = length(datfiles), ncol = 4)

for(i in 1:length(datfiles)){
  Results <- readRDS(datfiles[i])
  BICmat[i,] <- c(Results$X_estimates$s_opt, Results$M_estimates$m_opt, Results$Y_estimates$t_opt, Results$Y_estimates$Y_BIC)
}

# model with lowest BIC
optBIC <- which(BICmat[,4] == min(BICmat[,4]))
optmod <- readRDS(datfiles[optBIC])

# model with third lowest BIC but better # latent factor estimates, what we
# discussed in our meeting

BICmat <- BICmat[order(BICmat[,4]),]
optlats <- BICmat[3,1:3]
optmod3 <- readRDS(paste0("Paper3_Real_s", optlats[1], "m", optlats[2], "t", optlats[3], ".rds"))

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
varnames <- colnames(XMY)
XMY <- matrix(as.numeric(unlist(XMY)), nrow = 73)
XMY <- XMY[complete.cases(XMY),]
colnames(XMY) <- varnames

X <- XMY[,1:5]
M <- XMY[,6:305]
Y <- XMY[,306:605]

# regress M on X
ind_mod1 <- glmnet(X, M, family = "mgaussian", alpha = 1, lambda = seq(15, 0, -0.1))
# extract coefficient matrix for model with lowest BIC
ind_mod1_BICs <- deviance(ind_mod1) + ind_mod1$df*log(nrow(X))
opt_ind_mod1 <- lapply(ind_mod1$beta, function(x) x[,which(ind_mod1_BICs == min(ind_mod1_BICs))])
opt_ind_mod1_betas <- matrix(unlist(opt_ind_mod1), nrow = 5, ncol = 300)

# regress Y on M
ind_mod2 <- glmnet(M, Y, family = "mgaussian", alpha = 1, lambda = seq(15, 0, -0.1))
# extract coefficient matrix for model with lowest BIC
ind_mod2_BICs <- deviance(ind_mod2) + ind_mod2$df*log(nrow(X))
opt_ind_mod2 <- lapply(ind_mod2$beta, function(x) x[,which(ind_mod2_BICs == min(ind_mod2_BICs))])
opt_ind_mod2_betas <- matrix(unlist(opt_ind_mod2), nrow = 300, ncol = 300)

# regress Y on X and M
dir_mod <- glmnet(cbind(X, M), Y, family = "mgaussian", alpha = 1, lambda = seq(15, 0, -0.1))
# extract coefficient matrices for model with lowest BIC
dir_mod_BICs <- deviance(dir_mod) + dir_mod$df*log(nrow(X))
opt_dir_mod <- lapply(dir_mod$beta, function(x) x[,which(dir_mod_BICs == min(dir_mod_BICs))])
opt_dir_mod_betas <- matrix(unlist(opt_dir_mod), nrow = 305, ncol = 300)
# X coefficient matrix
opt_dir_mod_Xbetas <- opt_dir_mod_betas[1:5,]
# M coefficient matrix
opt_dir_mod_Mbetas <- opt_dir_mod_betas[6:305,]
