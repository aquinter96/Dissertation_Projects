## Load required packages
library(parallel)
library(readxl)
library(nloptr)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(rtracklayer)

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

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# cpg_coords <- ann450k[colnames(M), c("chr", "pos", "strand")]
# 
# hg19_granges <- GRanges(seqnames = cpg_coords$chr,
#                 ranges = IRanges(cpg_coords$pos, width = 1),
#                 strand = cpg_coords$strand, probeID = rownames(cpg_coords))
# 
# chain <- rtracklayer::import.chain("hg19ToHg38.over.chain")
# hg38_granges <- rtracklayer::liftOver(hg19_granges, chain)

# mapping cpg sites to genes
genes_map <- ann450k[colnames(M), c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
genes_map$UCSC_RefGene_Name[optmod$M_estimates$M_final_pars$A[,1] != 0]

# mapping RNA expression to genes
ensmbl_ids <- gsub("\\..*", "", colnames(Y))
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ensmbl_ids, column = "SYMBOL",
                keytype = "ENSEMBL", multiVals = "list")
gene_symbols[is.na(gene_symbols)] <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ensmbl_ids, column = "SYMBOL",
                                     keytype = "ALIAS", multiVals = "list")[is.na(gene_symbols)]