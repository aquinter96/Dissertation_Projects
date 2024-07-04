#Comparison_Code.R

setwd("~/Reformatted_Code")

# Initialize the vectors that, for each dataset replicate, store the number of elements that are equivalent between the parallel and non-parallel results
# for the B, A, and Gamma matrices respectively.

Bmats <- rep(0, 100)
Amats <- rep(0, 100)
Gmats <- rep(0, 100)

# As the results for all tuning parameters tested in the grid search were saved, initialize vectors that will store the BICs for each tuning parameter
# and subsequently used to extract the model with the lowest BIC

pBBICs <- rep(0, 11)
pABICs <- rep(0, 11)
nBBICs <- rep(0, 11)
nABICs <- rep(0, 11)

# For each replicate (once the optimal model is chosen), calculate the MSEs of B, A, and Gamma and store the results in these vectors. Note that whether
# we use the parallel or non-parallel results does not matter, since they are identical

BMSEs <- rep(0, 100)
AMSEs <- rep(0, 100)
GMSEs <- rep(0, 100)

# These vectors will store the runtimes for each replicate for the parallel and non-parallel versions, respectively.

paratime <- rep(0, 100)
normtime <- rep(0, 100)

for(i in 1:100){
  
  paraest <- readRDS(paste("pars100rep", i, "n50", ".rds", sep=""))
  normest <- readRDS(paste("norm100rep", i, "n50", ".rds", sep=""))
  
  for(j in 1:11){
    pBBICs[j] <- paraest[[1]][,j]$B_estimates$B_BIC
    pABICs[j] <- paraest[[1]][,j]$A_estimates$A_BIC
    nBBICs[j] <- paraest[[1]][,j]$B_estimates$B_BIC
    nABICs[j] <- paraest[[1]][,j]$A_estimates$A_BIC
    
  }
  
  paraB <- paraest[[1]][,which(pBBICs == min(pBBICs))][[1]]$B_final_pars$B
  paraA <- paraest[[1]][,which(pABICs == min(pABICs))][[2]]$A_final_pars$A
  paraG <- paraest[[1]][,which(pABICs == min(pABICs))][[2]]$A_final_pars$Gamma
  normB <- paraest[[1]][,which(nBBICs == min(nBBICs))][[1]]$B_final_pars$B
  normA <- paraest[[1]][,which(nABICs == min(nABICs))][[2]]$A_final_pars$A
  normG <- paraest[[1]][,which(nABICs == min(nABICs))][[2]]$A_final_pars$Gamma

 Bmats[i] <- sum(paraB == normB)
 Amats[i] <- sum(paraA == normA)
 Gmats[i] <- sum(paraG == normG)
 
 BMSEs[i] <- 1/294*sum((Model_Params$B - paraB)^2)
 AMSEs[i] <- 1/197*sum((Model_Params$A - paraA)^2)
 GMSEs[i] <- 1/6*sum((Model_Params$Gamma - paraG)^2)
 
 paratime[i] <- paraest[[2]][3]/60
 normtime[i] <- normest[[2]][3]/60
 
}

# end of code