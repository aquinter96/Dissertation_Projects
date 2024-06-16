#' Generate Simulated Dataset
#' 
#' Generate 2 datasets, an X data and Y dataset, given a
#' prespecified set of parameters
#' 
#' @param ak A px(m+1) matrix where the first column is the true A0 and the remaining m columns correspond to the true A
#' @param bk A qx(s+1) matrix where the first column is the true B0 and the remaining s columns correspond to the true B
#' @param gamma An mxs matrix goving the true Gamma matrix
#' @param n The desired sample size of the X and Y datasets
#' @param varvecB The q length vector giving the diagonal elements of Phi1
#' @param varvecA The p length vector giving the diagonal elements of Phi2
#' @param varvecG The m length vector giving the diagonal elements of Phi3
#' 
#' @return The generated X and Y datasets, as well as the datasets for the latent variables W and Z for debugging purposes
#' 
#' @examples
#' # data = datgen(trueA, trueB, trueGamma, sampsize, diag(truePhi1), diag(truePhi2), diag(truePhi3))
#' 
#' 
#' 
#' 
#' @export
datgen <- function(ak, bk, gamma, n,varvecB,varvecA,varvecG){
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    W <- matrix(rnorm(n, mean = 0, sd = 1), nrow = 1, n) ## generate matrix U N(0,I) m x n
  }
  else{
    W <- matrix(rnorm(n*ncol(gamma), mean = 0, sd = 1), nrow = ncol(gamma), n) ## generate matrix U N(0,I) m x n
  }
  xlin.term  <-  bk %*% rbind(rep(1, ncol(W)),W)
  xlin.term <- t(xlin.term)
  xx <- matrix(0, nrow = n, ncol = nrow(bk))
  for(jj in 1: (ncol(xx))){
    xx[, jj] <- xlin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecB[jj])))
  }
  colnames(xlin.term) <- paste0('X', 1:(ncol(xx)))
  colnames(xx) <- paste0('X', 1:(ncol(xx)))
  X = xx
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    Z <- gamma*W + matrix(rnorm(n, mean = 0, sqrt(varvecG)), nrow = 1) ## generate latent matrix Z
  }
  else{
    Z <- matrix(0, nrow=nrow(gamma), ncol = n)
    for(i in 1:nrow(gamma)){
      Z[i,] <- gamma[i,]%*%W + matrix(rnorm(n, mean = 0, sd = sqrt(varvecG[i])), nrow = 1) ## generate latent matrix Z
    }
  }
  ylin.term  <-  ak %*% rbind(rep(1, ncol(Z)),Z)
  ylin.term <- t(ylin.term)
  yy <- matrix(0, nrow = n, ncol = nrow(ak))
  for(jj in 1: (ncol(yy))){
    yy[, jj] <- ylin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecA[jj])))
  }
  colnames(ylin.term) <- paste0('Y', 1:(ncol(yy)))
  colnames(yy) <- paste0('Y', 1:(ncol(yy)))
  Y = yy
  output <- list("X" = X, "Y" = Y, "W" = W, "Z" = Z)
}
