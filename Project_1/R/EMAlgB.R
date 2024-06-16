#' Initial B EM Algorithm estimate
#' 
#' An EM algorithm that estimates an unpenalized initial value for B
#' and initial estimates of B0 and Phi1for a specified number of latent factors
#' 
#' @param X The X dataset used to obtain the initial B, B0, and Phi1
#' @param xk_sig The desired number of latent factors of X, corresponding
#' to the number of W variables and the number of columns of B
#' 
#' @return Initial estimates of B, B0, and Phi1 to be used when calculating
#' the penalized estimate of B using the EM algorithm, as well as the
#' total number of iterations until convergence
#' 
#' @examples
#' Binits <- EMAlgB(data$X, 3)
#' 
#' 
#' 
#' @export
EMAlgB <- function(X, xk_sig){
  
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  Xcent <- scale(X, scale = FALSE)
  B0 <- matrix(0, nrow = q, ncol = 1, byrow = T)
  B0new <-  matrix(colMeans(X), nrow = q, ncol = 1, byrow = T)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Bnew <- as.matrix((svd(Xcent)$v%*%diag(svd(Xcent)$d)/sqrt(n))[,1:xk_sig])
  Phi1 <- matrix(0, nrow = q, ncol = q)
  Phi1new <- diag(q)
  
  B31 <- Bnew%*%(eigen(1/(ncol(X))*crossprod(Bnew))$vectors)
  B51 <- B31%*%qr.Q(qr(t(B31[(1:xk_sig),(1:xk_sig)])))
  if(xk_sig == 1){
    Bnew <- B51 %*% solve(B51[(1:xk_sig),(1:xk_sig)])
    Bnew[1,1] <- 1
  }else{
    Bnew <- B51 %*% solve(diag(diag(B51[(1:xk_sig),(1:xk_sig)])))
    Bnew[(1:xk_sig),(1:(xk_sig))][upper.tri(Bnew[(1:xk_sig),(1:(xk_sig))], diag = F)] <- 0
    diag(Bnew[(1:xk_sig),(1:(xk_sig))]) <- 1
  }
  
  while(((norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001)) & (niter < 5000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    invmodv <-Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
    
    EW <- t(B)%*%invmodv%*%(t(X)-matrix(B0, nrow = ncol(X), ncol = n, byrow = F))
    B0new <- (1/n)*as.matrix(colSums(X - t(B%*%EW)))
    Bnew <- (crossprod(X,t(EW)) - matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(EW))%*%solve(n*condvar + tcrossprod(EW))
    Phi1new <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0new, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(Bnew) + matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0new, nrow = q, ncol = n, byrow = F)) + 
                                 2*matrix(B0new, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(Bnew) + n*Bnew%*%condvar%*%t(Bnew) + Bnew%*%EW%*%t(EW)%*%t(Bnew)))
    niter <- niter + 1
  }
  B0 <- B0new
  Phi1 <- Phi1new
  Bem <- Bnew
  B3 <- Bem%*%(eigen(1/(ncol(X))*crossprod(Bem))$vectors)
  B5 <- B3%*%qr.Q(qr(t(B3[(1:xk_sig),(1:xk_sig)])))
  if(xk_sig==1){
    B <- B5 %*% solve(B5[(1:xk_sig),(1:xk_sig)])
    B[1,1] <- 1
  }
  else{
    B <- B5 %*% solve(diag(diag(B5[(1:xk_sig),(1:xk_sig)])))
    B[(1:xk_sig),(1:(xk_sig))][upper.tri(B[(1:xk_sig),(1:(xk_sig))], diag = F)] <- 0
    diag(B[(1:xk_sig),(1:(xk_sig))]) <- 1
  }
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter)
  return(output)
}
