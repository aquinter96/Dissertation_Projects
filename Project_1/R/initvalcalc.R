#' Initial A EM algorithm estimates
#' 
#' A function that obtains initial estimates of A0, A, Gamma
#' Phi2, and Phi3 using a combination of the nlminb() function
#' and method of moments estimates to extract them from the
#' observed covariance matrix cov(X,Y) and the estimate of B
#' obtained from EMAlgBAdLassoCV
#' 
#' @param X The observed X dataset
#' @param Y The observed Y dataset
#' @param Binit The estimate of B obtained from EMAlgAdLassoCV
#' @param s The number of latent factors for X
#' @param m The number of latent factors for Y
#' @param n The sample size of X and Y
#' 
#' @return Initial estimates of A, Gamma, Phi2, and Phi3
#' 
#' @examples
#' initvalcalc(data$X, data$Y, Best, 3, 2, 400)
#' 
#' 
#' @export
initvalcalc <- function(X, Y, Binit, s, m, n){
  
  AGamma <- t(solve(t(Binit)%*%Binit)%*%t(Binit)%*%cov(X, Y))
  
  p <- ncol(Y)
  
  Aobj <-function(a){
    amat <- matrix(0, nrow = p, ncol = m)
    amat[(1:m),(1:m)][lower.tri(amat[(1:m),(1:m)], diag = T)] <- a[1:(m*(m+1)/2)]
    amat[(m+1):p,] <- a[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
    phi2mat <- diag(a[(p*m-m*(m-1)/2)+(1:p)])
    norm((var(Y) - AGamma%*%t(AGamma) - amat%*%t(amat) - phi2mat), type = "F")
  }
  
  Ainit <- matrix(0, nrow = p, ncol = m)
  if(m==1){
    Ainit[1,1] <- 1
  }
  else{
    diag(Ainit[(1:m),(1:m)]) <- 1
  }
  Phi2init <- diag(p)
  Ainitvec <- c(Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)], Ainit[(m+1):p,])
  Phi2initvec <- diag(Phi2init)
  Anewvec <- nlminb(c(Ainitvec, Phi2initvec), Aobj)
  Ainit <- matrix(0, nrow = p, ncol = m)
  Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)] <- Anewvec$par[1:(m*(m+1)/2)]
  Ainit[(m+1):p,] <- Anewvec$par[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
  
  Phi2init <- diag(Anewvec$par[(p*m-m*(m-1)/2)+(1:p)])
  A3 <- Ainit%*%(eigen(1/n*t(Ainit)%*%Ainit)$vectors)
  A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))
  
  if(m==1){
    Phi3init <- Ainit[1,1]^2
    Ainit <- A5 %*% solve(A5[(1:m),(1:m)])
    Ainit[1,1] <- 1
  }
  else{
    Phi3init <- diag(diag(A5[(1:m),(1:m)]))%*%diag(diag(A5[(1:m),(1:m)]))
    Ainit <- A5 %*% solve(diag(diag(A5[(1:m),(1:m)])))
    Ainit[(1:m),(1:(m))][upper.tri(Ainit[(1:m),(1:(m))], diag = F)] <- 0
    diag(Ainit[(1:m),(1:(m))]) <- 1
  }
  AAinv <- Matrix::solve(t(Ainit)%*%Ainit)
  Gammainit <- AAinv%*%t(Ainit)%*%AGamma
  output <- list("Ainit" = Ainit, "Ginit" = Gammainit, "Phi2init" = Phi2init, "Phi3init" = Phi3init)
  return(output)
}
