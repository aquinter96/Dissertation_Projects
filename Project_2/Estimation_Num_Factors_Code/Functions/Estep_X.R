## Estep_X.R

Estep_X <- function(Data, Old_Par){
  
  ests <- list()
  
  A1     <-   Old_Par$A1
  A2     <-   Old_Par$A2
  B1     <-   Old_Par$B1
  B2     <-   Old_Par$B2
  
  PhiR   <-   as.matrix(Old_Par$PhiR)
  PhiS1  <-   as.matrix(Old_Par$PhiS1)
  PhiS2  <-   as.matrix(Old_Par$PhiS2)
  Phi11  <-   Old_Par$Phi11
  Phi12  <-   Old_Par$Phi12
  
  m <- ncol(PhiR)
  u1 <- ncol(PhiS1)
  u2 <- ncol(PhiS2)   
  q1 <- ncol(Phi11)
  q2 <- ncol(Phi12)

  PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
  PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
  PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))
  
  
  ests$sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2)), cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2)), cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2)))
  
  #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
  
  ests$sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2)), cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12))
  
  #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
  # 
  ests$sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2)),
                        cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2)),
                        cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2))
  
  #calculate conditional covariance of latent factors given observed data, i.e. Cov(Ei|Fi)
  #and then extract the conditional covariances of each latent factor (Cov(Ri|Fi), Cov(Si1|Fi), etc)
  
  ests$condvar <- ests$sigma22 - ests$sigma21%*%solve(ests$sigma11)%*%t(ests$sigma21)
  
  # calculate conditional expectations of latent factors given observed data
  # 
  # ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
  # Xpanel = list(t(MyData$X1), t(MyData$X2)), Ypanel = list(t(MyData$Y1), t(MyData$Y2))  
  ests$ERS <- ests$sigma21%*%solve(ests$sigma11)%*%rbind( t(Data$X1), t(Data$X2))

  return(ests)
}

## end of code