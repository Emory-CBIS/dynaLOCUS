#' Preprocess the connectivity data through whitening
#' 
#' @param Ynew A demeaned concatenated group-level connectivity data before whitening
#' @param q The number of latent sources one would like to extract
#' 
#' @return whitened Y (Y tilde) and the dewhitening matrix H_star

dyna_preprocess = function(Ynew, q){
  # get the number of subjects
  N = dim (Ynew)[1]
  # calculate eigen values and vectors of YY'
  eigenA = eigen(Ynew %*% t(Ynew), T)
  # calculate H (whiten mat)
  H = diag((eigenA$values[1:q]-mean(eigenA$values[(q+1):N]))^(-0.5)) %*% t(eigenA$vectors[,1:q])
  Ynew = H %*% Ynew 
  # scale up
  multiplier = 5/sd(Ynew)
  Ynew = Ynew * multiplier
  
  # calculate dewhitening matrix
  H_star = eigenA$vectors[,1:q] %*% diag((eigenA$values[1:q]-mean(eigenA$values[(q+1):N]))^(0.5))/multiplier
  
  # return preprocessed Ynew and dewhitening matrix H
  return(list(Ynew = Ynew, H_star = H_star))
}
