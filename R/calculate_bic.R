#' Calculate the Bayesian Information Criterion (BIC) for a given model.
#' 
#' This function computes the BIC for a dynamic model represented by the matrices A and S.
#' 
#' @param Y The observed data matrix.
#' @param A The loading matrix.
#' @param S The latent source matrix.
#' 
#' @return The BIC value for the model.

calculate_bic = function(Y, A, S) {
  # calculate BIC
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  # define function to calculate norm
  norm_vec = function(x) sum(x^2)
  
  # calculate sigma
  sigma = sqrt(1 / (n * p) * sum(apply(Y - A %*% S, MARGIN = 1, norm_vec)))
  
  # calculate log likelihood
  loglike = 0
  for(i in 1:n){
    mean = as.vector(A[i,] %*% S)
    loglike = loglike - 2 * sum(log(dnorm(Y[i,], mean, sigma)))
  }
  
  # calculate log(N) * sum(||S||0)
  L11 = log(n) * sum(abs(S) > 1e-1)
  
  # calculate BIC
  bic = loglike + L11
  
  return(bic)
}

