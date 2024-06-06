#' Generate starting values for the iteration
#' 
#' @param Y A demeaned concatenated group-level connectivity data after whitening if preprocess = TRUE
#' @param q The number of latent sources one would like to extract
#' @param V The number of nodes in the study
#' @param rho A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks 
#' for modeling each connectivity trait. The value of rho represents the closeness of the connectivity traits estimated with and 
#' without the low-rank structure. A higher value of rho will lead to a higher rank. Defaults to 0.9.
#' @param R A numeric vector of pre-specifying ranks of latent sources. Defaults to NULL.
#' @param maxIter The maximum iteration number when doing ica
#' 
#' @import ica
#' 
#' @return Initial values of loading matrix \code{A} and latent source \code{S}

dyna_initial = function(Y, q, V, rho = 0.9, R = NULL, maxIter = 100){ 
  # initial decomposition through ica
  ICcorr = ica::icaimax(t(Y), nc = q, center = F, maxit = maxIter)
  # construct initial value of S and theta
  S_ini = matrix(0, ncol = dim(ICcorr$S)[1], nrow=q)
  theta_ini = list()
  # loop over q latent components
  for(i in 1:q){
    # for each component, perform eigen decomposition
    Sl = Ltrinv(ICcorr$S[,i], V, F)
    Sl = Sl + diag(rep(mean(ICcorr$S[,i]), V))
    eigenSl = eigen(Sl)
    # order eigen values from large to small
    orderEigen = order(abs(eigenSl$values),decreasing = T)
    # reconstruct each component Sl using eigen vectors and eigen values
    # stop when similarity to raw Sl larger than rho
    if(is.null(R)){
      Rl = 2
      while(TRUE){
        eigenset = orderEigen[1:Rl]
        imgeRL = eigenSl$vectors[,eigenset] %*% diag(eigenSl$values[eigenset]) %*% 
          t(eigenSl$vectors[,eigenset])
        if(cor(Ltrans(imgeRL, F), ICcorr$S[,i]) > rho) break
        Rl = Rl + 1
      }
    }else{
      Rl = R[i]
      eigenset = orderEigen[1:Rl]
    }
    
    # store information in theta
    theta_ini[[i]] = list()
    # store eigen values in lam_l
    theta_ini[[i]]$lam_l = eigenSl$values[eigenset]
    # change the sign of the largest eigen value to positive
    # sign of other eigen values are changed respectively
    if(theta_ini[[i]]$lam_l[1] < 0){
      theta_ini[[i]]$lam_l = -1*theta_ini[[i]]$lam_l
    }   
    # store corresponding eigen vectors
    theta_ini[[i]]$X_l = matrix(0, ncol = V, nrow = Rl)
    for(j in 1:Rl){
      theta_ini[[i]]$X_l[j,] = eigenSl$vectors[,eigenset[j]]
    }
    # store reconstructed Sl
    S_ini[i,] = Ltrans(t(theta_ini[[i]]$X_l) %*% diag(theta_ini[[i]]$lam_l) %*% 
                         theta_ini[[i]]$X_l, F)
  }
  
  # compute initial value of A
  A_ini = Y%*%t(S_ini)%*%solve(S_ini%*%t(S_ini))

  # scale up
  for(l in 1:q){
    # unit norm each column of A
    # scaleL = sqrt(sum(A_ini[,l]^2))
    scaleL = sd(A_ini[,l])
    A_ini[,l] = A_ini[,l] / scaleL
    S_ini[l,] = S_ini[l,] * scaleL
    # scale X_l correspondingly
    theta_ini[[l]]$X_l = theta_ini[[l]]$X_l * sqrt(scaleL)
  }
  
  # since after preprocessing, A_tilde is orthogonal
  # compute A_tilde transpose/inverse (g-inverse of A-tilde)
  # if X has full column rank, (X'X)^(-1)X' is its g-inverse
  # why use g-inverse here: since A_ini is N*q
  # afterwards, in the update process, we actually can just use A_tilde transpose(after scale), or g-inverse
  M_ini =  solve(t(A_ini)%*%A_ini)%*%t(A_ini) # g-inverse of A
  for(l in 1:q){
    theta_ini[[l]]$M_l = M_ini[l,]
  }
  
  return(list(A = A_ini, theta = theta_ini, S = S_ini))
}

