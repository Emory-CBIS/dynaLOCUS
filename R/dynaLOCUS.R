#' A Novel Decomposition Method for Brain Network Dynamic Connectivity Matrices
#' Using Low-Rank Structure with Uniform Sparsity and Temporal Smoothness Regularization
#'
#' @param Y Group-level dynamic connectivity data, represented as a matrix of dimension \code{NT} \times \code{p}, 
#' where \code{NT} denotes the number of subjects multiplied by the number of time points, 
#' and \code{p} represents the number of edges in the connectivity network.
#' @param q The number of connectivity traits to extract.
#' @param V The number of nodes. Note that \code{p} (the columns of \code{Y}) needs to be equal to \code{(V-1)V/2}.
#' @param n_subject The total number of subjects in the study.
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data will be preprocessed. Defaults to \code{TRUE}.
#' @param penalty_S The option for the penalization function for the sparsity regularization for the connectivity traits. 
#' Users can choose "NULL", "L1", or "SCAD". Defaults to "L1". 
#' @param phi A tuning parameter for the element-wise penalty on \code{S}. Defaults to 2.
#' @param rho A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks 
#' for modeling each connectivity trait. The value of \code{rho} represents the closeness of the connectivity traits estimated with and 
#' without the low-rank structure. A higher value of \code{rho} will lead to a higher rank. Defaults to 0.9.
#' @param penalty_A The option for the temporal smoothness regularization for the trait loadings. 
#' If \code{TRUE}, a temporal smoothness penalty is included in the optimization function to encourage similarity in the trait loadings 
#' in adjacent time windows. Defaults to \code{TRUE}.
#' @param lambda A numeric tuning parameter for the temporal smooth lasso penalty on \code{A}. Defaults to 0.1.
#' @param maxIteration The maximum number of iterations. Defaults to 100.
#' @param speed_up If \code{FALSE}, use the Node-rotation algorithm for learning dyna-LOCUS. If \code{TRUE}, use the Alternative algorithm to 
#' further speed up the computation for learning dyna-LOCUS. Defaults to \code{FALSE}.
#' @param espli1 A number describing the tolerance for change on \code{A}.
#' @param espli2 A number describing the tolerance for change on \code{S}.
#' @param demean If \code{TRUE}, demean each column of the group-level connectivity matrix \code{Y}. Defaults to \code{TRUE}.
#' @param silent If \code{TRUE}, suppress the output. Defaults to \code{TRUE}.
#'
#' @return The resulting loading matrix \code{A} and connectivity traits \code{S} based on the decomposition method.
#'
#' @export
#' 
#' @import ica
#' @import MASS

dynaLOCUS = function(Y, q, V, n_subject, 
                     preprocess = TRUE, penalty_S = "L1", phi = 2, rho = 0.9, penalty_A = TRUE, lambda = 0.1,
                     maxIteration = 100, speed_up = FALSE, 
                     espli1 = 0.01, espli2 = 0.05, demean = TRUE, silent = TRUE){
  # demean Y 
  if(demean){
    Y = sweep(Y, 2, apply(Y, 2, mean), "-") 
  }
  
  # total number of observations
  N = dim(Y)[1]
  # number of edges
  K = dim(Y)[2]
  # time length
  t_length = N/n_subject
  # verify number of nodes
  if(V != (sqrt(1+8*K)+1)/2){
    print("V is not correctly specified! Please double check the dimension of your input data.")
    stop()
  }
  
  # preprocess Y based on choice
  if(preprocess){ 
    Yraw = Y
    prep_rslt = dyna_preprocess(Y,q) 
    Y = prep_rslt$Ynew
    H_star = prep_rslt$H_star
  }else{
    H_star = NULL
  }
  
  # initial estimation based on ica
  theta_ini = dyna_initial(Y, q, V, rho = rho)
  A = theta_ini$A
  S = theta_ini$S
  theta = theta_ini$theta
  
  # whether add penalty on S (latent connectivity traitss)
  if(is.null(penalty_S)){
    if(!silent) cat(paste("dyna-LOCUS without penalty on latent connectivity traits \n"))
  }else{
    if(!silent) cat(paste("dyna-LOCUS with", penalty_S, "penalty on latent connectivity traits \n"))
  }
  
  # whether add smooth penalty on A (loading matrix)
  if(!penalty_A){
    if(!silent) cat(paste("dyna-LOCUS without penalty on loading matrix \n"))
  }else{
    if(t_length == 1){
      if(!silent) cat(paste("dyna-LOCUS without penalty on loading matrix \n"))
    }else{
      if(!silent) cat(paste("dyna-LOCUS with smooth penalty on loading matrix \n"))
    }
  }
  
  # initialize list to record sparseness 
  list_sparse = c()
  
  # update parameters
  Iter = 1
  while(Iter <= maxIteration){
    # update 
    theta_new = dyna_update(Y = Y, A = A, S = S, theta = theta, 
                            q = q, K = K, V = V, n_subject = n_subject, t_length = t_length,
                            penalty_S = penalty_S, phi = phi, 
                            penalty_A = penalty_A, H_star = H_star, lambda = lambda, list_sparse = list_sparse,
                            preprocess = preprocess, speed_up = speed_up, silent = silent)
    
    # estimation results from one iteration
    A_new = theta_new$A
    S_new = theta_new$S
    theta_new  = theta_new$theta
    list_sparse = theta_new$list_sparse
    
    # calculate error from ica based S
    errS = norm(as.matrix(S_new-S))/norm(as.matrix(S))
    errA = norm(as.matrix(A_new-A))/norm(as.matrix(A))
    
    # if any NA values generated, break the iteration
    if(sum(is.na(c(errS, errA))) > 0){
      print("Failed to finish!")
      return(list(A = A, S = S, theta = theta))
    }
    
    if(!silent){
      print(paste("Iter ", Iter, "; Percentage change on S: ", round(errS, 3),
                  "; Percentage change on A: ", round(errA, 3), "." , sep = ""))
    }
    
    A = A_new
    S = S_new
    theta = theta_new
    
    if((errA < espli1 & errS < espli2)){
      print("Finished!")
      if(preprocess){
        # transform back to the original scale
        A = H_star %*% A
      }
      
      return(list(A = A, S = S, theta = theta))
    }
    
    Iter = Iter + 1
  }

  print("Finished!")
  if(preprocess){
    # transform back to the original scale
    A = H_star %*% A
  }
  
  return(list(A = A, S = S, theta = theta))
}
