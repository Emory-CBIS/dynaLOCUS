#' Goodness-of-Fit Method for Hyper-Parameter Lambda Selection
#'
#' @param Y Group-level dynamic connectivity data represented as a matrix of dimension \code{NT \times p}, 
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
#' @param lambda_grid_search A vector of candidate values for the tuning parameter \code{lambda}. 
#' Defaults to \code{c(NA, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)}, where \code{NA} represents no temporal smoothness penalty.
#' @param maxIteration The maximum number of iterations. Defaults to 100.
#' @param espli1 A number describing the tolerance for change on \code{A}. 
#' @param espli2 A number describing the tolerance for change on \code{S}. 
#' @param demean If \code{TRUE}, demean each column of the group-level connectivity matrix \code{Y}. Defaults to \code{TRUE}.
#' @param save_output If \code{TRUE}, saves the output. Defaults to \code{FALSE}.
#' 
#' @return A list containing the following elements:
#' \item{selected_lambda}{The lambda value that achieves the minimal goodness-of-fit.}
#' \item{goodness_of_fit}{Goodness-of-fit values for each \code{lambda}.}
#' \item{results}{A list of dyna-LOCUS outputs, if \code{save_output} is set to \code{TRUE}.}
#' 
#' @export

lambda_selection = function(Y, q, V, n_subject,
                            preprocess = TRUE, penalty_S = "L1", phi = 2, rho = 0.9,
                            lambda_grid_search = c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1), 
                            maxIteration = 100, espli1 = 0.01, espli2 = 0.05, demean = TRUE, save_output = FALSE){
  
  # demean and preprocess Y for further usage
  if(demean){
    Y = sweep(Y, 2, apply(Y, 2, mean), "-") 
  }

  # run grid search  
  results = list()
  goodness_of_fit = data.frame(lambda = lambda_grid_search, goodness_of_fit = NA)
  k = 1
  for(lambda in lambda_grid_search){
    if(!is.na(lambda)){
      message('Running dyna-locus for lambda = ', lambda, '...')
      result = dynaLOCUS(Y, q, V, n_subject,
                          preprocess = preprocess, penalty_S = penalty_S, phi = phi, penalty_A = TRUE, lambda = lambda,
                          maxIteration = maxIteration, speed_up = TRUE,
                          espli1 = espli1, espli2 = espli2, rho = rho, silent = TRUE, demean = FALSE)
    }else{
      message('Running dyna-locus with no temporal smoothness penalty')
      result = dynaLOCUS(Y, q, V, n_subject,
                          preprocess = preprocess, penalty_S = penalty_S, phi = phi, penalty_A = FALSE, 
                          maxIteration = maxIteration, speed_up = TRUE,
                          espli1 = espli1, espli2 = espli2, rho = rho, silent = TRUE, demean = FALSE)
    }
    # calculate the log reconstruction err
    A = result$A
    S = result$S
    goodness_of_fit[k,] = c(lambda, sum((Y - A %*% S)^2))
    if(save_output){
      results[[k]] = list(result = result, lambda = lambda)
    }
    k = k + 1
  }
  
  # selected lambda
  selected_lambda = goodness_of_fit$lambda[which.min(goodness_of_fit$goodness_of_fit)]
  
  if(save_output){
    return(list(selected_lambda = selected_lambda, goodness_of_fit = goodness_of_fit, results = results))
  }
  
  return(list(selected_lambda = selected_lambda, goodness_of_fit = goodness_of_fit))
}
