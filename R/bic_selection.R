#' BIC method for selection of hyper-parameters of phi and rho
#'
#' @param Y Group-level dynamic connectivity data is represented as a matrix of dimension \code{NT}*\code{p}, 
#' where \code{NT} denotes the number of subjects multiplied by the number of time points, 
#' and \code{p} represents the number of edges in the connectivity network.
#' @param q The number of connectivity traits to extract.
#' @param V The number of nodes. Note that \code{p} (the columns of \code{Y}) needs to be equal to \code{(V-1)V/2}.
#' @param n_subject The total number of subjects in the study.
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data will be preprocessed. Defaults to \code{TRUE}.
#' @param penalty The option for the penalization function for the sparsity regularization for the connectivity traits. 
#' The users can choose "NULL"","L1", or "SCAD". Defaults to "SCAD". 
#' @param phi_grid_search Grid search candidates of tuning parameter \code{phi} for penalty on connectivity traits. Defaults to \code{seq(1, 2, 0.2)}. 
#' @param rho_grid_search Grid search candidates of tuning parameter for the adaptive selection method to determine the number of ranks 
#' for modeling each connectivity trait. Defaults to \code{c(0.9)}. 
#' @param maxIteration The maximum number of iterations. Defaults to 100.
#' @param espli1 A number describing the tolerance for change on \code{A}. 
#' @param espli2 A number describing the tolerance for change on \code{S}. 
#' @param demean If \code{TRUE}, performs demeaning on the input data. Defaults to \code{TRUE}.
#' @param save_output If \code{TRUE}, saves the output. Defaults to \code{FALSE}.
#' 
#' @return A list containing the following elements:
#' \item{bic_tab}{A dataframe containing BIC values for each combination of \code{phi} and \code{rho}.}
#' \item{results}{A list of dyna-LOCUS outputs, if \code{save_output} is set to \code{TRUE}.}
#' 
#' @export

bic_selection = function(Y, q, V, n_subject,
                         preprocess = TRUE, penalty = "SCAD", phi_grid_search = seq(1, 2, 0.2), rho_grid_search = c(0.9), 
                         maxIteration = 100, espli1 = 0.01, espli2 = 0.05, demean = TRUE, save_output = FALSE){
  
  # demean the data if specified
  if(demean){
    Y = sweep(Y, 2, apply(Y, 2, mean), "-") 
  }
  
  # preprocess the data if specified
  if(preprocess){
    prep_result = dyna_preprocess(Y, q) 
    Y = prep_result$Ynew
  }
  
  # initialize results and BIC table
  results = list()
  bic_tab = matrix(0, ncol = 3, nrow = length(rho_grid_search) * length(phi_grid_search))
  k = 1
  
  # grid search over rho and phi values
  for(rho in rho_grid_search){
    for(phi in phi_grid_search){
      message('Running dyna-LOCUS for rho = ', rho, ', phi = ', phi, '...')
      
      # run dyna-LOCUS algorithm
      output = dynaLOCUS(Y, q, V, n_subject,
                         preprocess = FALSE, penalty_S = penalty, phi = phi, penalty_A = FALSE, 
                         maxIteration = maxIteration, speed_up = TRUE,
                         espli1 = espli1, espli2 = espli2, rho = rho, silent = TRUE, demean = FALSE)
      
      # calculate BIC value
      bic_value = calculate_bic(Y, output$A, output$S)
      bic_tab[k, ] = c(rho, phi, bic_value)
      
      # save output if specified
      if(save_output){
        results[[k]] = list(output = output, rho = rho, phi = phi)
      }
      
      k = k + 1
    }
  }
  
  # set column names for BIC table
  bic_tab = as.data.frame(bic_tab)
  colnames(bic_tab) = c("rho", "phi", "bic_value")
  
  # return results
  if(save_output){
    return(list(bic_tab = bic_tab, results = results))
  } else {
    return(list(bic_tab = bic_tab))
  }
}
