#' Reliability Index Calculation
#'
#' This function calculates the reliability index and reproducibility (RI) for the dyna-LOCUS method.
#'
#' @param Y Group-level dynamic connectivity data, represented as a matrix of dimension \code{NT} \times \code{p}, 
#' where \code{NT} denotes the number of subjects multiplied by the number of time points, 
#' and \code{p} represents the number of edges in the connectivity network.
#' @param q_choices A vector of integers representing different choices for the number of connectivity traits to extract.
#' @param V The number of nodes. Note that \code{p} (the columns of \code{Y}) needs to be equal to \code{(V-1)V/2}.
#' @param n_subject The total number of subjects in the study.
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data will be preprocessed. Defaults to \code{TRUE}.
#' @param penalty_S The option for the penalization function for the sparsity regularization for the connectivity traits. 
#' Users can choose "NULL", "L1", or "SCAD". Defaults to "L1". 
#' @param phi A tuning parameter for the element-wise penalty on S. Defaults to 2.
#' @param rho A proportional tuning parameter ranging from 0 to 1 for the adaptive selection method to determine the number of ranks 
#' for modeling each connectivity trait. The value of \code{rho} represents the closeness of the connectivity traits estimated with and 
#' without the low-rank structure. A higher value of \code{rho} will lead to a higher rank. Defaults to 0.9.
#' @param penalty_A The option for the temporal smoothness regularization for the trait loadings. 
#' If \code{TRUE}, a temporal smoothness penalty is included in the optimization function to encourage similarity in the trait loadings 
#' in adjacent time windows. Defaults to \code{TRUE}.
#' @param lambda A numeric tuning parameter for the temporal smooth lasso penalty on \code{A}. Defaults to 0.1.
#' @param maxIteration The maximum number of iterations for the dyna-LOCUS algorithm.
#' @param espli1 A number describing the tolerance for change on A. 
#' @param espli2 A number describing the tolerance for change on S. 
#' @param silent Logical, if FALSE, print out the penalty added on A and S.
#' @param demean Logical, if TRUE, demean the data before processing.
#' @param seeds A vector of integers to set seeds for bootstrap sampling, default is 1:50.
#' 
#' @return A list containing:
#' \item{reliability}{A list where each component corresponds to a different choice of q. Each component is of length q, 
#' representing the reliability index for each trait based on that choice of q.}
#' \item{RI}{A list where each component corresponds to a different choice of q. Each component is of length q, 
#' representing the reproducibility (RI) for each trait based on that choice of q.}
#' 
#' @export

reliability_index = function(Y, q_choices, V, n_subject, 
                             preprocess = TRUE, penalty_S = "L1", phi = 2, rho = 0.9, penalty_A = TRUE, lambda = 0.1,
                             maxIteration = 100, espli1 = 0.01, espli2 = 0.05, silent = TRUE, demean = TRUE, 
                             seeds = 1:50){
  
  # total number of observations
  N = dim(Y)[1]
  # time length per subject
  t_length = N / n_subject
  
  # initialize lists for RI and reliability
  RI = list()
  reliability = list()
  
  # Loop over each choice of q
  for (q in q_choices) {
    
    # track progress
    message("Processing q = ", q)
    
    # Perform dyna-LOCUS decomposition using original data
    result_ref = dynaLOCUS(Y = Y, q = q, V = V, n_subject = n_subject,
                           preprocess = preprocess, penalty_S = penalty_S, phi = phi, rho = rho, penalty_A = penalty_A, lambda = lambda,
                           maxIteration = maxIteration, speed_up = TRUE,
                           espli1 = espli1, espli2 = espli2, silent = silent, demean = demean)
    
    # extract the reference S matrix from the result
    S_ref = result_ref$S
    
    # initialize matrices to record results
    matching = matrix(NA, nrow = q, ncol = length(seeds))
    indexing = matrix(NA, nrow = q, ncol = length(seeds))
    averaging = matrix(NA, nrow = q, ncol = length(seeds))
    
    # loop over each seed
    for(seed in seeds){
      
      # set seed for reproducibility
      set.seed(seed)
      
      # sample subjects with replacement
      index_subjects = sample(1:n_subject, replace = TRUE)
      
      # construct new Y matrix by combining sampled subjects' data
      Yraw_subj_new = lapply(index_subjects, function(idx) Y[((idx - 1) * t_length + 1):(idx * t_length), ])
      Y_new = do.call(rbind, Yraw_subj_new)
      
      # perform dyna-LOCUS decomposition
      result = dynaLOCUS(Y = Y_new, q = q, V = V, n_subject = n_subject,
                         preprocess = preprocess, penalty_S = penalty_S, phi = phi, rho = rho, penalty_A = penalty_A, lambda = lambda,
                         maxIteration = maxIteration, speed_up = TRUE,
                         espli1 = espli1, espli2 = espli2, silent = silent, demean = demean)
      
      # extract the S matrix from the result
      S_new = result$S
      
      # calculate correlations between the original S and the new S matrix
      dat = data.frame(t(rbind(S_ref, S_new)))
      Cor = abs(cor(dat)[1:q, (q + 1):(2 * q)])
      
      # sort and record the highest correlations
      comb = matrix(0, nrow = q, ncol = 3)
      for (i in 1:q) {
        comb[i, ] = c(i, which.max(Cor[i, ]), max(Cor[i, ]))
      }
      
      # record results
      index = which(seeds == seed)
      matching[, index] = comb[, 3]
      indexing[, index] = comb[, 2]
      averaging[, index] = apply(Cor, 1, mean)
    }
    
    # calculate mean values
    avg_matching = apply(matching, 1, mean)
    avg_averaging = apply(averaging, 1, mean)
    
    # calculate reliability index and reproducibility
    RI[[paste0("q = ", q)]] = (avg_matching - avg_averaging) / (1 - avg_averaging)
    reliability[[paste0("q = ", q)]] = (avg_matching - avg_averaging) / avg_matching
  }
  
  return(list(RI = RI, 
              reliability = reliability))
}
