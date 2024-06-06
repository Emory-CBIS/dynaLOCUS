#' Calculate energy and variation for connectivity traits
#'
#' This function calculates the energy and variation for connectivity traits. 
#'
#' @param A The loading matrix.
#' @param q The number of connectivity traits.
#' @param t_length The length of the time series per subject.
#'
#' @return A list containing:
#' \item{trait_energy}{A vector representing the mean log-transformed energy of each connectivity trait across subjects.}
#' \item{trait_variation}{A vector representing the mean variation of each connectivity trait across subjects.}
#'
#' @export

energy_var_cal = function(A, q, t_length){
  
  # total number of subjects
  N = dim(A)[1] / t_length
  
  # initialize matrices to record energy and variation
  energy_mat = matrix(nrow = N, ncol = q)
  var_mat = matrix(nrow = N, ncol = q)
  
  # calculate energy and variation for each subject and connectivity trait
  for(i in 1:N){
    for(j in 1:q){
      begin_index = (i - 1) * t_length
      # calculate energy
      energy_mat[i, j] = sum(A[(begin_index + 1):(begin_index + t_length), j]^2)
      # calculate variation
      vec = A[(begin_index + 1):(begin_index + t_length), j]
      var_mat[i, j] = mean(abs(diff(vec, lag = 1)))
    }
  }
  
  # log-transform energy due to high skewness and calculate the mean across subjects
  trait_energy = apply(log(energy_mat), MARGIN = 2, FUN = mean)
  
  # calculate the mean standard deviation across subjects
  trait_variation = apply(var_mat, MARGIN = 2, FUN = mean)
  
  return(list(trait_energy = trait_energy, 
              trait_variation = trait_variation))
}
