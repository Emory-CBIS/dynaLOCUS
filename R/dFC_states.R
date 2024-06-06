#' Retrieve Whole Brain dFC States Based on dyna-LOCUS Results
#'
#' This function retrieves whole-brain dynamic functional connectivity (dFC) states based on the results from the dyna-LOCUS method.
#'
#' @param A The loading matrix.
#' @param S The dyna-LOCUS result connectivity traits.
#' @param n_subject The total number of subjects in the study.
#' @param q The number of connectivity traits.
#' @param n_examplers The number of exemplars to select from each subject for clustering initialization. Default is 10.
#' @param gap The gap between selected exemplars in time points. Default is 5.
#' @param seed An integer for setting the random seed to ensure reproducibility. Default is 1.
#' @param n_centers The number of centers (or clusters) to use in k-means clustering.
#'
#' @return A matrix of reconstructed connectivity patterns based on the clustered loading matrix. Each row of the matrix represents the upper triangular part of a dFC state.

dFC_states = function(A, S, n_subject, q, n_examplers = 10, gap = 5, seed = 1, 
                      n_centers){
  
  # dimensions of A
  N = dim(A)[1]
  # time length per subject
  t_length = N / n_subject
  # split the loading matrix
  dat = matsplitter(A, t_length, q) 
  
  # initialize centers using exemplar loading series
  # select n_examplers from each subject
  index = seq(n_examplers, t_length, gap)
  A_sub = matrix(nrow = 0, ncol = q)
  for(i in 1:n_subject){
    A_sub = rbind(A_sub, dat[,,i][index,])
  }
  
  # # the following lines can be used to find the optimal number of clusters using the elbow method
  # set.seed(seed)
  # rslt_k_A = fviz_nbclust(A_sub, kmeans, method = "wss", iter.max = 50, k.max = 20)
  
  # perform k-means clustering on examplers to initialize centers
  rslt_A_sub = kmeans(A_sub, n_centers, iter.max = 50, nstart = 200)
  
  # perform k-means clustering on the loading matrix A based on initialized centers
  rslt_A = kmeans(A, centers = rslt_A_sub$centers, iter.max = 50, nstart = 200)
  
  # reconstruct the connectivity pattern
  centers = rslt_A$centers
  Y_construct = centers %*% S
  
  return(Y_construct = Y_construct)
}
