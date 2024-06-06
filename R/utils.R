#' Helper function to transform the upper triangular of a matrix to a vector
#' 
#' @param X A symmetric matrix
#' @param d If \code{True}, include the diagonal part of the matrix X. Defaults to \code{FALSE}.
#' 
#' @export

Ltrans = function(X, d = F){ 
  X[upper.tri(X, d)]  
} 

#' Helper function to transform a vector to a matrix
#' 
#' @param X A vector representing the upper triangular of a matrix
#' @param V The dimension of the target matrix
#' @param d If \code{TRUE}, the vector X include the diagonal part of the matrix. Defaults to \code{FALSE}.
#' 
#' @export

Ltrinv = function(x, V, d = F){ 
  Y = matrix(0, ncol = V, nrow = V)
  Y[upper.tri(Y, d)] = x
  
  return(Y + t(Y) - d*diag(diag(Y)))  
}


#' Helper function to calculate the scad penalty
#' 
#' @param yvpen A vector upon penalization
#' @param phi A tunning parameter
#' @param gamma A tunning parameter

SCAD_func = function(yvpen, phi = 2,gamma = 2.1){
  if(gamma <= 2){gamma = 2.01; print("Gamma needs > 2!!!")}
  ynew = sign(yvpen)*(abs(yvpen)-phi)*(abs(yvpen)>=phi)*(abs(yvpen)<=2*phi)+ 
    yvpen*(abs(yvpen) > gamma*phi) + 
    ((gamma-1)*yvpen-sign(yvpen)*gamma*phi)/(gamma-2)*(abs(yvpen)<=gamma*phi)*(abs(yvpen)>2*phi)
  if(sd(ynew) < 0.0000001){return(yvpen)}
  return(ynew)
}

#' Split a matrix into submatrices
#'
#' This function splits a given matrix \code{M} into smaller submatrices of dimensions \code{r} by \code{c}.
#'
#' @param M A matrix to be split.
#' @param r The number of rows in each submatrix.
#' @param c The number of columns in each submatrix.
#'
#' @return An array where each slice along the third dimension represents a submatrix of dimensions \code{r} by \code{c}.

matsplitter = function(M, r, c) {
  rg = (row(M) - 1) %/% r + 1
  cg = (col(M) - 1) %/% c + 1
  rci = (rg - 1) * max(cg) + cg
  N = prod(dim(M)) / (r * c)
  cv = unlist(lapply(1:N, function(x) M[rci == x]))
  dim(cv) = c(r, c, N)
  cv
}


