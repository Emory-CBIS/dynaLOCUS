#' Iterative Update Algorithm for dyna-LOCUS
#'
#' @param Y The (preprocessed) concatenated group-level connectivity data.
#' @param A A loading matrix from the previous iteration.
#' @param theta A list containing trait-specific parameters and channels.
#' @param q The number of latent connectivity traits to extract.
#' @param K The number of edges.
#' @param V The number of nodes.
#' @param n_subject The total number of subjects in the study.
#' @param t_length The number of sliding windows.
#' @param penalty_S The option for the penalization function for the sparsity regularization for the connectivity traits. 
#' Users can choose "NULL", "L1", or "SCAD". Defaults to "L1". 
#' @param phi A tuning parameter for the element-wise penalty on S. Defaults to 2.
#' @param penalty_A The option for the temporal smoothness regularization for the trait loadings. 
#' If \code{TRUE}, a temporal smoothness penalty is included in the optimization function to encourage similarity in the trait loadings 
#' in adjacent time windows. Defaults to \code{TRUE}.
#' @param lambda A tuning parameter for the temporal smooth penalty on \code{A}.
#' @param H_star The dewhitening matrix.
#' @param list_sparse A list to record the sparseness of connectivity traits.
#' @param preprocess If \code{TRUE}, the concatenated group-level connectivity data is preprocessed.
#' @param speed_up If \code{TRUE}, the updating process is sped up using singular value decomposition.
#' @param silent If \code{TRUE}, suppresses the printout of penalties applied on \code{A} and \code{S}.
#' 
#' @import MASS

dyna_update = function(Y, A, S, theta,
                       q, K, V, n_subject, t_length, 
                       penalty_S, phi, penalty_A, lambda, H_star, list_sparse,
                       preprocess, speed_up, silent){
  
  # preparations to update A
  if(t_length > 1){
    # construct R
    R = matrix(0, nrow = t_length-1, ncol = t_length)
    for(t in 1:(t_length-1)){
      R[t,t] = 1
      R[t, t+1] = -1
    }
    # construct R_star
    R_star = kronecker(diag(1, n_subject), R)
    if(preprocess){
      # construct W = R_star*H_star
      W = R_star %*% H_star
      # calculate eigen values and vectors of W'W
      eigenW = eigen(t(W) %*% W, T)
      # eigen vectors
      Q2 = eigenW$vectors
      # eigen values
      Lmbd2 = eigenW$values
    }else{
      # calculate eigen values and vectors of R_star'R_star
      eigenR_star = eigen(t(R_star) %*% R_star, T)
      # eigen vectors
      Q2 = eigenR_star$vectors
      # eigen values
      Lmbd2 = eigenR_star$values
    }
  }
  
  # node to node edge list 
  rmat = matrix(rep(1:V, V), ncol=V)
  cmat = t(rmat)
  Lcoorsym = matrix(0, ncol=2, nrow=V*(V-1)/2)
  Lcoorsym[,1] = Ltrans(rmat, F)
  Lcoorsym[,2] = Ltrans(cmat, F)
  # use Rls to store the value of Rl in each latent connectivity traitss
  Rls = numeric(q)
  
  # use theta to store source-specific parameters and channels
  theta_new = list()
  
  # for each ic, conditioning on others, estimate latent channels 
  newS = array(dim = c(q, K))
  # update Xl and Dl's
  for(curr_ic in 1:q){

    theta_new[[curr_ic]] = list()
    
    # current latent connectivity traits related theta
    theta_ic = theta[[curr_ic]]
    
    # lth row of A_tilde transpose (1*q) = lth column of A_tilde (1*q)
    S_lold = t(theta_ic$M_l %*% Y)
    
    # penalize S or Yic
    if(is.null(penalty_S)){
      S_new_0 = S_lold
    }else if(penalty_S == "SCAD") {
      S_new_0 = SCAD_func(S_lold, phi = phi, gamma = 2.1)
      S_new_0 = S_new_0 / sd(S_new_0) * sd(S_lold)
    }else if(penalty_S == "L1") {
      S_new_0 = sign(S_lold) * (abs(S_lold) - phi) * (abs(S_lold) >= phi)
      S_new_0 = S_new_0 / sd(S_new_0) * sd(S_lold)
    }else{
      print("No penalty available!")
      stop()
    }
    
    if(mean(abs(S[curr_ic,]) < 0.1) > 0.75){
      list_sparse = unique(c(list_sparse, curr_ic))
    }
  
    if(speed_up){
      # dimension of Rl in the current latent connectivity traits
      Rls[curr_ic] = dim(theta_ic$X_l)[1]
      Sl = Ltrinv(S_new_0, V, F) + diag(diag(t(theta_ic$X_l)%*%diag(theta_ic$lam_l)%*%theta_ic$X_l))
      
      # update and assign
      eigenSl = eigen(Sl)
      orderEigen = order(abs(eigenSl$values), decreasing = T)
      Rl = Rls[curr_ic]
      eigenset = orderEigen[1:Rl]
      
      theta_new[[curr_ic]]$lam_l = eigenSl$values[eigenset]
      if(theta_new[[curr_ic]]$lam_l[1] < 0){
        theta_new[[curr_ic]]$lam_l = -1 * theta_new[[curr_ic]]$lam_l
      }
      for(j in 1:Rl){
        theta_ic$X_l[j,]= eigenSl$vectors[,eigenset[j]]
      }
      theta_new[[curr_ic]]$X_l = theta_ic$X_l
    }else{
      
      if(curr_ic %in% list_sparse){
        theta_new[[curr_ic]] = theta[[curr_ic]]
      }else{
        # D inverse (of length Rl)
        Dinverse = 1/theta_ic$lam_l 
        
        v = 1 
        while(v <= V){
          # Xl with the vth column removed -> transpose -> (V-1)*Rl
          Hlv = t(theta_ic$X_l[,-v])  
          # penalized Yic(p*1 vector) with elements relate to node v -> (V-1)*1
          yvpen = S_new_0[which(Lcoorsym[,1] == v | Lcoorsym[, 2] == v),]              
          # inverse of t(Hlv)*Hlv
          Sigmalv = MASS::ginv( t(Hlv)%*%Hlv )
          
          if(sd(yvpen) == 0){
            theta[[curr_ic]]$X_l[,v] = 0
          }else{
            theta[[curr_ic]]$X_l[,v] = Dinverse*(Sigmalv%*%t(Hlv)%*%yvpen)
          }
          v = v+1
        }
        # if Xl is NA then set it to be 0
        theta[[curr_ic]]$X_l[is.na(theta[[curr_ic]]$X_l)] = 0
        
        theta_new[[curr_ic]]$X_l = theta[[curr_ic]]$X_l
        
        # calculate Zl with rth column being L(Xl%*%t(Xl)) in the paper
        Xstarstack = t(apply(theta[[curr_ic]]$X_l, 1, function(x){ 
          x = matrix(x, ncol=1)
          return(Ltrans(x%*%t(x),F))
        }))
        theta_new[[curr_ic]]$lam_l =  as.numeric(MASS::ginv(Xstarstack %*% t(Xstarstack)) %*% Xstarstack %*% S_new_0)  
        
        # scale each column of Xl 
        Rl = dim(theta_new[[curr_ic]]$X_l)[1]
        for(r in 1:Rl){
          scale = sqrt(sum(theta_new[[curr_ic]]$X_l[r,]^2))
          theta_new[[curr_ic]]$X_l[r,] = theta_new[[curr_ic]]$X_l[r,]/scale
          theta_new[[curr_ic]]$lam_l[r] = theta_new[[curr_ic]]$lam_l[r]*scale
        }
      }
    }
    
    # calculate Sl based on Xl and Dl
    newS[curr_ic,] = Ltrans(t(theta_new[[curr_ic]]$X_l)%*% 
                              diag(theta_new[[curr_ic]]$lam_l) %*%
                              theta_new[[curr_ic]]$X_l,F)  
  
  }
  
  # update A 
  if((!penalty_A)|(t_length == 1)){
    # if no penalty or length = 1
    newA = Y %*% t(newS) %*% MASS::ginv(newS %*% t(newS))
  }else{
    # calculate eigen values and vectors of SS'
    eigenS = eigen(newS %*% t(newS),T)
    # eigen vectors
    Q1 = eigenS$vectors
    # eigen values
    Lmbd1 = eigenS$values
    # construct D
    D = t(Q2) %*% Y %*% t(newS) %*% Q1 
    
    # calculate a scaler for maintaining lambda at 0-10 level
    n_digits = nchar(as.integer(floor(max(abs(D)))))
    scaler = 10^n_digits
    
    # calculate A_star
    A_star = matrix(nrow = dim(D)[1], ncol = dim(D)[2])
    for(i in 1:(dim(D)[1])){
      for(j in 1:(dim(D)[2])){
        A_star[i, j] = D[i,j]/(lambda*scaler*Lmbd2[i] + Lmbd1[j])
      }
    }
    # calculate A_tilde or newA
    newA = Q2 %*% A_star %*% t(Q1)
  }
  
  for(i in 1:q){
    ai = sd(newA[,i])
    theta_new[[i]]$lam_l = theta_new[[i]]$lam_l * ai
    newA[,i] = newA[,i] / ai
    newS[i,] = newS[i,] * ai
  }
  
  # save M_l, X_l into theta_new
  # g-inverse of newA
  Mnew = MASS::ginv(t(newA) %*% newA) %*% t(newA) 
  for(l in 1:q){
    theta_new[[l]]$M_l = Mnew[l,]
  }
  return(list(A = newA, theta = theta_new, S = newS, list_sparse = list_sparse))
}
