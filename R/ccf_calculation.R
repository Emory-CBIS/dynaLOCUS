#' Calculate cross-correlation function (CCF) for connectivity traits
#'
#' This function calculates the cross-correlation function (CCF) between pairs of connectivity traits.
#'
#' @param A The trait loading matrix.
#' @param q The number of connectivity traits.
#' @param n_subject The number of subjects.
#'
#' @return An array containing the final CCF values for each pair of connectivity traits.
#'
#' @details The function follows these steps:
#' - For each subject, calculate the lag which achieves the largest absolute CCF between trait pairs.
#' - Determine the sign of the CCF for each trait pair.
#' - Recalculate the maximum value lag between connectivity trait pairs using the determined sign.
#' - Determine the mode of lag for each connectivity trait pair.
#' - Calculate and return the final CCF values.
#'
#' @export

ccf_calculation = function(A, q, n_subject){
  
  # function to calculate the mode
  get_mode = function(v){
    uniqv = unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  # dimension of A
  N = dim(A)[1]
  # time length
  t_length = N / n_subject
  
  # initialize arrays to record results
  ccf_lag_initial = array(NA, c(q, q, n_subject))
  ccf_lag_intermediate = array(NA, c(q, q, n_subject))
  ccf_value_final = array(NA, c(q, q, n_subject))
  
  # split A
  dat = matsplitter(A, t_length, q)
  # lag_max
  lag_max = t_length - 1
  
  # step 1: calculate initial CCF lags for each subject
  for(k in 1:n_subject){
    for(i in 1:(q - 1)){
      for(j in (i + 1):q){
        # calculate CCF and lags
        fit = ccf(dat[, i, k], dat[, j, k], ylab = "Cross-correlation", lag.max = lag_max, plot = FALSE)
        mat = data.frame(acf = fit$acf, lag = fit$lag)
        
        # identify the lag with the largest absolute CCF value
        index = which.max(abs(mat$acf))
        ccf_lag_initial[i, j, k] = mat$lag[index]
      }
    }
  }
  
  # step 2: determine the sign of CCF for each trait pair
  lag_sign = matrix(NA, nrow = q, ncol = q)
  for(i in 1:(q - 1)){
    for(j in (i + 1):q){
      prop = c(mean(ccf_lag_initial[i, j, ] < 0), mean(ccf_lag_initial[i, j, ] == 0), mean(ccf_lag_initial[i, j, ] > 0))
      index = which.max(prop)
      
      if(index == 1){
        lag_sign[i, j] = -1
      }else if (index == 2){
        lag_sign[i, j] = 0
      }else{
        lag_sign[i, j] = 1
      }
    }
  }
  
  # step 3: recalculate the max value lag using determined sign
  for(k in 1:n_subject){
    for(i in 1:(q - 1)){
      for(j in (i + 1):q){
        fit = ccf(dat[, i, k], dat[, j, k], ylab = "Cross-correlation", lag.max = lag_max, plot = FALSE)
        mat = data.frame(acf = fit$acf, lag = fit$lag)
        
        if(lag_sign[i, j] == -1){
          index = which.max(abs(mat[mat$lag < 0, ]$acf))
        }else if(lag_sign[i, j] == 1) {
          index = which.max(abs(mat[mat$lag > 0, ]$acf)) + t_length
        }else{
          index = t_length
        }
        
        ccf_lag_intermediate[i, j, k] = mat$lag[index]
      }
    }
  }
  
  # step 4: determine the mode of lag for each connectivity trait pair
  lag_mode = matrix(NA, nrow = q, ncol = q)
  for(i in 1:(q - 1)){
    for(j in (i + 1):q){
      lag_mode[i, j] = get_mode(ccf_lag_intermediate[i, j, ])
    }
  }
  
  # step 5: calculate the final CCF values based on the mode of lag
  for(k in 1:n_subject){
    for(i in 1:(q - 1)){
      for(j in (i + 1):q){
        fit = ccf(dat[, i, k], dat[, j, k], ylab = "Cross-correlation", lag.max = lag_max, plot = FALSE)
        mat = data.frame(acf = fit$acf, lag = fit$lag)
        ccf_value_final[i, j, k] = mat[mat$lag == lag_mode[i, j], ]$acf
      }
    }
  }
  
  # assign to an array
  ccf_value_array = ccf_value_final
  
  return(ccf_value_array)
}
