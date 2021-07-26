Forward_Learning_fun_out = function(X, C = NULL, method = 'in.procedure') {
  
  require(MASS)
  X = as.matrix(X)
  p = ncol(X) ;   n = nrow(X)
  Ordering = rep(0, p) ; RemNodes = 1:p
  k = 1
  valid_list = list()
  score_valid_sample_temp = sapply(RemNodes,
                                   function(i) score_valid_sample_fun(X,i,predictor = NULL,C = C,valid_set = 1:n , method = method))
  Ordering[k] = which.min(unlist(score_valid_sample_temp[1,]))
  
  valid_list[[k]] =score_valid_sample_temp[2,][[Ordering[k] ]]
  RemNodes = setdiff(RemNodes, Ordering)
  
  while( length(RemNodes) > 1 ) {
    k = k + 1 ; Z = NULL
    
    Z = as.matrix(X[,Ordering[1:k-1]])
    score_valid_sample_temp = sapply(RemNodes,
                                     function(j) score_valid_sample_fun(X,j,Z,C= C,
                                                                        valid_set = valid_list[[k-1]], method = method))
    
    Ordering_j = which.min(unlist(score_valid_sample_temp[1,]))
    
    Ordering[k] = RemNodes[Ordering_j]
    
    valid_list[[k]] =score_valid_sample_temp[2,][[Ordering_j]]
    
    RemNodes = setdiff(RemNodes, Ordering)
  }
  Ordering[p] = RemNodes
  
  Z = as.matrix(X[,Ordering[1:p-1]])
  
  score_valid_sample_temp =score_valid_sample_fun(X , j = RemNodes , Z ,C ,valid_set = valid_list[[p-1]], method = method )
  
  valid_list[[p]] = score_valid_sample_temp[[2]]
  
  return( list(Ordering, valid_list ))
}


score_valid_sample_fun = function(X, j , predictor,C, valid_set , method) {
  n = nrow(X)
  
  X_idx = 1:n
  out_obs = NULL
  if(method == 'post.procedure') {
    X_idx = X_idx[valid_set] 
    out_obs = setdiff(1:n , valid_set)
  }
  
  
  if( is.null(predictor) ) {
    lm_fit = lm(X[X_idx,j] ~ 1) 
    # out_obs_temp = which(abs(rstudent(lm_fit)) >  qt(  1-out_alpha/n, n-2))
    out_obs_temp = which(abs(cooks.distance(lm_fit)) >  C/(length(X_idx)-1))
  }else {
    lm_fit = lm(X[X_idx,j] ~ predictor[X_idx,])
    # out_obs_temp = which(abs(rstudent(lm_fit)) >  qt(  1-out_alpha/n, n-ncol(predictor)-2))
    out_obs_temp = which(abs(cooks.distance(lm_fit)) >  C/(length(X_idx)-ncol(predictor)-1))
  }
  
  while(length(out_obs_temp) != 0)
  {
    out_obs_temp = X_idx[out_obs_temp]
    out_obs = union(out_obs, out_obs_temp)
    X_idx = setdiff(X_idx, out_obs_temp)
    
    if(is.null(predictor) ) {
      lm_fit = lm(X[X_idx,j] ~ 1) 
      # out_obs_temp = which(abs(rstudent(lm_fit)) >  qt(  1-out_alpha/(n-length(out_obs)), n-length(out_obs)-2))
      out_obs_temp = which(abs(cooks.distance(lm_fit)) >  C/(length(X_idx)-1))
    }else {
      lm_fit = lm(X[X_idx,j] ~ predictor[X_idx,])
      # out_obs_temp = which(abs(rstudent(lm_fit)) >  qt(  1-out_alpha/(n-length(out_obs)), n-length(out_obs)-ncol(predictor)-2))
      out_obs_temp = which(abs(cooks.distance(lm_fit)) >  C/(length(X_idx)-ncol(predictor)-1))
    }
    
  }
  if( is.null(predictor) ){ 
    return(list(var(X[X_idx,j]), X_idx))
  }else{
    return(list(sum( lm_fit$resid^2 )/length(X[X_idx,j]), X_idx))
  }
  
  
}


