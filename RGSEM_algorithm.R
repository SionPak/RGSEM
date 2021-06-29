
GSEM_Algorithm = function(data, method, alpha = 0.001,  direction ="forward", graph = NULL,
                          max_degree = 1 , C=NULL){
  set.seed(1)
  library(bnlearn)
  ###################
  X = as.matrix( data )
  p = ncol(X)
  n = nrow(X)
  RemNode = 1:p 
  pi_GSEM = NULL
  Estimated_G = Estimated_O = matrix(0, p ,p)
  evaluation_result_GSEM = evaluation_result_GSEM_MEC = evaluation_result_GSEM_Oracle = NULL
  ####
  Runtime = proc.time()[3]
  
  #### Step 1): Finding the Ordering ####
  

  result = Forward_Learning_fun_out(X, max_degree = max_degree,C, method = method)
  Ordering = result[[1]]
  valid_obs = result[[2]]
  
  #### Step 2): Finding the Parents ####
  
  used_ci_test = "zf"
  #used_ci_test = "smc-zf"
  #used_ci_test = "mi-g"
  #used_ci_test = "smc-mi-g"
  #used_ci_test = "cor"
  
  ####### Cook's
  for(m in 2:p){
    j = Ordering[m]
    for(k in Ordering[1:(m-1)]){
      if(m > 2){
        S = setdiff( Ordering[ 1:(m-1)], k )
        
        valid_idx = Reduce(intersect,valid_obs[Ordering[1:m]])
        
        parent_pvalue = ci.test(X[valid_idx, j], X[valid_idx, k], X[valid_idx, S], test = used_ci_test)$p.value
      }else{
        
        valid_idx = Reduce(intersect,valid_obs[Ordering[1:m]])
        parent_pvalue = ci.test(X[valid_idx, j], X[valid_idx, k], test = used_ci_test)$p.value
      }
      if(parent_pvalue < alpha){
        Estimated_G[j, k] = 1
      }
    }
  } 
  
  
  ####
  Runtime = proc.time()[3] - Runtime
  
  #print(paste("It takes: ", Runtime))
  
  Estimated_G = matrix(Estimated_G, ncol = p)
  
  ####
  est_MEC = dag2cpdagAdj(Estimated_G)
  
  ####
  if( !is.null(graph) ){
    B = graph
    B[B!=0] =1
    MEC = B + t(B)
    Oracle_DAG = estimated_graph_fun( MEC, Ordering)
    evaluation_result_GSEM = evaluation_fun( B, Estimated_G ) 
    evaluation_result_GSEM_MEC = evaluation_fun( dag2cpdagAdj(B), est_MEC)
    evaluation_result_GSEM_Oracle = evaluation_fun( B, Oracle_DAG ) 
  }
  
  return(
    list( DAG_Evaluation = evaluation_result_GSEM, 
          MEC_Evaluation = evaluation_result_GSEM_MEC,  
          Oracle_Evaluation = evaluation_result_GSEM_Oracle,
          DAG = Estimated_G, 
          Ordering = Ordering, 
          Time = Runtime)
  )
}

