This algorithm can learn Gaussian linear structural equation model with in-procedure and post-procedure outlier. 

paper : https://link.springer.com/article/10.1007/s42952-021-00160-2 (Robust estimation of Gaussian linear structural equation models with equal error variances)


## RGSEM_Algorithm

* data : n x p matrix or data.frame
* method : 'in.procedure', 'post.procedure' / Outlier type in graphical model
* alpha : a significant level for conditional independence tests
* C : a constant for Cook's distance 
* graph : a true graph matrix (e.g. 1 -> 2 mat[2,1] = 1) (If a true graph is entered as a input, estimated graph can be evaluated, otherwise it will not be evaluated.)
