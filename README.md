This algorithm can learn Gaussian linear structural equation model with in-procedure and post-procedure outlier. 

paper : https://www.springer.com/journal/42952 (Robust estimation of Gaussian linear structural equation models with equal error variances)


## RGSEM_Algorithm

* data : n x p matrix or data.frame
* method : 'in.procedure', 'post.procedure' / Outlier type in graphical model
* alpha : a significant level for conditional independence tests
* C : a constant for Cook's distance 
