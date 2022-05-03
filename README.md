# KFRBcomp
Comparison between Kalman Filter and Recursive Bayes methods for estimating Random walks in alpha parameter of the Ricker stock recruitment curve. 





These comparisons include estimates obtained with the Kalman filter method implemented within the R package dlm (the estimates are filtered and smoothed).
The recursive Bayes model is coded in TMB. The model is designed to be run with MCMC sampling, but MLE estimates are also produced considerinthe vector of alpha parameters as a random effect.


Comparison is based on empirical data compiled by Dan Greenberg

TODO:

-	Check the model fits for convergence. 
-	Estimate the recursive Bayes with the MCMC algorithms
-	Compare the confidence intervals of the estimates – I expect that the Kalman filter might give me narrower Confidence intervals, because it is an analytical solution (I think?) but I might be wrong.
-	I want to add to the comparison the Kalman filter algorithm coded by Carrie Holt . As I think that he methodology is more in line with the models commonly applied in fisheries. - This seems to fail to converge more often than not. 
