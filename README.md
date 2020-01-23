Comparing different matrix low-rank approximations with leverage scores

The dataset we tested our code on is Abalone Data set available from
https://archive.ics.uci.edu/ml/datasets/Abalone

compareAll generates the distance matrix from the data set and then compares to
Nystorm with uniform sampling
Nystrom with Exact Leverage scores
Nystrom with Frobenius Leverage scores
Nystrom with Spectral Leverage scores
Nystrom with Power method Leverage scores
for 20 different size of samplings. The desired rank of approximation is set to k=20. 