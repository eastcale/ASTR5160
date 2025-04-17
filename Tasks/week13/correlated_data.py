#!/usr/bin/env python3

import numpy as np

#CTE Taking the transpose of the data this time,
#CTE so that each row is a list of the y values
#CTE as this is what numpy expects with np.cov.
#CTE This is also the only way to get a 10x10 cov. matrix,
#CTE as it should be, since there are only 10 bins.
dat = np.loadtxt('/d/scratch/ASTR5160/week13/line.data').T

#CTE Covariance Matrix
cov = np.cov(dat)

#CTE Diagonal of the covariance (should be equal to the variances)
cov_diag = np.diag(cov)
var = np.array([np.var(i, ddof = 1) for i in dat])

#CTE As expected, all the diagonal terms give the variance for each bin
diffs = cov_diag - var

#CTE Finding the correlation matrix
corr = np.corrcoef(dat)

#CTE np.where() gives an ordered pair of
#CTE arrays with the first array corresponding to
#CTE row indices and the second corresponding to 
#CTE column indices.
#CTE The result of the following is (array([0]), array([3]))
#CTE which means the min is in row 0 and column 3
min_val = np.min([np.min(i) for i in corr])
min_in = np.where(corr == min_val)

#CTE Getting rid of the diagonal elements is tricky
#CTE Not sure if there is a more elegant solution to this.
max_val = np.max([np.max(np.delete(corr[i], i)) for i in range(len(corr))])
max_in = np.where(corr == max_val)

print(f'The most correlated columns are x{max_in[0]} and x{max_in[1]} with a correlation coefficient of {max_val}')
print(f'The most anti-correlated columns are x{min_in[0]} and x{min_in[1]} with a correlation coefficient of {min_val}')