rmd-runs
========

The RMD_RUNS function will log10 transform peptide peak intensity, that is, peptide abundance data and determine if any LC-MS analyses (ie,runs) in a N x P peptide data set are statistical outliers. The statistical analysis is based on summarizing each LC-MS run as a vector of q=5 summary statistics which describe the peptide abundance distribution for a specific run; a N x q matrix is then analyzed using robust PCA to compute a robust estimate of the covariance matrix used in a the calculation of a robut Mahalanobis distance.

Croux, C. and Ruiz-Gazen, A (2005), "High breakdown estimatros for principal componets: the Projection-pursuit approach revisited", J Multivariate Analysis,95,206-226.
