# MR

MR is a package that lets you extend standard dimensionality reduction algorithms such as ISOMAP, LLE, Laplacian EigenMap to handle the case when we have missing data.

The exact details about the algorithm can be found at https://arxiv.org/abs/1807.07610

The exmaples folder has interactive notebooks that steps you through some the major experiments from the paper. While doing so it also provides a useful demo for how to use this package. If you use this package please cite https://ieeexplore.ieee.org/document/8635955/metrics#metrics

## Usage:

`Isomap(X,k,d)` - Outputs the cooridnates of a d dimensional represntation of the data X. Here k is the Isomap paramter (i.e. will use k nearest neighbors). Data points should be represented as rows of the matrix X.

`mds(D,d;center=true,align=true)` - Outputs a dimensional representation for data points with distance given D. If center is true, then the data points have mean 0. If align is true, then it fixes a rotation of the data points. 

`procrustes(X,Y)` - Does a procrustes alginment of the data Y with the data X and returns the aligned data. Here each data point is a row.

`IOMR(D)` - returns the pertubration that needs to be added to D so that the metric has been fixed. Is an implementation of the IOMR-fixed algorithm. 

`MR_Missing(X,Q,k)` - Returns the adjacency matrix of the graph of the k nearest neighbors. Here X is the data, Q is a 0,1 matrix where a 0 entry represents the fact that the corresponding entry in the X matrix is missing. 

`apsp(A)` - Takes in the adjaceny matrix of graph and return a matrix with the ijth entry is the shorted disjance between points i and j
