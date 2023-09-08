# scPMP
Recent advances in single-cell technologies have enabled high-resolution characterization of tissue and cancer compositions. 
Although numerous tools for dimension reduction and clustering are available for single-cell data analyses, these methods often fail to simultaneously preserve local cluster structure and global data geometry. 

To address those challenges, we developed a novel analyses framework, **S**ingle-**C**ell **P**ath **M**etrics **P**rofiling (scPMP), which leverages the usefulness of $p$-powered path metrics for dimension reduction and clustering. 

Distances between cells are measured in a data-driven way which is both density sensitive (decreasing distances across high density regions) and respects the underlying data geometry. 

By combining path metrics with multidimensional scaling, a low dimensional embedding of the data is obtained which respects both the global geometry of the data and preserves cluster structure.

For the use of **scPMP**, please refer to our [R_tutorial](https://github.com/andrianamanousidaki/scPMP/tree/master/R_tutorial).

For the reproduction of our results:
1. Access and pre-process datasets by running the script [R_functions/Preparation_of_Data_Sets.R](https://github.com/andrianamanousidaki/scPMP/blob/master/R_functions/Preparation_of_Data_Sets.R) (Count_data folder contains only simulated manifold data sets)
2. Reproduce the results of a method for a particular dataset with the function [R_functions/clustering.R](https://github.com/andrianamanousidaki/scPMP/blob/master/R_functions/clustering.R)

