# Non-Uniform Sampling of Fixed Margin Binary Matrices 

This is the accompanying code to the paper submitted to 2020 Foundations of Data Science Conference in Seattle Wa, Oct 18-20, 2020.

## Abstract

Data sets in the form of binary matrices are ubiquitous across scientific domains, and researchers are often interested in identifying and quantifying noteworthy structure. 
One approach is to compare the observed data to that which might be obtained under a null model.
Here we consider sampling from the space of binary matrices which satisfy a set of marginal row and column sums.
Whereas existing sampling methods have focused on uniform sampling from this space, we introduce modified versions of two elementwise swapping algorithms which sample according to a non-uniform probability distribution defined by a weight matrix, which gives the relative probability of a one for each entry.
We demonstrate that values of zero in the weight matrix, i.e. _structural zeros_ are generally problematic for swapping algorithms, except when they have special monotonic structure.
We explore the properties of our algorithms through simulation studies, and illustrate the potential impact of employing a non-uniform null model using a classic bird habitation dataset.

Here is the link to the [pre-published paper]()
