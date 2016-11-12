# tensor-pipeline
Software to build and decompose tensors in Matlab from fMRI data

Software developed by Katharina Glomb

The pipeline is intended for building and decomposing 3-way-tensors and comparing resulting features to templates. Multiple subjects are needed to apply the parameter selection presented here.

Steps:
1)build tensors from time series
This method is intended for fMRI data with low spatial resolution which allows comparison with node and edge-style models. Tensors consist of adjacency matrices - in particular, functional connectivity (FC) matrices - which are computed for sliding windows. FC is computed using 1) Pearson correlation, 2) Mutual Information as introduced in Kraskov et al. (2004).

2)decompose tensors 
Tensor decomposition assumes that tensors are described by a small number of features that superimpose linearly (see Cichocki et al., 2009). We use two different algorithms to obtain these features: CPD (canonical polyadic decomposition, also known as PARAFAC) for correlation, and NTF (non-negative tensor factorization) for MI. Algorithms used here are described in Phan et al. (2013) and Kim and Park (2012), respectively.

3)estimate parameters for decomposition
It is impossible to analytically estimate the number of features that should be used in the decomposition. Furthermore, we found that reducing noise by thresholding tensors, removing all but the very biggest FC values, yields better results. In order to evaluate the goodness of the decompositions, we cluster extracted features across subjects and use the silhouette value to quantify how well the features are clustered. This means that our criterion is that the features generalize across subjects. Reconstruction fit is also provided. Importantly, these criteria have to be compared to surrogate data (see Hindriks et al. (2016) for an example of correct usage with fMRI data). The result is a set of "prototypical" features, or "templates".

The pipeline makes use of third party software: Tensor toolbox (Bader & Kolda, 2012), NTF toolbox (Kim & Park, 2008), TENSORBOX (Phan et al., 2013)

References: 
Tensor decomp in general:
Cichocki, Andrzej, et al. Nonnegative matrix and tensor factorizations: applications to exploratory multi-way data analysis and blind source separation. John Wiley & Sons, 2009.

MI:
Kraskov, Alexander, Harald Stögbauer, and Peter Grassberger. "Estimating mutual information." Physical review E 69.6 (2004): 066138.

Decomp algorithms:
Kim, Hyunsoo, and Haesun Park. "Nonnegative matrix factorization based on alternating nonnegativity constrained least squares and active set method." SIAM journal on matrix analysis and applications 30.2 (2008): 713-730.

Phan, Anh-Huy, Petr Tichavsky, and Andrzej Cichocki. "Low complexity damped Gauss--Newton algorithms for CANDECOMP/PARAFAC." SIAM Journal on Matrix Analysis and Applications 34.1 (2013): 126-147.

Surrogate data:
Hindriks, R., et al. "Can sliding-window correlations reveal dynamic functional connectivity in resting-state fMRI?." NeuroImage 127 (2016): 242-256.

Third party software
Tensor toolbox:
Bader, Brett W., and Tamara G. Kolda. "Matlab tensor toolbox version 2.5." Available online, January 7 (2012).

Anh-Huy Phan, Petr Tichavsky and Andrzej Cichocki, "TENSORBOX: a Matlab package for tensor decomposition", 2013, available online at http://www.bsp.brain.riken.jp/~phan/tensorbox.php
