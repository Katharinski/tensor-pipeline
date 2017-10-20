# tensor-pipeline
Software to build and decompose tensors in Matlab from fMRI and simulated data, as described in Glomb, Katharina, et al. "Resting state networks in empirical and simulated dynamic functional connectivity." NeuroImage 159 (2017): 388-402.

Software developed by Katharina Glomb

The pipeline is intended for building and decomposing 3-way-tensors representing dynamic functional connectiviy. Multiple subjects are needed to apply the parameter selection presented here. The TR used is 2 seconds (can be adapted).

All steps can be seen in main_script.m. Summary with function names:

0) load your data, choose the FC measure (see point 2) and, if sliding windows used (correlation-based measures, mutual information), choose a window width (in frames)

optional: simulate data - sim_script, Get_balanced_weights, DMF_sim - as in Deco et al. (2014)

1) build tensors from time series (make_tensor)
This method is intended for fMRI data with low spatial resolution which allows comparison with node and edge-style models. Tensors consist of adjacency matrices - in particular, functional connectivity (FC) matrices. FC is computed using 1) Pearson correlation ('corr'), 2) absolute values of Pearson correlation ('abs'), 3) Mutual Information as introduced in Kraskov et al. (2004) ('MI'), 4) phase differences (phase differences between 0 and pi are mapped to adjacency values between 1 and 0, respectively). This latter measure is the only instantaneous one and was not used in the reference mentioned above.

surrogates_cov: for creating phase-shuffled (stationary) surrogate time series; decomposition results will be compared to real data

make_tensor: 
  for methods 'abs' and 'corr': calls prepdata_fcms
  for method 'MI': calls MI_kraskov 
  for method 'ph': calls get_phases (-->filter_fMRI), phase_diffs_adj

2) decompose tensors (decomp_tens) - requires third party software (freely available)
Tensor decomposition assumes that tensors are described by a small number of features that superimpose linearly (see Cichocki et al., 2009). We use two different algorithms to obtain these features: CPD (canonical polyadic decomposition, also known as PARAFAC) for correlation, and NCP (non-negative CP) for the other measures, because they are nonnegative. Algorithms used here are described in Phan et al. (2013) and Kim and Park (2012), respectively. References and links are given in decomp_tens.m.
Tensors are composed using different numbers of features (F), and better results are achieved when tensors are binarized, keeping only the biggest FC pairs (different thresholds). As a result, one decomposition is run for each value of F and threshold value and resulting features as well as decomposition errors are saved. For binarized tensors, the Hamming distance is used to compute these errors (get_error_hamming.m).

3) cluster resulting features to get templates - cluster_spfeats, make_templates
In order to evaluate the quality of the decompositions, we K-means-cluster extracted spatial features across subjects and use the silhouette value to quantify how well the features are clustered. This means that our criterion is that the features generalize across subjects. For each value of F, K, and threshold, clustering quality is assessed taking the difference between surrogate and real data. The result is a set of "prototypical" features, or "templates".

optional: match templates from real data with features extracted from simulated/surrogate/single subject data - match_templates
Vectors are quantized and overlap is computed using confusion matrices and Cohen's kappa (see Glomb et al. 2017 for details).

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

Simulations: 
Deco G, et al. "How local excitation-inhibition ratio impacts the whole brain dynamics." JNeurosci 34 (2014):7886–7898.
