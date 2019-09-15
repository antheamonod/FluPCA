# FluPCA

This repository provides the code and data to implement the real data example in Monod et al., which performs principal component analysis (PCA) on sets of phylogenetic trees with respect to the tropical and BHV metrics.

PCA is a fundamental technique in descriptive and exploratory statistics that visualizes relationships within the data by reducing their dimensionality.  Classical PCA inherently relies on an assumption of a vector space, which does not apply to phylogenetic tree space.  As such, alternative interpretations are needed and have been proposed by [Nye et al.](https://academic.oup.com/biomet/article/104/4/901/4259146) and [Yoshida et al.](https://link.springer.com/article/10.1007/s11538-018-0493-4) to accommodate the non-Euclidean nature of phylogenetic tree space.

Tree PCA is implemented as a set of R routines following the methodology developed in the aforementioned references.

### The R Environment

R is a widely used, free, and open source software environment for statistical computing and graphics.  The most recent version of R can be downloaded from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).  CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).
