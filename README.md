# FluPCA

This repository provides the code and data to implement the real data example in Monod et al., which performs principal component analysis (PCA) on sets of phylogenetic trees arising from influenza data with respect to the BHV and tropical metrics.

PCA is a fundamental technique in descriptive and exploratory statistics that visualizes relationships within the data by reducing their dimensionality.  Classical PCA inherently relies on an assumption of a vector space, which does not apply to phylogenetic tree space.  As such, alternative interpretations are needed and have been proposed by [Nye et al.](https://academic.oup.com/biomet/article/104/4/901/4259146) with respect to the BHV metric (which we refer to as BHV PCA), and [Yoshida et al.](https://link.springer.com/article/10.1007/s11538-018-0493-4) with respect to the tropical metric (which we refer to as tropical PCA), to accommodate the non-Euclidean nature of phylogenetic tree space.

Tree PCA in our repository is implemented as a set of R routines following the methodology developed in the aforementioned references.  In particular, BHV PCA in our repository is implemented via code freely available at the [GeoPhytter+](http://www.mas.ncl.ac.uk/~ntmwn/geophytterplus/index.html) repository, which is written in java.

### The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics.  The most recent version of R can be downloaded from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).  CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  For specific details on how to compile, install, and manage R and R packages, refer to the manual [R Installation and Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).

### R Packages Required for Running FluPCA
Implementation of tree PCA on the provided influenza data requires the installation of the following R libraries:
*[ape](https://cran.r-project.org/web/packages/ape/index.html)
