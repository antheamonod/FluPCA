# FluPCA

This repository provides the code and data to implement the real data example in [Monod et al.](https://arxiv.org/abs/1805.12400), which performs principal component analysis (PCA) on sets of phylogenetic trees arising from influenza data with respect to the BHV and tropical metrics.

PCA is a fundamental technique in descriptive and exploratory statistics that visualizes relationships within the data by reducing their dimensionality.  Classical PCA inherently relies on an assumption of a vector space, which does not apply to phylogenetic tree space.  As such, alternative interpretations are needed and have been proposed by [Nye et al.](https://academic.oup.com/biomet/article/104/4/901/4259146) with respect to the BHV metric (which we refer to as BHV PCA), and [Yoshida et al.](https://link.springer.com/article/10.1007/s11538-018-0493-4) with respect to the tropical metric (which we refer to as tropical PCA), to accommodate the non-Euclidean nature of phylogenetic tree space.

Tree PCA in our repository is implemented as a set of R routines following the methodology developed in the aforementioned references.  In particular, BHV PCA in our repository is implemented via code freely available at the [GeoPhytter+](http://www.mas.ncl.ac.uk/~ntmwn/geophytterplus/index.html) repository, which is written in [java](https://go.java/index.html?intcmp=gojava-banner-java-com).

For more detail on the underlying theory, please see [Monod et al.](https://arxiv.org/abs/1805.12400) for an overview (especially in reference to the data provided in this repository), or [Nye et al.](https://academic.oup.com/biomet/article/104/4/901/4259146) and [Yoshida et al.](https://link.springer.com/article/10.1007/s11538-018-0493-4) for full descriptions.

## The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics.  The most recent version of R can be downloaded from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).  CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  For specific details on how to compile, install, and manage R and R packages, refer to the manual [R Installation and Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).

## R Packages Required for Running FluPCA
Implementation of tree PCA on the provided influenza data requires the installation of the following R libraries:
* [ape](https://cran.r-project.org/web/packages/ape/index.html)
* [geophyttertools](https://github.com/grady/geophyttertools)
* [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html)

These packages may be installed either by typing the following command in an R shell:
```
install.packages("ape", dependencies=TRUE)
```
or by [installing R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### Parallel Computing
Given the computational intensity of the implementation on the data example, we use parallel computation.  In order to run this, the library [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) is also required.  This package is integrated within the R core, however must still be loaded by typing the following command in an R shell:
```
  require(parallel)
```
Currently, only a Windows version is available, however, users may modify the source code to create a Linux version accordingly.

## Data Availability and Preprocessing
The data provided in this repository are available in the `Data` directory are preprocessed in the following manner.  Genomic data for 1089 full length sequences of hemagglutinin (HA) for influenza A H3N2 from 1993 to 2017 in the state of New York were obtained from the [GI-SAID EpiFlu&trade;](https://www.gisaid.org/) database and aligned with [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) using default settings.  HA sequences from each season were related to those of the preceding season.  [Tree dimensionality reduction](https://arxiv.org/abs/1607.07503) using temporal windows of 5 consecutive seasons to create 21 datasets.  The date of each dataset corresponds to the first season (for example, the dataset dated 2013 consists of 5-leaved trees where the leaves come from seasosn 2013 through 2017).  Each unrooted tree in these datasets was constructed using the [neighbor-joining method](https://academic.oup.com/mbe/article/4/4/406/1029664) with Hamming distance.  Outliers were then removed from each season using [KDETREES](http://vps.fmvz.usp.br/CRAN/web/packages/kdetrees/vignettes/kdetrees.pdf).  On average, there are approximately 20'000 remaining trees in each dataset.  The trees are stored in [Newick format](https://en.wikipedia.org/wiki/Newick_format).  

## Description of R Scripts
The `Software` directory contains the following set of R scripts:
* `func_ssh.R` is a set of supplemental functions called in running BHV and tropical PCA, and computing the proportion of variance explained (R^2) for each PCA method

## Running Tree PCA and Reproducing Projective Visualizations

The `Figures` directory contains (edited) figures of the output of the BHV and tropical PCA routines applied to the 21 datasets provided in the `Data` directory.  Some of these figures are referenced in the main text.

## Questions and Feedback
For questions and comments on this repository, please contact [Anthea Monod](mailto:antheam@tauex.tau.ac.il) or [Qiwen Kang](mailto:qiwen.kang@uky.edu).  We appreciate any feedback you may have on our repository.
