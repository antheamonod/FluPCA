# FluPCA

This repository provides software and data to implement the real data example in [Monod et al.](https://arxiv.org/abs/1805.12400), which performs principal component analysis (PCA) on sets of phylogenetic trees arising from influenza data with respect to the BHV and tropical metrics.

PCA is a fundamental technique in descriptive and exploratory statistics that visualizes relationships within the data by reducing their dimensionality.  Classical PCA inherently relies on an assumption of a vector space, which does not apply to phylogenetic tree space.  As such, alternative interpretations are needed and have been proposed by [Nye et al.](https://academic.oup.com/biomet/article/104/4/901/4259146) with respect to the BHV metric (which we refer to as BHV PCA), and [Yoshida et al.](https://link.springer.com/article/10.1007/s11538-018-0493-4) with respect to the tropical metric (which we refer to as tropical PCA), to accommodate the non-Euclidean nature of phylogenetic tree space.

Tree PCA in our repository is implemented as a set of R routines following the methodology developed in the aforementioned references.  In particular, BHV PCA in our repository is implemented via code freely available at the [GeoPhytter+](http://www.mas.ncl.ac.uk/~ntmwn/geophytterplus/index.html) repository, which is written in [Java](https://go.java/index.html?intcmp=gojava-banner-java-com).  Tropical PCA is implemented by code freely available in the [tropPCA](https://github.com/QiwenKang/tropPCA) repository, which is written in R.

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
The data provided in this repository are available in the `Data` directory are preprocessed and cleaned in the following manner.  Genomic data for 1089 full length sequences of hemagglutinin (HA) for influenza A H3N2 from 1993 to 2017 in the state of New York were obtained from the [GI-SAID EpiFlu&trade;](https://www.gisaid.org/) database and aligned with [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) using default settings.  HA sequences from each season were related to those of the preceding season.  [Tree dimensionality reduction](https://arxiv.org/abs/1607.07503) using temporal windows of 5 consecutive seasons to create 21 datasets.  The date of each dataset corresponds to the first season (for example, the dataset dated 2013 consists of 5-leaved trees where the leaves come from seasosn 2013 through 2017).  Each unrooted tree in these datasets was constructed using the [neighbor-joining method](https://academic.oup.com/mbe/article/4/4/406/1029664) with Hamming distance.  Outliers were then removed from each season using [KDETREES](http://vps.fmvz.usp.br/CRAN/web/packages/kdetrees/vignettes/kdetrees.pdf).  On average, there are approximately 20'000 remaining trees in each dataset.  The trees are stored in [Newick format](https://en.wikipedia.org/wiki/Newick_format).  

## Description of R Scripts
The `Software` directory contains the following set of R scripts:
* `func_ssh.R` is a set of supplemental functions called in running tropical PCA, and computing the proportion of variance explained (R^2) for each PCA method

## Running Tree PCA and Reproducing Projective Visualizations
Please ensure that the files in the `Data` and `Software` directories are downloaded, and set your working directory to this location.

The tree PCA routines essentially produce the same output for their respective metrics: a second-order principal component as a surface (more precisely: a locus, in the BHV setting; and a tropical triangle, in the tropical setting), whose boundaries are given by three points (trees); and projections of the dataset onto this component.  The software we provide in this repository creates the plots for the analyses of each of our 21 influenza datasets, given in the `Figures` directory.

BHV PCA is implemented by GeoPhytter+ on our datasets in a Linux terminal as follows:
```
for file in N_NYh3n2_HA_20000_5_*.txt
do 
java -classpath ".:*" geophytterplus.FitLFMTriangle $file -> ${file/.txt/.col}
sed "20020,20026d" ${file/.txt/.col} | sed "1,17d"  > ${file/.txt/.colc}
sed "12,20026d" ${file/.txt/.col} | sed "1,8d"  > ${file/.txt/.tree}
done
```
For each dataset, the first line creates a `.col` file which is made up of a `.colc` component, which gives the projected trees, and a `.tree` component, which gives the three points making up the boundaries of the BHV loci (second principal components).  The second and third lines separate the first `.col` file into these separate components.  Then, running
```
java -classpath ".:*" geophytterplus.DecomposeLFMTriangle $file 152 10 -> ${file/.tree/.tri}
```
gives a `.tri` file, which is used to create the BHV loci for the datasets.  The `.col`, `.tree`, and `.tri` files outputted by running GeoPhytter+ on our data is given in the `GeoPhytter+ Output` directory.

The `BHV_PCA_plot.R` script runs with data provided in the `GeoPhytter+ Output` directory, and creates figures for each of the datasets depicting the BHV loci and projected points, and figures for the tree topologies for the projected points as well as the vertices of the BHV locus.

The `Tropical_PCA_plot.R` script runs tropical PCA and creates figures for each of the datasets depicting the tropical triangles (as second principal components) and projected points, and figures for the tree topologies for the projected points as well as the vertices of the BHV locus.  For convenience, the output from tropPCA on our data required to plot the figures in the main text is incorporated directly into the script.

The `Figures` directory contains (edited) figures of the output of the BHV and tropical PCA routines applied to the 21 datasets provided in the `Data` directory.  They display the second-order principal components; the projections of the dataset onto these components; as well as the tree topologies of the three points defining the components and the projected points themselves.  Some of these figures are referenced in the main text.

## Calculating Proportion of Variance Explained
The `Trees_R_squared.R` script calculates and compares the percentage of variance explained (R^2) between the two methods.  These values are reported in Table 1 in the main text.

