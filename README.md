CepoManuscript <img src="Cepo_logo.png" align="right" width="225" height="250"/>
================

Defining the identity of a cell is fundamental to understand the heterogeneity of cells to various environmental signals and perturbations. We present Cepo, a new method to explore cell identities from single-cell RNA-sequencing data using differential stability as a new metric to define cell identity genes. Cepo computes cell-type specific gene statistics pertaining to differential stable gene expression. This repository contains all code required to reproduce the results in [our preprint](link). 

## Installation

You can install the development version of _Cepo_ that can be installed from GitHub
using the `remotes` package:

``` r
# install.packages("remotes")
remotes::install_github("PYangLab/Cepo")
```

To also build the vignettes use:

``` r
# install.packages("remotes")
remotes::install_github("PYangLab/Cepo", dependencies = TRUE,
                         build_vignettes = TRUE)
```

**NOTE:** Building the vignettes requires the installation of additional packages.

## Documentation

The documentation for _Cepo_ is available from http://github.io/PYangLab/Cepo

To view the vignette and all the package documentation for the development version visit http://github.io/HaniJieunKim/Cepo.

## Citing _Cepo_

If you use _Cepo_ in your work please cite
our preprint ["Kim H.J., Wang K., Yang P. Cepo uncovers cell identity through differential stability. bioRxiv DOI:][paper].

## Developers

The following individuals were involved in developing the Cepo package:

* [@HaniJieunKim](https://github.com/HaniJieunKim)
* [@kevinwang09](https://github.com/kevinwang09)
* [@PYangLab](https://github.com/PYangLab) 

## Contact us

If you have any enquiries, especially about using Cepo to analyse your data, please contact hani.kim@sydney.edu.au. We actively welcome any feedback and suggestions! 
