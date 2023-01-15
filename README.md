# trajr

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/trajr)](https://cran.r-project.org/package=trajr)
![CRAN download count](https://cranlogs.r-pkg.org/badges/trajr)
[![R-CMD-check](https://github.com/JimMcL/trajr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JimMcL/trajr/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/JimMcL/trajr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/JimMcL/trajr?branch=master)
<!-- badges: end -->

<!-- To display total CRAN downloads, use
https://cranlogs.r-pkg.org/badges/grand-total/trajr
-->

Trajectory Analysis in R

An R package for analysing 2-dimensional animal trajectories, which represent the movements of animals through space and time. It provides functions to easily calculate various quantitative measures of trajectories such as speed and straightness or tortuosity. It also provides a range of other trajectory operations such as resampling to fixed step lengths (_rediscretization_), determining turning angles and step lengths, and so on. If you use `trajr` in your research, please cite [McLean, D. J., & Skowron Volponi, M. A. (2018). trajr: An R package for characterisation of animal trajectories. Ethology, 124(6), 440-448. https://doi.org/10.1111/eth.12739](https://doi.org/10.1111/eth.12739). 

`trajr` does not provide functionality to create trajectories; it operates on existing trajectories that are sequences of `(x, y, time)` coordinates. It does, however, provide some functionality to generate random or constrained trajectories for simulations or for testing of analyses.

Refer to [NEWS.md](NEWS.md) for a list of changes in each version. The article [trajr: An R package for characterisation of animal trajectories](https://doi.org/10.1111/eth.12739) described `trajr` version 1.0.0. Please refer to [NEWS.md](NEWS.md) for a brief summary of what has changed since the article was written.

## Installation
    $ install.packages("trajr")

Or to install the latest development version directly from Github:

    $ install.packages("devtools")
    $ devtools::install_github("JimMcL/trajr")
    
Or else, to install from Github **and** get the user documentation (vignette):

    $ install.packages("devtools")
    $ devtools::install_github("JimMcL/trajr", build_vignettes = TRUE)
    
Unfortunately, installing vignettes can be time-consuming.

## Using trajr

The best way to start is to read [our paper](https://doi.org/10.1111/eth.12739) and the package vignette. After installation from CRAN (or from github including vignettes), the vignette can be accessed by running `browseVignettes(package = "trajr")`. The vignette for the CRAN version is available on CRAN at https://cran.r-project.org/web/packages/trajr/vignettes/trajr-vignette.html. 

A vignette demonstrating trajectory simulation is available online at https://cran.r-project.org/web/packages/trajr/vignettes/simulations.html. 

## Environment
<!-- Output from devtools::session_info() -->
```
- Session info --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_Australia.1252      
 ctype    English_Australia.1252      
 tz       Australia/Sydney            
 date     2020-12-16                  

- Packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 package     * version date       lib source        
 assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
 backports     1.1.7   2020-05-13 [1] CRAN (R 4.0.0)
 callr         3.4.3   2020-03-28 [1] CRAN (R 4.0.2)
 cli           2.0.2   2020-02-28 [1] CRAN (R 4.0.2)
 crayon        1.3.4   2017-09-16 [1] CRAN (R 4.0.2)
 desc          1.2.0   2018-05-01 [1] CRAN (R 4.0.2)
 devtools      2.3.0   2020-04-10 [1] CRAN (R 4.0.2)
 digest        0.6.25  2020-02-23 [1] CRAN (R 4.0.2)
 ellipsis      0.3.1   2020-05-15 [1] CRAN (R 4.0.2)
 fansi         0.4.1   2020-01-08 [1] CRAN (R 4.0.2)
 fs            1.4.2   2020-06-30 [1] CRAN (R 4.0.2)
 glue          1.4.1   2020-05-13 [1] CRAN (R 4.0.2)
 magrittr      1.5     2014-11-22 [1] CRAN (R 4.0.2)
 memoise       1.1.0   2017-04-21 [1] CRAN (R 4.0.2)
 packrat       0.5.0   2018-11-14 [1] CRAN (R 4.0.2)
 pkgbuild      1.0.8   2020-05-07 [1] CRAN (R 4.0.2)
 pkgload       1.1.0   2020-05-29 [1] CRAN (R 4.0.2)
 prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.0.2)
 processx      3.4.3   2020-07-05 [1] CRAN (R 4.0.2)
 ps            1.3.3   2020-05-08 [1] CRAN (R 4.0.2)
 R6            2.4.1   2019-11-12 [1] CRAN (R 4.0.2)
 remotes       2.1.1   2020-02-15 [1] CRAN (R 4.0.2)
 rlang         0.4.7   2020-07-09 [1] CRAN (R 4.0.2)
 rprojroot     1.3-2   2018-01-03 [1] CRAN (R 4.0.2)
 rstudioapi    0.11    2020-02-07 [1] CRAN (R 4.0.2)
 sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
 testthat      2.3.2   2020-03-02 [1] CRAN (R 4.0.2)
 tinytex       0.24    2020-06-20 [1] CRAN (R 4.0.2)
 usethis       1.6.1   2020-04-29 [1] CRAN (R 4.0.2)
 withr         2.2.0   2020-04-20 [1] CRAN (R 4.0.2)
 xfun          0.15    2020-06-21 [1] CRAN (R 4.0.2)
 ```
