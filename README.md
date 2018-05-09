# trajr

[![Travis-CI Build Status](https://travis-ci.org/JimMcL/trajr.svg?branch=master)](https://travis-ci.org/JimMcL/trajr)
[![codecov](https://codecov.io/gh/JimMcL/trajr/branch/master/graph/badge.svg)](https://codecov.io/gh/JimMcL/trajr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/trajr)](https://cran.r-project.org/package=trajr)
![CRAN download count](http://cranlogs.r-pkg.org/badges/grand-total/trajr)

Trajectory Analysis in R

An R package for analysing 2-dimensional animal trajectories, which represent the movements of animals through space and time. It provides functions to easily calculate various quantitative measures of trajectories such as speed and straightness or tortuosity. It also provides a range of other trajectory operations such as resampling to fixed step lengths (_rediscretization_), determining turning angles and step lengths, and so on. If you use `trajr` in your research, please cite [McLean, D. J., & Skowron Volponi, M. A. (2018). trajr: An R package for characterisation of animal trajectories. Ethology, 124(6), 440-448. https://doi.org/10.1111/eth.12739](https://doi.org/10.1111/eth.12739). 

`trajr` does not provide functionality to create trajectories; it operates on existing trajectories which are sequences of `(x, y, time)` coordinates.

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

The best way to start is to read the package vignette. After installation from CRAN (or from github including vignettes), run `browseVignettes(package = "trajr")`.

## Environment
<!-- Output from devtools::session_info() -->
```
Session info ----------------------------------------------------------------------------------------------------------------
 setting  value                       
 version  R version 3.5.0 (2018-04-23)
 system   x86_64, mingw32             
 ui       RStudio (1.1.447)           
 language (EN)                        
 collate  English_Australia.1252      
 tz       Australia/Sydney            
 date     2018-05-09                  

Packages --------------------------------------------------------------------------------------------------------------------
 package   * version    date       source        
 base      * 3.5.0      2018-04-23 local         
 compiler    3.5.0      2018-04-23 local         
 datasets  * 3.5.0      2018-04-23 local         
 devtools    1.13.5     2018-02-18 CRAN (R 3.5.0)
 digest      0.6.15     2018-01-28 CRAN (R 3.5.0)
 graphics  * 3.5.0      2018-04-23 local         
 grDevices * 3.5.0      2018-04-23 local         
 memoise     1.1.0      2017-04-21 CRAN (R 3.5.0)
 methods   * 3.5.0      2018-04-23 local         
 stats     * 3.5.0      2018-04-23 local         
 tools       3.5.0      2018-04-23 local         
 trajr     * 1.1.0      2018-05-07 local         
 utils     * 3.5.0      2018-04-23 local         
 withr       2.1.2      2018-03-15 CRAN (R 3.5.0)
 yaml        2.1.18     2018-03-08 CRAN (R 3.5.0)```
