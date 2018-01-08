# trajr

[![Travis-CI Build Status](https://travis-ci.org/JimMcL/trajr.svg?branch=master)](https://travis-ci.org/JimMcL/trajr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/trajr)](https://cran.r-project.org/package=trajr)
[![codecov](https://codecov.io/gh/JimMcL/trajr/branch/master/graph/badge.svg)](https://codecov.io/gh/JimMcL/trajr)

Trajectory Analysis in R

An R package for analysing 2-dimensional animal trajectories, which represent the movements of animals through space and time. It provides functions to easily calculate various quantitative measures of trajectories such as speed and straightness or tortuosity. It also provides a range of other trajectory operations such as resampling to fixed step lengths (_rediscretization_), determining turning angles and step lengths, and so on.

`Trajr` does not provide support for creating trajectories; it operates on existing trajectories which are sequences of (x, y, time) coordinates.

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
```
Session info ----------------------------------------------------------------------------
 setting  value                       
 version  R version 3.4.2 (2017-09-28)
 system   x86_64, mingw32             
 ui       RStudio (1.1.383)           
 language (EN)                        
 collate  English_United States.1252  
 tz       Australia/Sydney            
 date     2017-12-24                  

Packages --------------------------------------------------------------------------------
 package   * version date       source        
 base      * 3.4.2   2017-09-28 local         
 compiler    3.4.2   2017-09-28 local         
 datasets  * 3.4.2   2017-09-28 local         
 devtools    1.13.4  2017-11-09 CRAN (R 3.4.3)
 digest      0.6.12  2017-01-27 CRAN (R 3.4.2)
 graphics  * 3.4.2   2017-09-28 local         
 grDevices * 3.4.2   2017-09-28 local         
 memoise     1.1.0   2017-04-21 CRAN (R 3.4.2)
 methods   * 3.4.2   2017-09-28 local         
 stats     * 3.4.2   2017-09-28 local         
 tools       3.4.2   2017-09-28 local         
 trajr     * 1.0.0   2017-12-21 local         
 utils     * 3.4.2   2017-09-28 local         
 withr       2.1.0   2017-11-01 CRAN (R 3.4.2)
 yaml        2.1.14  2016-11-12 CRAN (R 3.4.2)
```
