# trajr

[![Travis-CI Build Status](https://travis-ci.org/JimMcL/trajr.svg?branch=master)](https://travis-ci.org/JimMcL/trajr)

Trajectory Analysis in R

An R package for analysing 2-dimensional animal trajectories. It provides functions 
to easily calculate various measures of speed and straightness or tortuosity.

`Trajr` does not provide a mechanism to extract trajectories from videos, 
it requires trajectories to already exist as a sequence of x,y coordinates.


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
