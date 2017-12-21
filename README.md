# trajr

[![Travis-CI Build Status](https://travis-ci.org/JimMcL/trajr.svg?branch=master)](https://travis-ci.org/JimMcL/trajr)

Trajectory Analysis in R

An R package for analysing 2-dimensional animal trajectories. It provides functions 
to easily calculate various measures of speed and straightness or tortuosity.

`Trajr` does not provide a mechanism to extract trajectories from videos, 
it requires trajectories to already exist as a sequence of x,y coordinates.


## Installation

    $ install.packages("devtools")
    $ devtools::install_github("JimMcL/trajr")
    
or else 

    $ devtools::install_github("JimMcL/trajr", build_vignettes = TRUE)
    
to install trajr from github with vignettes (which are user guides). Unfortunately, installing vignettes can be time-consuming.

## Using trajr

The best way to start is to read the package vignette. 

