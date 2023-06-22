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
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.2 (2022-10-31 ucrt)
 os       Windows 10 x64 (build 22621)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_Australia.utf8
 ctype    English_Australia.utf8
 tz       Australia/Sydney
 date     2023-06-22
 rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
 pandoc   3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version    date (UTC) lib source
 cachem        1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
 callr         3.7.3      2022-11-02 [1] CRAN (R 4.2.2)
 cli           3.6.0      2023-01-09 [1] CRAN (R 4.2.2)
 crayon        1.5.2      2022-09-29 [1] CRAN (R 4.2.2)
 devtools      2.4.5      2022-10-11 [1] CRAN (R 4.2.2)
 digest        0.6.31     2022-12-11 [1] CRAN (R 4.2.2)
 ellipsis      0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
 evaluate      0.20       2023-01-17 [1] CRAN (R 4.2.2)
 fastmap       1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
 fs            1.6.1      2023-02-06 [1] CRAN (R 4.2.3)
 glue          1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
 htmltools     0.5.4      2022-12-07 [1] CRAN (R 4.2.2)
 htmlwidgets   1.6.1      2023-01-07 [1] CRAN (R 4.2.2)
 httpuv        1.6.8      2023-01-12 [1] CRAN (R 4.2.2)
 knitr         1.42       2023-01-25 [1] CRAN (R 4.2.2)
 later         1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
 lifecycle     1.0.3      2022-10-07 [1] CRAN (R 4.2.2)
 magrittr      2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
 MASS          7.3-58.1   2022-08-03 [2] CRAN (R 4.2.2)
 memoise       2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
 mime          0.12       2021-09-28 [1] CRAN (R 4.2.0)
 miniUI        0.1.1.1    2018-05-18 [1] CRAN (R 4.2.2)
 pkgbuild      1.4.0      2022-11-27 [1] CRAN (R 4.2.2)
 pkgload       1.3.2      2022-11-16 [1] CRAN (R 4.2.2)
 prettyunits   1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
 processx      3.8.0      2022-10-26 [1] CRAN (R 4.2.2)
 profvis       0.3.7      2020-11-02 [1] CRAN (R 4.2.1)
 promises      1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
 ps            1.7.2      2022-10-26 [1] CRAN (R 4.2.2)
 purrr         1.0.1      2023-01-10 [1] CRAN (R 4.2.2)
 R6            2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
 Rcpp          1.0.10     2023-01-22 [1] CRAN (R 4.2.2)
 remotes       2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
 rlang         1.1.0      2023-03-14 [1] CRAN (R 4.2.3)
 rmarkdown     2.20       2023-01-19 [1] CRAN (R 4.2.2)
 rstudioapi    0.14       2022-08-22 [1] CRAN (R 4.2.2)
 sessioninfo   1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
 shiny         1.7.4      2022-12-15 [1] CRAN (R 4.2.2)
 signal        0.7-7      2021-05-25 [1] CRAN (R 4.2.0)
 stringi       1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
 stringr       1.5.0      2022-12-02 [1] CRAN (R 4.2.2)
 urlchecker    1.0.1      2021-11-30 [1] CRAN (R 4.2.2)
 usethis       2.1.6      2022-05-25 [1] CRAN (R 4.2.2)
 vctrs         0.5.2      2023-01-23 [1] CRAN (R 4.2.2)
 xfun          0.36       2022-12-21 [1] CRAN (R 4.2.2)
 xtable        1.8-4      2019-04-21 [1] CRAN (R 4.2.2)
 yaml          2.3.7      2023-01-23 [1] CRAN (R 4.2.2)
 ```
