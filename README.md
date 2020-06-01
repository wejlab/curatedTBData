# curatedTBData

The curatedTBData package collects 42 transcriptome studies for patients with turberculosis (TB) and other clinical information.

## Installation

curatedTBData is under development. You can install the devel version via
GitHub:

``` r
install.packages("devtools")
devtools::install_github("compbiomed/curatedTBData")
```

## Load Data
To get a listing of all datasets:

``` r
library(curatedTBData)
data(package="curatedTBData")
```
Load specific datasets
``` r
data(GSE39939_sobject)
data(GSE107991_mobject)
```
