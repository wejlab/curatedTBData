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

View summarized table for all the available data
``` r
data("DataSummary", package = "curatedTBData")
```

Load specific datasets:

``` r
library(curatedTBData)

# Load Single Study
get_curatedTBData(geo_access = "GSE39939")

# Load Multiple Studies
get_curatedTBData(geo_access = c("GSE39939","GSE107993"))
```

Load all available datasets in the pakcage:

``` r
library(curatedTBData)
get_curatedTBData(geo_access = "All")
```
