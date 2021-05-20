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
View(DataSummary$`GEO accession`)
```

Load specific studies:

``` r
library(curatedTBData)

# Load Single Study
curatedTBData("GSE39939")

# Load Multiple Studies
curatedTBData(c("GSE39939","GSE107993"))
```

Load all available datasets in the pakcage:

``` r
library(curatedTBData)
get_curatedTBData(geo_access = "All")
```
