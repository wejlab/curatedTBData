# curatedTBData

The curatedTBData package collects 49 transcriptomic studies, with ongoing additions, focusing on patients with tuberculosis (TB) and other clinical conditions.

## Installation

curatedTBData is now available on Bioconductor. You can install the package via:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedTBData")
```


## Load Data

View summarized table for all the available data:
``` r
data("DataSummary", package = "curatedTBData")
View(DataSummary)
```

Load studies:

``` r
library(curatedTBData)

# List of all available resources with dryrun = TRUE
curatedTBData("", dryrun = TRUE)

# Load full version of single study
curatedTBData("GSE39939", dryrun = FALSE, curated.only = FALSE)

# Load curated version of multiple studies
curatedTBData(c("GSE39939","GSE107993"), dryrun = FALSE, curated.only = TRUE)
```
