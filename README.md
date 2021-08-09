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

Load studies:

``` r
library(curatedTBData)

# Load all studies
curatedTBData(DataSummary$`GEO accession`)

# Load single study
curatedTBData("GSE39939")

# Load multiple studies
curatedTBData(c("GSE39939","GSE107993"))
```
