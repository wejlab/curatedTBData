if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
##### Read in raw data #####
geo <- "GSE42831"
sequencePlatform <- "GPL10558"
GSE42831_Non_normalized_counts <- GSE42831_Non_pvalue <- readRawData(geo, sequencePlatform)

gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
colnames(GSE42831_Non_pvalue) <- colnames(GSE42831_Non_normalized_counts) <-
  names(GEOquery::GSMList(gse))
