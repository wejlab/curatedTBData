if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE62525"
sequencePlatform <- "GPL16951"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_cel <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_cel,temp)
utils::untar(temp,exdir = tempd)

gprFiles <- list.files(path = tempd, pattern = '*.gpr', full.names = TRUE)
# Label cy5
GSE62525_RG <- limma::read.maimages(gprFiles, source = "genepix",
                                    columns = list(R = "F635 Median", G = "B635 Median"))
GSE62525_RG_background <- limma::backgroundCorrect(GSE62525_RG, method = "normexp")
GSE62525_RG_background
