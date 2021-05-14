if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE25534"
sequencePlatform <- "GPL1708"
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_raw <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
# Large file 644.7 MB
utils::download.file(url_raw, temp)
utils::untar(temp, exdir = tempd)
files <- list.files(tempd, pattern = "GSM.*", full.names = TRUE)
GSE25534_raw <- limma::read.maimages(files, source = "agilent")
MA <- normalizeWithinArrays(GSE25534_raw)


