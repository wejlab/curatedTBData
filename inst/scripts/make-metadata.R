# Make meta-data
if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}

base::read.dcf("DESCRIPTION", "Suggests") %>%
  base::gsub("\n", "", x = .) %>%
  base::strsplit(",") %>%
  base::unlist() %>%
  for (i in .) {
    if (!require(i, character.only = TRUE)) {
      BiocManager::install(i)
      require(i, character.only = TRUE)
    }
  }

####################################################

Title1 <- list.files("data/", pattern="GSE")
Title <- gsub("\\..*", "",Title1)

####################################################

Description <- rep(0,length(Title))
RDataClass <- rep(0,length(Title))
for (i in seq_len(length(Title))){
  tt_split <- strsplit(Title[i], "_")
  GSE <- tt_split[[1]][1]
  tt <- Title[i]

  if(length(grep("assay_raw", tt)) == 1){
    Description[i] <- paste("Raw transcriptome data for", GSE)
    RDataClass[i] <- "matrix"
  }
  if(length(grep("assay_reprocess", tt)) == 1){
    Description[i] <- paste("Reprocessed transcriptome data for", GSE)
    RDataClass[i] <- "matrix"
  }
  if(length(grep("column_data", tt)) == 1){
    Description[i] <- paste("Clinical information for samples from ", GSE)
    RDataClass[i] <- "DFrame"
  }
  if(length(grep("meta_data", tt) == 1)){
    Description[i] <- paste("Meta data information (MIAME object) for", GSE)
    RDataClass[i] <- "MIAME"
  }
  if(length(grep("row_data", tt) == 1)){
    Description[i] <- paste("Probe set information for samples from", GSE)
    RDataClass[i] <- "DFrame"
  }
  if(length(grep("SCAN_counts", tt) == 1)){
    Description[i] <- paste("SCAN-normalized transcriptome data for", GSE)
    RDataClass[i] <- "matrix"
  }
}

####################################################

BiocVersion <-
  BiocManager::version() %>%
  base::as.character()

####################################################

Genome <- base::as.character("9606")

####################################################

SourceType <- base::as.character("TXT")
series_accession <- lapply(strsplit(Title,"_"), function(x) x[1]) %>% unlist() %>% unique()
series_accession_reduce <- series_accession[-grep("GSEXXXXX",series_accession)]
accession_nchar <-
  stringr::str_length(series_accession_reduce)

accession_short <-
  stringr::str_trunc(series_accession_reduce, 5, ellipsis = "") %>%
  stringr::str_pad(accession_nchar, side = "right", pad = "n")

SourceUrl_temp <- stringr::str_c("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                            accession_short, "/", series_accession_reduce, "/matrix/",
                            series_accession_reduce, "_series_matrix.txt.gz")
SourceUrl_sub <- rep(SourceUrl_temp, times =
                   lapply(strsplit(Title[-grep("GSEXXXX", Title)],"_"),
                          function(x) x[1]) %>% unlist() %>%
                   table() %>% as.vector())
SourceUrl <- c(SourceUrl_sub, rep(NA, times = length(grep("GSEXXXX", Title))))


####################################################
SourceVersion  <- base::as.character(NA)

####################################################

Species <- base::as.character("Homo Sapiens")

####################################################

TaxonomyId <- base::as.character("9606")

####################################################

Coordinate_1_based <- base::as.logical(NA)

####################################################

DataProvider <- base::as.character("The National Center for Biotechnology Information")

####################################################

Maintainer <- base::as.character("Xutao Wang <xutaow@bu.edu>")

####################################################

RDataClass <- RDataClass

####################################################

DispatchClass <- base::as.character("Rda")

####################################################

RDataPath <- paste0("curatedTBData/", Title1)

####################################################
metadata <- data.frame(Title, Description, BiocVersion, Genome, SourceType,
                   SourceUrl, SourceVersion, Species, TaxonomyId,
                   Coordinate_1_based, DataProvider, Maintainer, RDataClass,
                   DispatchClass, RDataPath)
utils::write.csv(metadata, "inst/extdata/metadata.csv", row.names = FALSE)

