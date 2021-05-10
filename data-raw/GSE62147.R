if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE62147"
sequencePlatform <- "GPL6480"
temp <- tempfile()
tempd <- tempdir()
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
# Large data file 466.8MB
utils::download.file(as.character(urls$url[2]), temp)
utils::untar(temp, exdir = tempd)
filesPath <- list.files(tempd, pattern = "GSM.*", full.names = TRUE)
GSE62147_raw <- limma::read.maimages(filesPath, source="agilent", green.only = TRUE)
GSE62147_Non_normalized_data <- GSE62147_raw$E
row.names(GSE62147_Non_normalized_data) <- GSE62147_raw$genes$ProbeName
col_name1 <- gsub(".*/", "", colnames(GSE62147_Non_normalized_data))
colnames(GSE62147_Non_normalized_data) <- gsub("_.*", "", col_name1)
GSE62147_Non_pvalue <- GSE62147_Non_normalized_data

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Gender", "Tissue")
characteristic_data_frame$Tissue <- "PBMC"
characteristic_data_frame$MeasurementTime <- rep(c("post", "recruit"), nrow(characteristic_data_frame)/2)
characteristic_data_frame$Treatment <- "standard anti-TB chemotherapy"
title <- sapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title)
characteristic_data_frame$PatientID <- gsub("_.*", "", title)
characteristic_data_frame$GeographicalRegion <- "The Gambia"
indexOD <- grep("Maf", characteristic_data_frame$PatientID)
sputum_culture <- rep("M.tuberculosis", nrow(characteristic_data_frame))
sputum_culture[indexOD] <- " M.africanum"
characteristic_data_frame$sputum_culture <- sputum_culture
TBStatus <- rep("PTB", nrow(characteristic_data_frame))
TBStatus[indexOD] <- "OD"
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
gpl6480 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
gpl6480_annotate <- gpl6480@dataTable@table

GSE62147_probeNames <- data.frame(row.names(GSE62147_Non_pvalue))
colnames(GSE62147_probeNames) <- c("ID_REF")

new_row_data <- GSE62147_probeNames %>%
  dplyr::left_join(gpl6480_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$GENE_SYMBOL

##### Create Metadata #####
GSE62147_experimentData <- methods::new("MIAME",
                                        name = "Jeroen Maertzdorf",
                                        lab = "MPIIB",
                                        contact = "maertzdorf@mpiib-berlin.mpg.de",
                                        title = "Differential transcriptomic and metabolic profiles of M. africanum- and M. tuberculosis-infected patients after, but not before, drug treatment.",
                                        abstract = "Gene expression analysis of whole blood RNA from tuberculosis patients prior to and following treatment.",
                                        url = "10.1038/gene.2010.51",
                                        pubMedIds = "26043170",
                                        other = list(Platform = "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version) (GPL6480)"))
GSE62147_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE62147_Non_normalized_data = as.matrix(GSE62147_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE62147_experimentData));GSE62147_sobject
save_raw_files(GSE62147_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)


