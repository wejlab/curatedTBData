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
# Large data file 433.2MB
utils::download.file(as.character(urls$url[1]), temp)
utils::untar(temp, exdir = tempd)
indexStart <- grep("GSM851850", list.files(tempd))
filesRNA_reduce <- list.files(tempd)[indexStart:length(list.files(tempd))]
filesRNA <- list.files(tempd, full.names = TRUE)[indexStart:length(list.files(tempd))]
# Reading Single-Channel Agilent Intensity Data
GSE34608_RNA <- limma::read.maimages(filesRNA, source="agilent", green.only = TRUE)
GSE34608_Non_normalized_data <- GSE34608_RNA$E
row.names(GSE34608_Non_normalized_data) <- GSE34608_RNA$genes$ProbeName
colnames(GSE34608_Non_normalized_data) <- gsub("_.*", "", filesRNA_reduce)
GSE34608_Non_pvalue <- GSE34608_Non_normalized_data

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame_full <- readRawColData(gse)
index1 <- grep("GSM851850", row.names(characteristic_data_frame_full))
characteristic_data_frame <- characteristic_data_frame_full[index1: nrow(characteristic_data_frame_full),]
colnames(characteristic_data_frame) <- c("Gender", "Tissue", "TBStatus", "PatientID")
characteristic_data_frame$Tissue <- "PBMC"
TBStatus <- TBStatus_temp <- characteristic_data_frame$TBStatus
for (i in 1:length(TBStatus_temp)) {
  if (TBStatus_temp[i] == "control") {
    TBStatus[i] <- "Control"
  } else if (TBStatus_temp[i] == "sarcoidosis") {
    TBStatus[i] <- "OD"
  } else {
    TBStatus[i] <- "PTB"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
SarcoidosisStatus <- ifelse(TBStatus == "OD", "Positive", "Negative")
characteristic_data_frame$SarcoidosisStatus <- SarcoidosisStatus
characteristic_data_frame$GeographicalRegion <- "Germany"
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
gpl6480 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
gpl6480_annotate <- gpl6480@dataTable@table

GSE34608_probeNames <- data.frame(row.names(GSE34608_Non_pvalue))
colnames(GSE34608_probeNames) <- c("ID_REF")

new_row_data <- GSE34608_probeNames %>%
  dplyr::left_join(gpl6480_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$GENE_SYMBOL

##### Create Metadata #####
GSE34608_experimentData <- methods::new("MIAME",
                                        name = "Jeroen Maertzdorf",
                                        lab = "MPIIB",
                                        contact = "maertzdorf@mpiib-berlin.mpg.de",
                                        title = "Common patterns and disease-related signatures in tuberculosis and sarcoidosis.",
                                        abstract = "Gene and microRNA expression analysis of whole blood RNA from tuberculosis and sarcoidosis patients and healthy controls",
                                        url = "10.1073/pnas.1121072109",
                                        pubMedIds = "22547807",
                                        other = list(Platform = "Agilent-019118 Human miRNA Microarray 2.0 G4470B (Feature Number version) (GPL6480)"))
GSE34608_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE34608_Non_normalized_data = as.matrix(GSE34608_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE34608_experimentData));GSE34608_sobject
save_raw_files(GSE34608_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

curatedExprs <- makeCuratedExprs(row_data = new_row_data,
                                 data_Non_normalized = GSE34608_Non_normalized_data,
                                 dataType = "Microarray", platform = "Agilent",
                                 method = "quantile", FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
