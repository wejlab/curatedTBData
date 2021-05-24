if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE83892"
sequencePlatform <- "GPL10558"
# Used the normalized data, cannot identify description ID to sample ID in the non-normalized data
GSE83892_data_list <- readRawData(geo, sequencePlatform)
GSE83892_Non_pvalues <- GSE83892_data_list$data_Non_normalized
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_ID <- as.character(urls$url[1])
temp_ID <- tempfile()
download.file(url_ID, temp_ID)
GSE83892_ID <- read.delim(gzfile(temp_ID), header = FALSE)
colnames(GSE83892_ID) <- c("SampleID", "DescriptionID")
GSE83892_ID$DescriptionID <- paste0("X", GSE83892_ID$DescriptionID)

indx <- match(GSE83892_ID$DescriptionID, colnames(GSE83892_Non_pvalues))
GSE83892_Non_pvalues <- GSE83892_Non_pvalues[,indx]
colnames(GSE83892_Non_pvalues) <- GSE83892_ID$SampleID
GSE83892_Non_normalized_data <- GSE83892_Non_pvalues

##### Create column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- "TBStatus"
IRIS_Status <- characteristic_data_frame$TBStatus
IRIS_Status <- ifelse(IRIS_Status == "TBM-IRIS", "Positive", "Negative")
characteristic_data_frame$IRIS_Status <- IRIS_Status
characteristic_data_frame$TBStatus <- ifelse(characteristic_data_frame$TBStatus == "Control",
                                             "Control", "PTB")
MeningitisStatus <- ifelse(characteristic_data_frame$TBStatus == "Control",
                           "Negative", "Positive")
characteristic_data_frame$MeningitisStatus <- MeningitisStatus
characteristic_data_frame$HIVStatus <- "Positive"
characteristic_data_frame$GeographicalRegion <- "South Africa"
characteristic_data_frame$Tissue <- "Whole Blood"
data_title <- sapply(seq_len(length(GEOquery::GSMList(gse))), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title)
MeasurementTime <- rep("Control Time", nrow(characteristic_data_frame))
MeasurementTime[grep("2wpART", data_title)] <- "2 Weeks Post antiretroviral therapy (Second time point)"
MeasurementTime[grep("2wpTBRx", data_title)] <- "2 Weeks Post TB treatment (First time point)"
MeasurementTime[grep("Diagnosis", data_title)] <- "TB meningitis diagnosis (Baseline)"
PatientID_raw <- sapply(strsplit(data_title, "_"), "[[", 2)
characteristic_data_frame$PatientID <- substr(PatientID_raw, start = 1,
                                              stop = nchar(PatientID_raw) - 1)
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE83892_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

GSE83892_experimentData <- methods::new("MIAME",
                                        name = "Rachel Lai",
                                        lab = "The Francis Crick Institute",
                                        contact = "rachel.lai@crick.ac.uk, rachel.lai@imperial.ac.uk",
                                        title = "Blood transcriptional signature of HIV-associated tuberculous meningitis.",
                                        abstract = "Patients with HIV-associated TB meningitis (TBM) are known to experience systemic hyperinflammation, clinically known as immune reconstitution inflammatory syndrome (IRIS), following the commencement of antiretroviral therapy (ART).",
                                        url = "10.1093/infdis/jiw561",
                                        pubMedIds = "27932622",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE83892_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE83892_Non_normalized_data = as.matrix(GSE83892_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE83892_experimentData));GSE83892_sobject
save_raw_files(GSE83892_sobject, path = "data-raw/", geo = geo)

##### Create normalized curated assay #####
GSE83892_normed <- GSE83892_data_list$data_normalized[,indx]
colnames(GSE83892_normed) <- GSE83892_ID$SampleID
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE83892_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)



