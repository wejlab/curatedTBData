if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE42830"
sequencePlatform <- "GPL10558"
GSE42830_data_list <- readRawData(geo, sequencePlatform)
GSE42830_Non_normalized_data <- GSE42830_Non_pvalue <- GSE42830_data_list$data_Non_normalized

gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
colnames(GSE42830_Non_pvalue) <- colnames(GSE42830_Non_normalized_data) <-
  names(GEOquery::GSMList(gse))

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Gender", "Ethnicity", "TBStatus", "Tissue")

TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)

for(i in 1:length(TBStatus)){
  if(TBStatus_temp[i] == "TB"){
    TBStatus[i] <- "PTB"
  }
  if(TBStatus_temp[i] != "TB" & TBStatus_temp[i] != "Control"){
    TBStatus[i] <- "OD"
  }
}
characteristic_data_frame$TBStatus <- TBStatus

characteristic_data_frame$Tissue <- "Whole Blood"

SarcoidosisStatus <- rep("Negative", nrow(characteristic_data_frame))
SarcoidosisStatus[which(TBStatus_temp == "Non-active sarcoidosis")] <- "Positive (Non-active)"
SarcoidosisStatus[which(TBStatus_temp == "Active Sarcoid")] <- "Positive (Active)"
characteristic_data_frame$SarcoidosisStatus <- SarcoidosisStatus
LungCancerStatus <- rep("Negative", nrow(characteristic_data_frame))
LungCancerStatus[which(TBStatus_temp == "lung cancer")] <- "Positive"
characteristic_data_frame$LungCancerStatus <- LungCancerStatus
PneumoniaStatus <- rep("Negative", nrow(characteristic_data_frame))
PneumoniaStatus[which(TBStatus_temp == "Baseline")] <- "Positive"
characteristic_data_frame$PneumoniaStatus <- PneumoniaStatus

characteristic_data_frame$Gender <- ifelse(characteristic_data_frame$Gender == "M",
                                           "Male", "Female")
characteristic_data_frame$GeographicalRegion <- "Germany"

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create row data #####
row_data <- map_gene_symbol(GSE42830_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create meta data #####
GSE42830_experimentData <- methods::new("MIAME",
                                        name = "Chole Bloom",
                                        lab = "MRC National Institute for Medical Research",
                                        contact = "cbloom@nimr.mrc.ac.uk",
                                        title = "Transcriptional blood signatures distinguish pulmonary tuberculosis, pulmonary sarcoidosis, pneumonias and lung cancers.",
                                        abstract = "An Interferon-inducible neutrophil-driven blood transcriptional signature was present in both sarcoidosis and tuberculosis, with a higher abundance and expression in tuberculosis.
                                        Heterogeneity of the sarcoidosis signature correlated significantly with disease activity.
                                        Transcriptional profiles in pneumonia and lung cancer revealed an over-abundance of inflammatory transcripts. After successful treatment the transcriptional activity in tuberculosis and pneumonia patients was significantly reduced.
                                        However the glucocorticoid-responsive sarcoidosis patients showed a significant increase in transcriptional activity.
                                        144-blood transcripts were able to distinguish tuberculosis from other lung diseases and controls.",
                                        url = "10.1371/journal.pone.0070630",
                                        pubMedIds = "23940611",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))

GSE42830_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE42830_Non_normalized_data = as.matrix(GSE42830_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE42830_experimentData));GSE42830_sobject
save_raw_files(GSE42830_sobject, path = "data-raw/", geo = geo)
##### Create normalized curated assay #####
GSE42830_normed <- GSE42830_data_list$data_normalized
colnames(GSE42830_normed) <- names(GEOquery::GSMList(gse))
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE42830_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)



