if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
##### Read in raw data #####
geo <- "GSE42825"
sequencePlatform <- "GPL10558"
GSE42825_Non_normalized_counts <- GSE42825_Non_pvalue <- readRawData(geo, sequencePlatform)

gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
colnames(GSE42825_Non_pvalue) <- colnames(GSE42825_Non_normalized_counts) <-
  names(GEOquery::GSMList(gse))

##### Create Column Data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c('Gender','Ethnicity','TBStatus','Tissue')

TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
for (i in 1:length(TBStatus)) {
  if (TBStatus_temp[i] == "TB") {
    TBStatus[i] <- "PTB"
  }
  if (TBStatus_temp[i] == "Active sarcoidosis") {
    TBStatus[i] <- "OD"
  }
  if (TBStatus_temp[i] == "Non-active sarcoidosis") {
    TBStatus[i] <- "OD"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Tissue <- "Whole Blood"
od <- rep("Negative", nrow(characteristic_data_frame))
od_status <- TBStatus_temp
for(i in 1:length(od_status)) {
  if (od_status[i] == "Active sarcoidosis") {
    od[i] <- "Positive (Active)"
  } else if (od_status[i] == "Non-active sarcoidosis") {
    od[i] <- "Positive (Non-active)"
  }
}
characteristic_data_frame$SarcoidosisStatus <- od
gender <- characteristic_data_frame$Gender
characteristic_data_frame$Gender <- ifelse(gender == "M", "Male", "Female")

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE42825_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create Meta Data #####
GSE42825_experimentData <- methods::new("MIAME",
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

GSE42825_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE42825_Non_normalized_counts = as.matrix(GSE42825_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE42825_experimentData));GSE42825_sobject
save_raw_files(GSE42825_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)







