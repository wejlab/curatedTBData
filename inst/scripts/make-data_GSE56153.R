if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE56153"
sequencePlatform <- "GPL6883"
GSE56153_data_list <- readRawData(geo, sequencePlatform)
GSE56153_Non_pvalue <-
  GSE56153_Non_normalized_data <- GSE56153_data_list$data_Non_normalized

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
colnames(GSE56153_Non_pvalue) <-
  colnames(GSE56153_Non_normalized_data) <- names(GEOquery::GSMList(gse))
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue", "TBStatus", "TreatmentResult")
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$GeographicalRegion <- "Indonesia"
characteristic_data_frame$HIVStatus <- "Negative"
characteristic_data_frame$DiabetesStatus <- "Negative"
characteristic_data_frame$TBStatus <- ifelse(characteristic_data_frame$TreatmentResult == "Control",
                                             "Control", "PTB")
Treatment <- rep(NA, nrow(characteristic_data_frame))
Treatment[grep("PTB", characteristic_data_frame$TBStatus)] <- "2HRZE/4H3R3"
characteristic_data_frame$Treatment <- Treatment
MeasurementTime <- rep(NA, nrow(characteristic_data_frame))
MeasurementTime[grep("Active", characteristic_data_frame$TreatmentResult)] <- "0 weeks"
MeasurementTime[grep("Treatment", characteristic_data_frame$TreatmentResult)] <- "8 weeks"
MeasurementTime[grep("Recover", characteristic_data_frame$TreatmentResult)] <- "28 weeks"
characteristic_data_frame$MeasurementTime <- MeasurementTime
TreatmentResult <- TreatmentResult_temp <- characteristic_data_frame$TreatmentResult
for (i in 1:length(TreatmentResult)) {
  if (TreatmentResult_temp[i] == "Active") {
    TreatmentResult[i] = "Pre-treatment"
  } else if (TreatmentResult_temp[i] == "Treatment") {
    TreatmentResult[i] = "Good-response"
  } else {
    TreatmentResult[i] = "Definite Cure"
  }
}
characteristic_data_frame$TreatmentResult <- TreatmentResult
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE56153_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create Metadata #####
GSE56153_experimentData <- methods::new("MIAME",
                                        name = "Ahmad Nazri Mohamed Naim",
                                        lab = "Genome institute of Singapore",
                                        contact = "mohamedan@gis.a-star.edu.sg",
                                        title = "Genome-wide expression profiling identifies type 1 interferon response pathways in active tuberculosis.",
                                        abstract = "Whole-genome expression profiling on a cohort of TB patients. Tuberculosis patients above 15 years of age were recruited from an outpatient tuberculosis clinic in central Jakarta (Indonesia). Randomly selected control subjects with the same sex and age (+/-10%) were recruited from neighboring households, with first degree relatives excluded.",
                                        url = "10.1371/journal.pone.0045839",
                                        pubMedIds = "23029268",
                                        other = list(Platform = "Illumina HumanRef-8 v3.0 expression beadchip (GPL6883)"))
GSE56153_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE56153_Non_normalized_data = as.matrix(GSE56153_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE56153_experimentData));GSE56153_sobject
save_raw_files(GSE56153_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
GSE56153_normed <- GSE56153_data_list$data_normalized
colnames(GSE56153_normed) <- names(GEOquery::GSMList(gse))
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE56153_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))




