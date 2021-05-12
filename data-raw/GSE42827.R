if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE42827"
sequencePlatform <- "GPL10558"
GSE42827_data_list <- readRawData(geo, sequencePlatform)
GSE42827_Non_normalized_data <- GSE42827_Non_pvalue <- GSE42827_data_list$data_Non_normalized

gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
colnames(GSE42827_Non_pvalue) <- colnames(GSE42827_Non_normalized_data) <-
  names(GEOquery::GSMList(gse))

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Treatment", "PatientID", "Tissue", "TBStatus")
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$PneumoniaStatus <- "Positive"
characteristic_data_frame$TBStatus <- "OD"

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create row data #####
row_data <- map_gene_symbol(GSE42827_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create meta data #####
GSE42827_experimentData <- methods::new("MIAME",
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

GSE42827_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE42827_Non_normalized_data = as.matrix(GSE42827_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE42827_experimentData));GSE42827_sobject
save_raw_files(GSE42827_sobject, path = "data-raw/", geo = geo)
##### Create normalized curated assay #####
GSE42827_normed <- GSE42827_data_list$data_normalized
colnames(GSE42827_normed) <- names(GEOquery::GSMList(gse))
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE42827_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
