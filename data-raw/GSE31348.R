if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

geo <- "GSE31348"
sequencePlatform <- "GPL570"
GSE31348_normalized_rma <- readRawData(geo, sequencePlatform)
colnames(GSE31348_normalized_rma) <- gsub("_.*", "", colnames(GSE31348_normalized_rma))
GSE31348_normalized_data <- GSE31348_normalized_rma

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "PatientID",
                                         "MeasurementTime", "Tissue")
characteristic_data_frame$TBStatus <- "PTB"
characteristic_data_frame$Tissue <- "Whole Blood"

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE31348_normalized_rma, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

GSE31348_experimentData <- new("MIAME",
                               name = "Jackie Cliff",
                               lab = "London School of Hygiene & Tropical Medicine",
                               contact = "jackie.cliff@lshtm.ac.uk",
                               title = "Distinct phases of blood gene expression pattern through tuberculosis treatment reflect modulation of the humoral immune response.",
                               abstract = "Ex vivo blood samples analysed during TB treatment. These samples are from 9 successfully cured patients at diagnosis and end-of-treatment at 26 weeks.",
                               url = "10.1093/infdis/jis499",
                               pubMedIds = "22872737",
                               other=list(Platform = "Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)"))
GSE31348_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE31348_normalized_data = as.matrix(GSE31348_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE31348_experimentData));GSE31348_sobject
save_raw_files(GSE31348_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE31348_normalized_data,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))





