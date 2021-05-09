if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
geo <- "GSE54992"
sequencePlatform <- "GPL570"
GSE54992_normalized_rma <- readRawData(geo, sequencePlatform)
colnames(GSE54992_normalized_rma) <- gsub("_.*", "", colnames(GSE54992_normalized_rma))
GSE54992_normalized_data <- GSE54992_normalized_rma

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "Tissue")
TBStatus <- TBStatus_temp <- characteristic_data_frame$TBStatus
MeasurementTime <- rep("Baseline", nrow(characteristic_data_frame))
Treatment <- rep("Untreated", nrow(characteristic_data_frame))
for (i in seq_len(length(TBStatus_temp))) {
  if (TBStatus_temp[i] == "tuberculosis") {
    TBStatus[i] <- "PTB"
    Treatment[i] <- "Pre-treatment"
  } else if (TBStatus_temp[i] == "TB patient after anti-tuberculosis treatment for 3 months") {
    TBStatus[i] <- "PTB"
    Treatment[i] <- "anti-tuberculosis treatment"
    MeasurementTime[i] <- "3 Month(s)"
  } else if (TBStatus_temp[i] == "latent tuberculosis infection") {
    TBStatus[i] <- "LTBI"
  } else if (TBStatus_temp[i] == "healthy donor") {
    TBStatus[i] <- "Control"
  } else if (TBStatus_temp[i] == "TB patient after anti-tuberculosis treatment for 6 months"){
    TBStatus[i] <- "PTB"
    Treatment[i] <- "anti-tuberculosis treatment"
    MeasurementTime[i] <- "6 Month(s)"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Treatment <- Treatment
characteristic_data_frame$MeasurementTime <- MeasurementTime
data_title <- sapply(seq_len(length(GEOquery::GSMList(gse))), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title)
characteristic_data_frame$PatientID <- sapply(strsplit(data_title, " "),
                                              function(x) paste0(x[1:2], collapse = ""))

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE54992_normalized_rma, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create Metadata #####
GSE54992_experimentData <- new("MIAME",
                               name = "Xinchun Chen",
                               lab = "Shenzhen Third People's Hospital",
                               contact = "chenxinchun@hotmail.com",
                               title = "Increased complement C1q level marks active disease in human tuberculosis.",
                               abstract = "Complement gene expression in peripheral blood mononuclear cells of tuberculosis patients and controls were determined using whole genome transcriptional microarray assays. The mRNA and protein levels of three C1q components, C1qA, C1qB, and C1qC, were further validated by qRT-PCR and enzyme-linked immunosorbent assay, respectively. The percentages of C1q expression in CD14 positive cells were determined by flow cytometry. Finally, C1qC protein level was quantified in the pleural fluid of tuberculosis and non-tuberculosis pleurisy.",
                               url = "10.1371/journal.pone.0092340",
                               pubMedIds = "24647646",
                               other=list(Platform = "Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)"))

GSE54992_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE54992_normalized_data= as.matrix(GSE54992_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE54992_experimentData));GSE54992_sobject
save_raw_files(GSE54992_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
