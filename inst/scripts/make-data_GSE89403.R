if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

#### Read in raw data ####
geo <- "GSE89403"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[3]), temp)

GSE89403_Non_normalized <- read.csv(gzfile(temp), row.names = 1, header = TRUE)
temp2 <- tempfile()
download.file(as.character(urls$url[2]), temp2)
GSE89403_normalized_log2 <- read.csv(gzfile(temp2), row.names = 1, header = TRUE)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
data_list <- list(data_Non_normalized = GSE89403_Non_normalized,
                  data_normalized = GSE89403_normalized_log2)
# all(colnames(GSE89403_Non_normalized) == colnames(GSE89403_normalized_log2)) TRUE
# all(row.names(GSE89403_Non_normalized) == row.names(GSE89403_normalized_log2)) TRUE
# Match column names to sample name
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
description_id_raw_new <- sapply(description_id_raw, function(x) unlist(strsplit(x,"-"))[1])
names(description_id_raw_new) <- NULL
# Most of the samples have 2 technical replicates, 2 samples have 4 replicates
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id_raw_new)
ID_table_merge <- dplyr::distinct(ID_table, DescriptionID, .keep_all = TRUE)
data_list_matched <- lapply(data_list, function(x) {
  data_counts <- as.matrix(x[, -1])
  indx <- match(ID_table_merge$DescriptionID, colnames(data_counts))
  data_counts <- data_counts[, indx]
  colnames(data_counts) <- ID_table_merge$SampleID
  data_counts
})
GSE89403_Non_normalized_counts <- data_list_matched$data_Non_normalized
GSE89403_normalized_counts <- data_list_matched$data_normalized

#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue","SampleCode", "PatientID",
                                         "TBStatus","TreatmentResult","MeasurementTime",
                                         "TimeToNegativity", "mgit", "xpert", "tgrv")
characteristic_data_frame$mgit <- as.numeric(characteristic_data_frame$mgit)
characteristic_data_frame$xpert <- as.numeric(characteristic_data_frame$xpert)
characteristic_data_frame$tgrv <- as.numeric(characteristic_data_frame$tgrv)
characteristic_data_frame$Tissue <- "Whole Blood"
index <- grep("DX", characteristic_data_frame$MeasurementTime)
characteristic_data_frame$MeasurementTime[index] <- "Baseline"
TBStatus_temp <- TBStatus <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
LungDxStatus <- rep("Negative", nrow(characteristic_data_frame))
MTPStatus <- rep("Negative", nrow(characteristic_data_frame))
for(i in 1:length(TBStatus_temp)) {
  if(TBStatus_temp[i] == "Healthy Controls") {
    TBStatus[i] <- "Control"
  } else if(TBStatus_temp[i] == "Lung Dx Controls") {
    TBStatus[i] <- "OD"
    LungDxStatus[i] <- "Positive"
  } else if (TBStatus_temp[i] == "MTP Controls") {
    TBStatus[i] <- "OD"
    MTPStatus[i] <- "Positive"
  } else if(TBStatus_temp[i] == "TB Subjects") {
    TBStatus[i] <- "PTB"
  } else if(TBStatus_temp[i] == "NA") {
    TBStatus[i] <- NA
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$LungDxStatus <- LungDxStatus
characteristic_data_frame$MTPStatus <- MTPStatus

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(GSE89403_Non_normalized),
                                     SYMBOL_NEW = GSE89403_Non_normalized$symbol)
##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Daniel Zak",
                      lab = "Center for Infectious Disease Research",
                      contact = "danielzak000@gmail.com",
                      title = "Host blood RNA signatures predict the outcome of tuberculosis treatment.",
                      abstract = "Biomarkers for tuberculosis treatment outcome will assist in guiding individualized treatment and evaluation of new therapies. To identify candidate biomarkers, RNA sequencing of whole blood from a well-characterized TB treatment cohort was performed. Application of a validated transcriptional correlate of risk for TB revealed symmetry in host gene expression during progression from latent TB infection to active TB disease and resolution of disease during treatment, including return to control levels after drug therapy. The symmetry was also seen in a TB disease signature, constructed from the TB treatment cohort, that also functioned as a strong correlate of risk. Both signatures identified patients at risk of treatment failure 1-4 weeks after start of therapy. Further mining of the transcriptomes revealed an association between treatment failure and suppressed expression of mitochondrial genes before treatment initiation, leading to development of a novel baseline (pre-treatment) signature of treatment failure. These novel host responses to TB treatment were integrated into a five-gene real-time PCR-based signature that captures the clinically relevant responses to TB treatment and provides a convenient platform for stratifying patients according to their risk of treatment failure. Furthermore, this 5-gene signature is shown to correlate with the pulmonary inflammatory state (as measured by PET-CT) and can complement sputum-based Gene Xpert for patient stratification, providing a rapid and accurate alternative to current methods.",
                      url = "10.1016/j.tube.2017.08.004",
                      pubMedIds = "29050771",
                      other=list(Platform = "Illumina HiSeq 2000 (Homo sapiens) (GPL11154)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE89403_Non_normalized_counts)),
  #  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(new_col_info, paste0("data-raw/", geo, "_column_data.RDS"))

##### Create curated matched assay #####
GSE89403_normalized_counts_matched <- probesetsToGenes(new_row_data, GSE89403_normalized_counts,
                                                       median)
saveRDS(GSE89403_normalized_counts_matched, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))
