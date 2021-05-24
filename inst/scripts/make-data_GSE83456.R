if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE83456"
sequencePlatform <- "GPL10558"
GSE83456_data_raw <- readRawData(geo, sequencePlatform, matrix.only = TRUE)
# Remove p-value manually
indexPvalue <- c(2, grep("\\.1", colnames(GSE83456_data_raw)))

# Create normalized data
GSE83456_Non_pvalues <- GSE83456_data_raw[, -indexPvalue]
colnames(GSE83456_Non_pvalues)[1] <- "X9370786041_A"
xr <- new("EListRaw", list(E = GSE83456_Non_pvalues,
                           other = list(Detection = GSE83456_data_raw[, indexPvalue])))
yr <- limma::neqc(xr)
##### Match colnames to sample ID #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(seq_len(length(names(GEOquery::GSMList(gse)))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$description[1])
description_id <- paste0("X", description_id_raw)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)
indx <- match(ID_table$DescriptionID, colnames(GSE83456_Non_pvalues))
GSE83456_Non_pvalues <- GSE83456_Non_pvalues[,indx]
colnames(GSE83456_Non_pvalues) <- ID_table$SampleID
GSE83456_Non_normalized_data <- GSE83456_Non_pvalues

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Gender", "Ethnicity", "Age", "TBStatus",
                                         "Tissue")
characteristic_data_frame$Gender <- ifelse(characteristic_data_frame$Gender == "M",
                                           "Male", "Female")
characteristic_data_frame$Tissue <- "Whole Blood"

TBStatus <- TBStatus_temp <- characteristic_data_frame$TBStatus
for(i in 1:length(characteristic_data_frame$TBStatus)) {
  if(TBStatus_temp[i] == "EPTB"){TBStatus[i] <-  "EPTB"}
  if(TBStatus_temp[i] == "HC"){TBStatus[i] <- "Control"}
  if(TBStatus_temp[i] == "PTB"){TBStatus[i] <- "PTB"}
  if(TBStatus_temp[i] == "Sarcoid"){TBStatus[i] <- "OD"}
}

characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
characteristic_data_frame$SarcoidosisStatus <- ifelse(characteristic_data_frame$TBStatus == "OD",
                                                      "Positive","Negative")
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE83456_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

GSE83456_experimentData <- methods::new("MIAME",
                                        name = "Krzysztof Potempa",
                                        lab = "National Institute for Medical Research (NIMR)",
                                        contact = "potempak@gmail.com",
                                        title = "The Transcriptional Signature of Active Tuberculosis Reflects Symptom Status in Extra-Pulmonary and Pulmonary Tuberculosis.",
                                        abstract = "Novel cohort comprising 61 healthy human controls, 47 human with EPTB, 45 human with PTB, 49 human with Sarcoid was recruited.",
                                        url = "10.1371/journal.pone.0162220",
                                        pubMedIds = "27706152",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE83456_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE83456_Non_normalized_data = as.matrix(GSE83456_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE83456_experimentData));GSE83456_sobject
save_raw_files(GSE83456_sobject, path = "data-raw/", geo = geo)

##### Create normalized curated assay #####
GSE83456_normed <- yr$E[,indx]
colnames(GSE83456_normed) <- ID_table$SampleID
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE83456_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)


