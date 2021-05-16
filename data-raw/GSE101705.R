if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
geo <- "GSE101705"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[1]), temp)
GSE101705_normalized_raw <- readxl::read_excel(temp)
colnames(GSE101705_normalized_raw)[1] <- "ID_REF"
GSE101705_normalized_counts <- as.matrix(GSE101705_normalized_raw[, -1])
row.names(GSE101705_normalized_counts) <- GSE101705_normalized_raw$ID_REF

row.names(GSE101705_normalized_counts) <- GSE101705_normalized_raw$ID_REF
GSE101705_normalized_counts_avg <- GSE101705_normalized_raw %>%
  dplyr::group_by(ID_REF) %>%
  dplyr::summarise_all(mean) %>% data.frame()
GSE101705_normalized_counts_avg1 <- as.matrix(GSE101705_normalized_counts_avg[, -1])
row.names(GSE101705_normalized_counts_avg1) <- GSE101705_normalized_counts_avg$ID_REF

description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID=description_id_raw)
indx <- match(ID_table$DescriptionID, colnames(GSE101705_normalized_counts))
GSE101705_normalized_counts <- GSE101705_normalized_counts[, indx]
colnames(GSE101705_normalized_counts) <- ID_table$SampleID
GSE101705_normalized_counts_avg1 <- GSE101705_normalized_counts_avg1[, indx]
colnames(GSE101705_normalized_counts_avg1) <- ID_table$SampleID

#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- "TBStatus"
characteristic_data_frame$TBStatus <- ifelse(characteristic_data_frame$TBStatus == "TB",
                                             "PTB","LTBI")
characteristic_data_frame$GeographicalRegion <- "South India"
characteristic_data_frame$PreviousTB <- "No"
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(GSE101705_normalized_counts),
                                     SYMBOL_NEW = row.names(GSE101705_normalized_counts))
##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Padmini Salgame",
                      lab = "Rutgers New Jersey Medical School",
                      contact = "padmini.salgame@rutgers.edu",
                      title = "Existing blood transcriptional classifiers accurately discriminate active tuberculosis from latent infection in individuals from south India.",
                      abstract = "Transcriptome profiles via RNAseq of TB and LTBI individuals in a South Indian cohort were used to assess the accuracy of established TB disease signatures.",
                      url = "10.1016/j.tube.2018.01.002",
                      pubMedIds = "29559120",
                      other=list(Platform = "Illumina NextSeq 500 (Homo sapiens) (GPL18573)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE101705_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(GSE101705_normalized_counts_avg1, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))


