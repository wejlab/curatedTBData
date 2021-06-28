if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
geo <- "GSE107104"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[1]), temp)
GSE107104_normalized_raw <- readxl::read_excel(temp)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

GSE107104_normalized_counts <- as.matrix(GSE107104_normalized_raw[,-c(1,2)])
row.names(GSE107104_normalized_counts) <- GSE107104_normalized_raw$id
colnames(GSE107104_normalized_counts) <- names(GEOquery::GSMList(gse))

###### Create Row Data #####
new_row_data <- GSE107104_normalized_raw %>%
  dplyr::select(c("Symbol","id")) %>%
  S4Vectors::DataFrame()
colnames(new_row_data) <- c("SYMBOL_NEW","ID_REF")

#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue", "TBStatus","Treatment")
HIVStatus <- rep("Positive", nrow(characteristic_data_frame))
TBStatus <- TBStatus_temp <- characteristic_data_frame$TBStatus
unique(TBStatus_temp)
for (i in 1:length(TBStatus)) {
  if (TBStatus_temp[i] == "HIV and TB") {
    TBStatus[i] <- "PTB"
    HIVStatus[i] <- "Positive"
  } else if (TBStatus_temp[i] == "HIV") {
    TBStatus[i] <- "Control"
    HIVStatus[i] <- "Positive"
  } else if (TBStatus_temp[i] == "TB-HIV") {
    TBStatus[i] <- "PTB"
    HIVStatus[i] <- "Negative"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$HIVStatus <- HIVStatus


characteristic_data_frame$HIVStatus <- rep("Positive", nrow(characteristic_data_frame))
characteristic_data_frame$TBStatus <- c(rep("PTB", 15), rep("Control", 16), rep("PTB", 2))
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$GeographicalRegion <- "Uganda"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)
##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Padmini Salgame",
                      lab = "Rutgers New Jersey Medical School",
                      contact = "padmini.salgame@rutgers.edu",
                      title = "Tuberculosis in advanced HIV infection is associated with increased expression of IFNγ and its downstream targets.",
                      abstract = "Tuberculosis (TB) is the major cause of death in HIV-infected individuals. However, diagnosis of TB in HIV remains challenging. The aim of this study was to determine the performance of published gene signatures in distinguishing active TB in advanced HIV. We concluded that gene expression of FcGR1A and BATF2, and plasma protein levels of IFNγ and CXCL10 can singly classify active TB in HIV-infected with advanced HIV.",
                      url = "10.1016/j.tube.2020.101898",
                      pubMedIds = "29764370",
                      other=list(Platform = "Illumina HiSeq 2000 (Homo sapiens) (GPL11154)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE107104_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)

##### Create curated matched assay #####
GSE107104_normalized_counts_matched <- probesetsToGenes(new_row_data,
                                                        GSE107104_normalized_counts,
                                                        median)
saveRDS(GSE107104_normalized_counts_matched, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))



