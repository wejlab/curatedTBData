if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

#### Read in raw data ####
geo <- "GSE107994"
sequencePlatform <- "GPL20301"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
data_list <- readDataGSE107995(geo)
probeInfo <- c("Genes", "Gene_name", "Gene_biotype")
data_list_processed <- lapply(data_list, function(x) {
  x_counts <- x %>%
    dplyr::select(-probeInfo) %>% as.matrix()
  colnames(x_counts) <- names(GEOquery::GSMList(gse))
  row.names(x_counts) <- x$Genes
  x_counts
})

data_Non_normalized_counts <- data_list_processed$data_Non_normalized
data_normalized_counts <- data_list_processed$data_normalized
#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue", "TBStatus", "Outlier",
                                         "PatientID", "TBType", "SputumSmearStatus",
                                         "Ethnicity", "BirthRegion", "uk_arrival_year",
                                         "Gender", "Age", "Visit_date", "MeasurementTime")
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$GeographicalRegion <- "UK"
characteristic_data_frame$Gender <- ifelse(characteristic_data_frame$Gender == "M",
                                           "Male","Female")
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)

MeasurementTime <- as.character(characteristic_data_frame$MeasurementTime)
index <- grep("Baseline|Active", MeasurementTime) # No space between |
MeasurementTime[-index] <- paste(MeasurementTime[-index],"months")
index1 <- grep("Active", MeasurementTime)
timepoint <- paste(sapply(strsplit(MeasurementTime[index1],"_"), function(x) x[1]), "months")
timepoint <- paste(timepoint,"(Active_TB_diagnosed)")
MeasurementTime[index1] <- timepoint
characteristic_data_frame$MeasurementTime <- MeasurementTime

TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus[grep("Active_TB", TBStatus)] <- "PTB"
TBStatus[grep("LTBI_Progressor", TBStatus)] <- "LTBI"
characteristic_data_frame$TBStatus <- TBStatus

Progression <- c()
for (i in 1:length(TBStatus_temp)){
  if(TBStatus_temp[i] == "Control"){
    Progression[i] <- "Control"
  }
  if(TBStatus_temp[i] == "Active_TB"){
    Progression[i] <- "TB"
  }
  if(TBStatus_temp[i] == "LTBI_Progressor"){
    Progression[i] <- "Positive"
  }
  if(TBStatus_temp[i] == "LTBI"){
    Progression[i] <- "Negative"
  }
}
characteristic_data_frame$Progression <- Progression
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- data_list$data_Non_normalized %>%
  dplyr::select(probeInfo) %>%
  S4Vectors::DataFrame()
colnames(new_row_data)[1:2] <- c("ID_REF", "SYMBOL_NEW")

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Akul Singhania",
                      lab = "The Francis Crick Institute",
                      contact = "akul.singhania@crick.ac.uk",
                      title = "A modular transcriptional signature identifies phenotypic heterogeneity of human tuberculosis infection.",
                      abstract = "Whole blood transcriptional signatures distinguishing active tuberculosis patients from asymptomatic latently infected individuals exist. Consensus has not been achieved regarding the optimal reduced gene sets as diagnostic biomarkers that also achieve discrimination from other diseases. Here we show a blood transcriptional signature of active tuberculosis using RNA-Seq, confirming microarray results, that discriminates active tuberculosis from latently infected and healthy individuals, validating this signature in an independent cohort. Using an advanced modular approach, we utilise the information from the entire transcriptome, which includes overabundance of type I interferon-inducible genes and underabundance of IFNG and TBX21, to develop a signature that discriminates active tuberculosis patients from latently infected individuals or those with acute viral and bacterial infections. We suggest that methods targeting gene selection across multiple discriminant modules can improve the development of diagnostic biomarkers with improved performance. Finally, utilising the modular approach, we demonstrate dynamic heterogeneity in a longitudinal study of recent tuberculosis contacts.",
                      url = "10.1038/s41467-018-04579-w",
                      pubMedIds = "29921861",
                      other=list(Platform = "Illumina HiSeq 4000 (Homo sapiens) (GPL20301)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Non_normalized_counts = as.matrix(data_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(data_normalized_counts, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg38"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))

