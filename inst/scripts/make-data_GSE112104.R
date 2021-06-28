if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

#### Read in raw data ####
geo <- "GSE112104"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)

url_raw <- as.character(urls$url[1])
temp <- tempfile()
utils::download.file(url_raw, temp)

GSE112104_normalized_raw <- read.csv(temp,row.names = 1)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
title <- sapply(GEOquery::GSMList(gse), function(x) x@header$title)
indx <- match(gsub(":.*","", title), colnames(GSE112104_normalized_raw))
GSE112104_normalized_counts <- GSE112104_normalized_raw[,indx]
colnames(GSE112104_normalized_counts) <- names(GEOquery::GSMList(gse))

#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <-  c("Gender", "Age", "TimeToTB", "Tissue")
characteristic_data_frame$GeographicalRegion <- "Brazil"
title_sub <- gsub(".*: ", "", title)
Progression <- rep("Negative", nrow(characteristic_data_frame))
Progression[grep("Progressor", title_sub)] <- "Positive"
Progression[grep("TB patient", title_sub)] <- NA
characteristic_data_frame$Progression <- Progression

TBStatus <- rep("LTBI", nrow(characteristic_data_frame))
TBStatus[grep("TB patient", title_sub)] <- "PTB"
characteristic_data_frame$TBStatus <- TBStatus
Gender <- ifelse(characteristic_data_frame$Gender == "female", "Female", "Male")
characteristic_data_frame$Gender <- Gender
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)

timeToTB <- as.character(characteristic_data_frame$TimeToTB)
timeToTB[which(timeToTB == "NA")] <- NA
timeToTB[-which(is.na(timeToTB))] <- paste(timeToTB[-which(is.na(timeToTB))], "days")
characteristic_data_frame$TimeToTB <- timeToTB
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(GSE112104_normalized_counts),
                                     SYMBOL_NEW = row.names(GSE112104_normalized_counts))

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Padmini Salgame",
                      lab = "Rutgers New Jersey Medical School",
                      contact = "padmini.salgame@rutgers.edu",
                      title = "Cross-validation of existing signatures and derivation of a novel 29-gene transcriptomic signature predictive of progression to TB in a Brazilian cohort of household contacts of pulmonary TB.",
                      abstract = "The goal of this study was to identify individuals at risk of progression and reactivation among household contacts (HHC) of pulmonary TB cases in Vitoria, Brazil. We first evaluated the predictive performance of six published signatures on the transcriptional dataset obtained from peripheral blood mononuclear cell samples from HHC that either progressed to TB disease or not (non-progressors) during a five-year follow-up. The area under the curve (AUC) values for the six signatures ranged from 0.670 to 0.461, and the PPVs did not reach the WHO published target product profiles (TPPs). We therefore used as training cohort the earliest time-point samples from the African cohort of adolescents (GSE79362) and applied an ensemble feature selection pipeline to derive a novel 29-gene signature (PREDICT29). PREDICT29 was tested on 16 progressors and 21 non-progressors. PREDICT29 performed better in segregating progressors from non-progressors in the Brazil cohort with the area under the curve (AUC) value of 0.911 and PPV of 20%. This proof of concept study demonstrates that PREDICT29 can predict risk of progression/reactivation to clinical TB disease in recently exposed individuals at least 5 years prior to disease development. Upon validation in larger and geographically diverse cohorts, PREDICT29 can be used to risk-stratify recently infected for targeted therapy.",
                      url = "10.1016/j.tube.2020.101898",
                      pubMedIds = "32090859",
                      other=list(Platform = "Illumina HiSeq 2500 (Homo sapiens) (GPL16791)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE112104_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(GSE112104_normalized_counts, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))

