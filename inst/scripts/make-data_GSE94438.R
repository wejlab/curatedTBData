if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

#### Read in raw data ####
geo <- "GSE94438"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[2]), temp)
GSE94438_Non_normalized <- read.csv(gzfile(temp), row.names = 1, header = TRUE)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

GSE94438_Non_normalized_counts <- GSE94438_Non_normalized %>%
  dplyr::select(-symbol) %>% as.matrix()
colnames(GSE94438_Non_normalized_counts) <- gsub(".*X", "",
                                                 colnames(GSE94438_Non_normalized_counts))
#### Create Column Data and modify raw reads ####
# only 418 out of 434 samples were included in the raw data obtained from online
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <-  c("Tissue", "Code", "PatientID", "GeographicalRegion",
                                          "Age", "Gender", "TBStatus", "TimeFromExposure",
                                          "TimeToTB")
################
characteristic_data_frame_sub <- characteristic_data_frame %>%
  data.frame() %>%
  dplyr::distinct(Code, .keep_all = TRUE)

indx <- match(characteristic_data_frame_sub$Code, colnames(GSE94438_Non_normalized_counts))
GSE94438_Non_normalized_counts <- GSE94438_Non_normalized_counts[, indx]
colnames(GSE94438_Non_normalized_counts) <- row.names(characteristic_data_frame_sub)
##################
characteristic_data_frame$Tissue <- "Whole Blood"

TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus <- rep(NA, nrow(characteristic_data_frame))
TBStatus[grep("case", TBStatus_temp)] <- "PTB"
TBStatus[grep("Control", TBStatus_temp)] <- "Control"
characteristic_data_frame$TBStatus <- TBStatus

Progression <- TBStatus
Progression[grep("PTB", TBStatus)] <- "Positive"
Progression[grep("Control", TBStatus)] <- "Negative"
characteristic_data_frame$Progression <- Progression

geo_region <- rep(NA, nrow(characteristic_data_frame))
unique(characteristic_data_frame$GeographicalRegion)
geo_region[grep("MRC", characteristic_data_frame$GeographicalRegion)] <- "The Gambia"
geo_region[grep("SUN", characteristic_data_frame$GeographicalRegion)] <- "South Africa"
geo_region[grep("AHRI", characteristic_data_frame$GeographicalRegion)] <- "Ethiopia"
characteristic_data_frame$GeographicalRegion <- geo_region
characteristic_data_frame$Gender <- ifelse(characteristic_data_frame$Gender == "M",
                                           "Male", "Female")
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)

timeToTB <- characteristic_data_frame$TimeToTB
timeToTB[timeToTB == "NA"] <- NA
index <- which(is.na(timeToTB))
timeToTB[-index] <- paste(timeToTB[-index], "month(s)")
characteristic_data_frame$TimeToTB <- timeToTB

timeFromEx <- characteristic_data_frame$TimeFromExposure
timeFromEx[timeFromEx == "NA"] <- NA
timeFromEx[-which(is.na(timeFromEx))] <- paste(timeFromEx[-which(is.na(timeFromEx))], "month(s)")
characteristic_data_frame$TimeFromExposure <- timeFromEx
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(GSE94438_Non_normalized),
                                     SYMBOL_NEW = GSE94438_Non_normalized$symbol)

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Daniel Zak",
                      lab = "Center for Infectious Disease Research",
                      contact = "danielzak000@gmail.com",
                      title = "Four-gene Pan-African Blood Signature Predicts Progression to Tuberculosis.",
                      abstract = "Rationale: Contacts of patients with tuberculosis (TB) constitute an important target population for preventive measures because they are at high risk of infection with Mycobacterium tuberculosis and progression to disease.Objectives: We investigated biosignatures with predictive ability for incident TB.Methods: In a case-control study nested within the Grand Challenges 6-74 longitudinal HIV-negative African cohort of exposed household contacts, we employed RNA sequencing, PCR, and the pair ratio algorithm in a training/test set approach. Overall, 79 progressors who developed TB between 3 and 24 months after diagnosis of index case and 328 matched nonprogressors who remained healthy during 24 months of follow-up were investigated.Measurements and Main Results: A four-transcript signature derived from samples in a South African and Gambian training set predicted progression up to two years before onset of disease in blinded test set samples from South Africa, the Gambia, and Ethiopia with little population-associated variability, and it was also validated in an external cohort of South African adolescents with latent M. tuberculosis infection. By contrast, published diagnostic or prognostic TB signatures were predicted in samples from some but not all three countries, indicating site-specific variability. Post hoc meta-analysis identified a single gene pair, C1QC/TRAV27 (complement C1q C-chain / T-cell receptor-α variable gene 27) that would consistently predict TB progression in household contacts from multiple African sites but not in infected adolescents without known recent exposure events.Conclusions: Collectively, we developed a simple whole blood-based PCR test to predict TB in recently exposed household contacts from diverse African populations. This test has potential for implementation in national TB contact investigation programs.",
                      url = "10.1164/rccm.201711-2340OC",
                      pubMedIds = "29624071",
                      other=list(Platform = "Illumina HiSeq 2000 (Homo sapiens) (GPL11154)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE94438_Non_normalized_counts)),
#  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(new_col_info, paste0("data-raw/", geo, "_column_data.RDS"))

##### Create curated matched assay #####
counts <- GSE94438_Non_normalized_counts
counts[counts<10] <- 10
counts <- log(counts, base=2)
NormFactor <- edgeR::calcNormFactors(counts, method = "TMM")
ScaleFactors <- colSums(counts) * NormFactor
counts_normalized <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))

GSE94438_normalized_counts_matched <- probesetsToGenes(new_row_data, counts_normalized,
                                                        median)
saveRDS(GSE94438_normalized_counts_matched, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))
