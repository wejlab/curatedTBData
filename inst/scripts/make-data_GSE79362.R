if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
#### Read in raw data ####
geo <- "GSE79362"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[1]), temp)
GSE79362_Non_normalized1 <- readxl::read_excel(temp, sheet = 1) # Training
GSE79362_Non_normalized2 <- readxl::read_excel(temp, sheet = 2) # Testing
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
row_info <- c("entry", "strand","chr","start","end","gene")
GSE79362_Non_normalized2_reduce <- GSE79362_Non_normalized2 %>%
  dplyr::select(-row_info[-1])

GSE79362_Non_normalized <- GSE79362_Non_normalized1 %>%
  dplyr::inner_join(GSE79362_Non_normalized2_reduce, by = c("entry" = "entry"))
GSE79362_Non_normalized_counts <- GSE79362_Non_normalized %>%
  dplyr::select(-row_info) %>%
  as.matrix()
row.names(GSE79362_Non_normalized_counts) <- GSE79362_Non_normalized$entry
# convert "X81_L7.LB16" to "81_L7.LB16"
colnames(GSE79362_Non_normalized_counts)[ncol(GSE79362_Non_normalized_counts)] <- sub(".*?X", "",
                                                colnames(GSE79362_Non_normalized_counts)[ncol(GSE79362_Non_normalized_counts)])
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
# Convert Non-progressor (82_L2.LB19) to 82_L2.LB19
description_id <- gsub(".*[(] *(.*?) *[)].*", "\\1", description_id_raw)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)
indx <- match(ID_table$DescriptionID, colnames(GSE79362_Non_normalized_counts))
GSE79362_Non_normalized_counts_final <- GSE79362_Non_normalized_counts[, indx]
colnames(GSE79362_Non_normalized_counts_final) <- ID_table$SampleID

#### Create Column Data ####
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
indx_diff <- which(sapply(data_characteristic,length) != length(data_characteristic[[1]]))
data_characteristic[[indx_diff]][8] <- c("qft: NA")
data_characteristic[[indx_diff]] <- c(data_characteristic[[indx_diff]], "tissue: blood")
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, "[[", x))
characteristic_data_frame <- sub("(.*?): ", "",characteristic_table) %>%
  S4Vectors::DataFrame()
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))
colnames(characteristic_data_frame) <- c("TBStatus", "Bin", "Age", "Gender",
                                         "Ethnicity", "PreviousTB", "TST", "QFT_GIT",
                                         "Tissue")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
TBStatus <- ifelse(TBStatus_temp == "case (TB progressor)", "PTB", "LTBI")
characteristic_data_frame$TBStatus <- TBStatus
Progression <- ifelse(TBStatus_temp == "case (TB progressor)", "Positive", "Negative")
characteristic_data_frame$Progression <- Progression
characteristic_data_frame$Gender <- ifelse(characteristic_data_frame$Gender == "male", "Male", "Female")
qft_temp <- qft <- as.character(characteristic_data_frame$QFT_GIT)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
QFT_GIT <- firstup(qft)
QFT_GIT[QFT_GIT == "NA"] <- NA
characteristic_data_frame$QFT_GIT <- QFT_GIT
characteristic_data_frame$GeographicalRegion <- "South Africa"
characteristic_data_frame$Tissue <- "Whole Blood"
hist_TB <- ifelse(characteristic_data_frame$PreviousTB == "n", "No", "Yes")
characteristic_data_frame$PreviousTB <- hist_TB
Gender <- characteristic_data_frame$Gender
characteristic_data_frame$Gender = ifelse(Gender == "female", "Female", "Male")
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
characteristic_data_frame$TST <- as.numeric(characteristic_data_frame$TST)
# Download patient information online
urlMeta <- c("https://ars.els-cdn.com/content/image/1-s2.0-S0140673615013161-mmc2.xlsx")
tempMeta <- tempfile()
download.file(urlMeta, tempMeta)
GSE79362Metadata1 <- readxl::read_excel(tempMeta, sheet = "SupTab6_RNASeqMetadata")
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
colnames(GSE79362Metadata1) <- GSE79362Metadata1[2,]
GSE79362Metadata_sub <- GSE79362Metadata1[-c(1:2), c(1:3, 5:6)]
colnames(GSE79362Metadata_sub) <- c("PatientID", "MeasurementTime", "TimeToTB","DemographicalBins",
                                    "ACS_cohort")
# The bin number in the metadata is not the same with bin on the GEO website, use GEO as the reference
MeasurementTime <- GSE79362Metadata_sub$MeasurementTime
MeasurementTime[grep("IC", MeasurementTime)] <- "End"
MeasurementTime <- gsub("D", " Day(s)", MeasurementTime)
GSE79362Metadata_sub$MeasurementTime <- MeasurementTime

TimeToTB <- TimeToTB_temp <- GSE79362Metadata_sub$TimeToTB
for (i in 1:length(TimeToTB_temp)) {
  if (TimeToTB_temp[i] != "---" && !is.na(TimeToTB_temp[i])) {
    TimeToTB[i] <- paste0(TimeToTB[i], " Day(s)")
  }
}
GSE79362Metadata_sub$TimeToTB <- TimeToTB
characteristic_data_frame_new <- cbind(characteristic_data_frame, GSE79362Metadata_sub)
col_info <- create_standard_coldata(characteristic_data_frame_new)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Rwo Data #####
new_row_data <- GSE79362_Non_normalized %>%
  dplyr::select(row_info)
colnames(new_row_data)[1] <- "ID_REF"
colnames(new_row_data)[6] <- "SYMBOL_NEW"

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Daniel Edward Zak",
                      lab = "Seattle Biomedical Research Institute",
                      contact = "Willem.hanekom@gatesfoundation.org",
                      title = "A blood RNA signature for tuberculosis disease risk: a prospective cohort study.",
                      abstract = "In this prospective cohort study, we followed up healthy, South African adolescents aged 12–18 years from the adolescent cohort study (ACS) who were infected with M tuberculosis for 2 years. We collected blood samples from study participants every 6 months and monitored the adolescents for progression to tuberculosis disease.",
                      url = "10.1016/S0140-6736(15)01316-1",
                      pubMedIds = "27017310",
                      other=list(Platform = "Illumina HiSeq 2000 (Homo sapiens) (GPL11154)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE79362_Non_normalized_counts_final)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)

##### Create curated matched assay #####
GSE79362_normalized_counts <- probesetsToGenes(new_row_data,
                                               GSE79362_Non_normalized_counts_final,
                                               median)
saveRDS(GSE79362_normalized_counts, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))
