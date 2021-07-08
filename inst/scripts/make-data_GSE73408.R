if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE73408"
sequencePlatform <- "GPL11532"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_cel <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_cel, temp)
utils::untar(temp, exdir = tempd)

celFiles <- list.files(path = tempd, pattern = "*.CEL", full.names = TRUE)
GSE73408_data.oligo <- oligo::read.celfiles(celFiles)

GSE73408_normalized_rma <- Biobase::exprs(oligo::rma(GSE73408_data.oligo)) # A matrix
colnames(GSE73408_normalized_rma) <- gsub("_.*", "", colnames(GSE73408_normalized_rma))
GSE73408_normalized_data <- GSE73408_normalized_rma

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue", "Gender", "TBStatus",
                                         "Ethnicity", "Age", "BirthRegion",
                                         "DiabetesStatus", "SmokingStatus",
                                         "BCG")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)

TBStatus[which(TBStatus_temp == "TB")] <- "PTB"
TBStatus[which(TBStatus_temp == "PNA")] <- "OD"
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$GeographicalRegion <- "USA"
characteristic_data_frame$Tissue <- "Whole Blood"
BirthRegion <- characteristic_data_frame$BirthRegion
BirthRegion[grep("mexico", BirthRegion)] <- "Mexico"
characteristic_data_frame$BirthRegion <- BirthRegion
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
DiabetesStatus <- rep("Negative", nrow(characteristic_data_frame))
DiabetesStatus[grep("yes",characteristic_data_frame$DiabetesStatus)] <- "Positive"
characteristic_data_frame$DiabetesStatus <- DiabetesStatus
SmokingStatus <- ifelse(as.character(characteristic_data_frame$SmokingStatus) == "yes",
                        "Positive", "Negative")
characteristic_data_frame$SmokingStatus <- SmokingStatus
bcg <- as.character(characteristic_data_frame$BCG)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
characteristic_data_frame$BcgVaccinated <- firstup(bcg)
characteristic_data_frame <- characteristic_data_frame[,-which(colnames(characteristic_data_frame) == "BCG")]
PneumoniaStatus <- rep("Negative", nrow(characteristic_data_frame))
PneumoniaStatus[which(characteristic_data_frame$Notes == "PNA")] <- "Positive"
characteristic_data_frame$PneumoniaStatus <- PneumoniaStatus

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
GPL11532 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
GPL11532_annotate <- GPL11532@dataTable@table
GPL11532_annotate$ID <- as.character(GPL11532_annotate$ID)
GSE73408_probeNames <- data.frame(row.names(GSE73408_normalized_data))
colnames(GSE73408_probeNames) <- c("ID_REF")

new_row_data <- GSE73408_probeNames %>%
  dplyr::left_join(GPL11532_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
SYMBOL_ori <- as.character(new_row_data$gene_assignment)
SYMBOL_ori <- ifelse(SYMBOL_ori == "---", "NA", SYMBOL_ori)
index <- which(SYMBOL_ori != "NA")
SYMBOL <- strsplit(SYMBOL_ori[index], "///")

SYMBOL_new <- lapply(1:length(SYMBOL),
                     function(i) lapply(strsplit(SYMBOL[[i]], "//"), "[[", 2) %>%
                       unlist() %>% unique())
SYMBOL_collapse <- sapply(SYMBOL_new, function(x) gsub(" ", "", paste(x, collapse = "///")))
SYMBOL_ori[index] <- SYMBOL_collapse

new_row_data$SYMBOL_NEW <- SYMBOL_ori

##### Create Metadata #####
GSE73408_experimentData <- methods::new("MIAME",
                                        name = "Nicholas Walter",
                                        lab = "University of Colorado",
                                        contact = "nicholas.walter@ucdenver.edu",
                                        title = "Blood Transcriptional Biomarkers for Active Tuberculosis among Patients in the United States: a Case-Control Study with Systematic Cross-Classifier Evaluation.",
                                        abstract = "Whole blood from US patients was collected to accurately diagnose active TB in racially and ethnically diverse populations for RNA extraction and hybridization on Affymetrix microarrays.",
                                        url = "10.1128/JCM.01990-15",
                                        pubMedIds = "26582831",
                                        other = list(Platform = "Affymetrix Human Gene 1.1 ST Array (GPL11532)"))
GSE73408_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE73408_Normalized_data = as.matrix(GSE73408_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE73408_experimentData));GSE73408_sobject
save_raw_files(GSE73408_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE73408_normalized_data,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
