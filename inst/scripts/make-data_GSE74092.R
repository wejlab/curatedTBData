if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE74092"
sequencePlatform <- "GPL21040"
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_sub <- as.character(urls$url[1])
temp <- tempfile()
utils::download.file(url_sub, temp)

GSE74092_Non_normalized <- read.delim(gzfile(temp), stringsAsFactors = FALSE)
# Remove columns with unknown sample names, remove rows with empty string
GSE74092_raw <- GSE74092_Non_normalized %>%
  dplyr::filter(GSE74092_Non_normalized[,1] != "") %>%
  dplyr::select_if(!is.na(GSE74092_Non_normalized[2, ]) == TRUE)
colnames(GSE74092_raw) <- as.character(GSE74092_raw[1, ])
GSE74092_raw_final <- GSE74092_raw[-1,]
GSE74092_Non_normalized_data <- GSE74092_raw_final %>%
  dplyr::select(-c("ID_REF", "SPOT_ID")) %>% as.matrix()
mode(GSE74092_Non_normalized_data) <- "numeric"
row.names(GSE74092_Non_normalized_data) <- GSE74092_raw_final$SPOT_ID

gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id_raw)
indx <- match(ID_table$DescriptionID, colnames(GSE74092_Non_normalized_data))
GSE74092_Non_normalized_data <- GSE74092_Non_normalized_data[,indx]
colnames(GSE74092_Non_normalized_data) <- ID_table$SampleID

##### Create Column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("SkinTest", "TBStatus", "Tissue",
                                         "ActivationAgent", "Time")
# ActivationAgent and Time are all NA, remove them
characteristic_data_frame <- characteristic_data_frame[, -c(4,5)]
characteristic_data_frame$Tissue <- "Whole Blood"

TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus <- ifelse(TBStatus_temp=='TB', 'PTB', "Control")
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$GeographicalRegion <- "India"
SkinTest <- SkinTest_temp <- as.character(characteristic_data_frame$SkinTest)
SkinTest[grep("GeneXpert-pos",SkinTest_temp)] <- NA
SkinTest[grep("TSTneg",SkinTest_temp)] <- "Negative"
SkinTest[grep("TSTpos",SkinTest_temp)] <- "Positive"
characteristic_data_frame$SkinTest <- SkinTest

xpert <- rep(NA, nrow(characteristic_data_frame))
xpert[grep("GeneXpert-pos",SkinTest_temp)] <- "Positive"
characteristic_data_frame$GenXpert <- xpert
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
GPL21040 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
GPL21040_annotate <- GPL21040@dataTable@table
GSE74092_probeNames <- data.frame(row.names(GSE74092_Non_normalized_data))
colnames(GSE74092_probeNames) <- c("ID_REF")

new_row_data <- GSE74092_probeNames %>%
  dplyr::left_join(GPL21040_annotate, by = c("ID_REF" = "WELL")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$ORF

##### Create Metadata #####
GSE74092_experimentData <- methods::new("MIAME",
                                        name = "Jeroen Maertzdorf",
                                        lab = "MPIIB",
                                        contact = "maertzdorf@mpiib-berlin.mpg.de",
                                        title = "Concise gene signature for point-of-care classification of tuberculosis.",
                                        abstract = "There is an urgent need for new tools to combat the ongoing tuberculosis (TB) pandemic. Gene expression profiles based on blood signatures have proved useful in identifying genes that enable classification of TB patients, but have thus far been complex. Using real-time PCR analysis, we evaluated the expression profiles from a large panel of genes in TB patients and healthy individuals in an Indian cohort. Classification models were built and validated for their capacity to discriminate samples from TB patients and controls within this cohort and on external independent gene expression datasets. A combination of only four genes distinguished TB patients from healthy individuals in both cross-validations and on separate validation datasets with very high accuracy. An external validation on two distinct cohorts using a real-time PCR setting confirmed the predictive power of this 4-gene tool reaching sensitivity scores of 88% with a specificity of around 75%. Moreover, this gene signature demonstrated good classification power in HIV(+) populations and also between TB and several other pulmonary diseases. Here we present proof of concept that our 4-gene signature and the top classifier genes from our models provide excellent candidates for the development of molecular point-of-care TB diagnosis in endemic areas.",
                                        url = "10.15252/emmm.201505790",
                                        pubMedIds = "26682570",
                                        other = list(Platform = "Qiagen custom RT2 Profiler Array"))
GSE74092_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE74092_Non_normalized_data = as.matrix(GSE74092_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE74092_experimentData));GSE74092_sobject
save_raw_files(GSE74092_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
GSE74092_normalized_data <- limma::normalizeBetweenArrays(GSE74092_Non_normalized_data,
                                                          method = "quantile")
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE74092_normalized_data,
                                 FUN = median)
# Impute missing values in the exprs
GSE74092_impute <- impute::impute.knn(curatedExprs)
GSE74092_impute_dat <- GSE74092_impute$data
saveRDS(GSE74092_impute_dat, paste0("data-raw/", geo, "_assay_curated.RDS"))


