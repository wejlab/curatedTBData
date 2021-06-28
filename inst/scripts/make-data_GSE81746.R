if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

geo <- "GSE81746"
sequencePlatform <- "GPL17077"
temp <- tempfile()
tempd <- tempdir()
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
utils::download.file(as.character(urls$url[1]), temp)
utils::untar(temp, exdir = tempd)
filesPath <- list.files(tempd, pattern = "GSM.*", full.names = TRUE)
GSE81746_raw <- limma::read.maimages(filesPath, source = "agilent", green.only = TRUE)
GSE81746_Non_normalized_data <- GSE81746_raw$E
row.names(GSE81746_Non_normalized_data) <- GSE81746_raw$genes$ProbeName
col_name1 <- gsub(".*/", "", colnames(GSE81746_Non_normalized_data))
colnames(GSE81746_Non_normalized_data) <- gsub("_.*", "", col_name1)
GSE81746_Non_pvalue <- GSE81746_Non_normalized_data

##### Create curated assay ####
curatedExprs <- norm_probeToGenes_Agilent(GSE81746_raw, FUN = median)
colnames(curatedExprs) <- colnames(GSE81746_Non_normalized_data)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
data_characteristic[[5]] <- c("subject status: tuberculosis (TB) patient",
                              "age: NA", "gender: NA", "tissue: whole blood")
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, "[[", x))

characteristic_data_frame <- sub("(.*?): ", "", characteristic_table) %>%
  S4Vectors::DataFrame()
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))
colnames(characteristic_data_frame) <- c("TBStatus", "Age", "Gender", "Tissue")
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$TBStatus <- ifelse(characteristic_data_frame$TBStatus == "tuberculosis (TB) patient",
                                             "PTB", "Control")
age <- characteristic_data_frame$Age
age1 <- unlist(strsplit(age,"y") %>% unlist())
age1[age1 == "NA"] <- NA
characteristic_data_frame$Age <- as.numeric(age1)

sex <- characteristic_data_frame$Gender
for(i in 1:length(sex)){
  if(sex[i] == "male"){sex[i] <-  "Male"}
  if(sex[i] == "female"){sex[i] <-  "Female"}
  if(sex[i] == "NA"){sex[i] <-  NA}
}
characteristic_data_frame$Gender <- sex
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
GPL17077 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
GPL17077_annotate <- GPL17077@dataTable@table

GSE81746_probeNames <- data.frame(row.names(GSE81746_Non_pvalue))
colnames(GSE81746_probeNames) <- c("ID_REF")

new_row_data <- GSE81746_probeNames %>%
  dplyr::left_join(GPL17077_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$GENE_SYMBOL

##### Create Metadata #####
GSE81746_experimentData <- methods::new("MIAME",
                                        name = "Nagasuma Chandra",
                                        lab = "Indian Institute of Science",
                                        contact = "nchandra@biochem.iisc.ernet.in",
                                        title = "Gene expression profiling of tuberculosis patients from India",
                                        abstract = "Whole genome microarray expression profiling was employed to identify differential gene expression profiles characteristic of tuberculosis patients in the South-Indian cohort.",
                                        url = "NA",
                                        pubMedIds = "NA",
                                        other = list(Platform = "Agilent-039494 SurePrint G3 Human GE v2 8x60K Microarray 039381 (Probe Name version) (GPL17077)"))
GSE81746_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE81746_Non_normalized_data = as.matrix(GSE81746_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE81746_experimentData));GSE81746_sobject
save_raw_files(GSE81746_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
