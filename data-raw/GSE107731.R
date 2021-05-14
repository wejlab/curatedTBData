if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

geo <- "GSE107731"
sequencePlatform <- "GPL15207"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_cel <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_cel,temp)
utils::untar(temp,exdir = tempd)

celFiles <- list.files(path = tempd, pattern = "*.CEL", full.names=TRUE)
data.affy <- affy::ReadAffy(filenames = celFiles)
GSE107731_normalized_rma <- Biobase::exprs(affy::rma(data.affy))
colnames(GSE107731_normalized_rma) <- gsub("_.*", "", colnames(GSE107731_normalized_rma))
GSE107731_normalized_data <- GSE107731_normalized_rma

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "Tissue", "PatientID",
                                         "PatientRelationship")
characteristic_data_frame$TBStatus <- ifelse(characteristic_data_frame$TBStatus == "normal",
                                             "Control", "PTB")

characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
GPL15207 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
GPL15207_annotate <- GPL15207@dataTable@table
GSE107731_probeNames <- data.frame(row.names(GSE107731_normalized_data))
colnames(GSE107731_probeNames) <- c("ID_REF")

new_row_data <- GSE107731_probeNames %>%
  dplyr::left_join(GPL15207_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$Gene.Symbol

##### Create Metadata #####
GSE107731_experimentData <- methods::new("MIAME",
                                        name = "Batudeligen Duan",
                                        lab = "Institute of clinical pharmacology of traditional Mongolian Medicine",
                                        contact = "bt8151290@sohu.com",
                                        title = "Expression profiling of blood from pulmonary tuberculosis patient and healthy controls",
                                        abstract = "Pulmonary tuberculosis is a multigene disease, and some of the genes affect the development of Pulmonary tuberculosis. The study wants to find different expression genes in blood from Pulmonary tuberculosis patient and normal people who have genetic relationship wtih each other.",
                                        url = "NA",
                                        pubMedIds = "NA",
                                        other = list(Platform = "[PrimeView] Affymetrix Human Gene Expression Array (GPL15207)"))
GSE107731_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE107731_Normalized_data = as.matrix(GSE107731_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE107731_experimentData));GSE107731_sobject
save_raw_files(GSE107731_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE107731_normalized_data,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))



