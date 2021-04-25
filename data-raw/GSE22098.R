if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
##### Read in raw data #####
geo <- "GSE22098"
sequencePlatform <- "GPL6947"
GSE22098_Non_normalized_list_noPvalue1 <- readRawData(geo, sequencePlatform)
# Take long time
GSE22098_Non_normalized_list_noPvalue <- lapply(GSE22098_Non_normalized_list_noPvalue1, function(x) {
    # Remove duplicates in ID_REF
    x %>% dplyr::as_tibble() %>%
          dplyr::group_by(ID_REF) %>%
          dplyr::summarise_all(median)})
# Merge list to a matrix based on the ID_REF
GSE22098_Non_normalized <- Reduce(function(x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE22098_Non_normalized_list_noPvalue, function(x) {x}))

row.names(GSE22098_Non_normalized) <- GSE22098_Non_normalized$ID_REF
GSE22098_Non_normalized_counts <- GSE22098_Non_pvalue <- GSE22098_Non_normalized[-1]

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Age", "Gender", "Ethnicity", "HealthControl")
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))

characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
characteristic_data_frame$DiabetesStatus <- "Negative"
characteristic_data_frame$HIVStatus <- "Negative"
data_title <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title) %>% unlist
Notes <- sapply(strsplit(data_title,"-"),function(x) x[2])
# Based on information of different diseases
characteristic_data_frame$TBStatus <-  ifelse(Notes == "H", "Control", "OD")
StillStatus <- ifelse(Notes == "Still", "Positive", "Negative")
AdultSLE_Status <- ifelse(Notes == "ASLE", "Positive", "Negative")
PaediatricSLE_Status <- ifelse(Notes == "PSLE", "Positive", "Negative")
StaphStatus <- ifelse(Notes == "Staph", "Positive", "Negative")
StrepStatus <- ifelse(Notes == "Strep", "Positive", "Negative")

characteristic_data_frame$StillStatus <- StillStatus
characteristic_data_frame$AdultSLE_Status <- AdultSLE_Status
characteristic_data_frame$PaediatricSLE_Status <-  PaediatricSLE_Status
characteristic_data_frame$StaphStatus <- StaphStatus
characteristic_data_frame$StrepStatus <- StrepStatus
characteristic_data_frame$Tissue <- "Whole Blood"

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE22098_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)
##### Create Metadata #####
GSE22098_experimentData <- methods::new("MIAME",
                                        name = "Damien Chaussabel",
                                        lab = "Baylor Institute for Immunology Research",
                                        contact = "DChaussabel@benaroyaresearch.org",
                                        title = "An interferon-inducible neutrophil-driven blood transcriptional signature in human tuberculosis",
                                        abstract = "Three milliliters of whole blood was collected in Tempus tubes from 12 pediatric streptococcus, 40 pediatric staphylococcus, 31 still’s disease, 82 pediatric systemic lupus erythematosus (SLE) and 28 adult SLE patients. RNA was extracted and globin reduced.",
                                        url = "10.1038/nature09247",
                                        pubMedIds = "20725040",
                                        other = list(Platform = "Illumina HumanHT-12 V3.0 expression beadchip (GPL6947)"))
GSE22098_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE22098_Non_normalized_counts= as.matrix(GSE22098_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE22098_experimentData));GSE22098_sobject
save_raw_files(GSE22098_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
