if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE19443"
sequencePlatform <- "GPL6947"
GSE19443_Non_normalized_list_noPvalue <- readRawData(geo, sequencePlatform)
# Merge list to a matrix based on the ID_REF
GSE19443_Non_normalized <- Reduce(function(x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE19443_Non_normalized_list_noPvalue, function(x) {x}))
row.names(GSE19443_Non_normalized) <- GSE19443_Non_normalized$ID_REF

GSE19443_Non_normalized_counts <- GSE19443_Non_pvalue <- GSE19443_Non_normalized[-1]
# dim(GSE19443_Non_normalized_counts) should be 48803*44

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Age", "Gender", "Ethnicity", "TBStatus",
                                         "GeographicalRegion", "BcgVaccinated",
                                         "BirthRegion", "TST", "exposure_latent",
                                         "index_case_disease_site", "smear_of_index_case",
                                         "modal_x_ray_grade", "SputumSmearStatus",
                                         "sputum_culture", "bal_smear", "bal_culture",
                                         "isolate_sensitivity")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus <- ifelse(TBStatus_temp == "Control (BCG+)", "Control", TBStatus_temp)
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
# Replace n.a. or empty string with NA
TST_temp <- characteristic_data_frame$TST
TST <- as.numeric(ifelse(TST_temp == "n.a", NA, TST_temp))
characteristic_data_frame$TST <- TST
characteristic_data_frame$exposure_latent <- NA
info_sub <- c("exposure_latent", "index_case_disease_site", "smear_of_index_case",
              "modal_x_ray_grade", "SputumSmearStatus", "sputum_culture", "bal_smear",
              "bal_culture", "isolate_sensitivity")
for(col_name in info_sub) {
  info <- characteristic_data_frame[, col_name]
  info1 <- ifelse(info == "n.a", NA, info)
  info2 <- ifelse(info1 == "", NA, info1)
  characteristic_data_frame[, col_name] <- info2
}

# Add addtional metadata information from the paper
characteristic_data_frame$DiabetesStatus <- "Negative"
characteristic_data_frame$HIVStatus <- "Negative"
sputum_culture <- gsub("tuberculosis", "M.tuberculosis",
                       characteristic_data_frame$sputum_culture)
characteristic_data_frame$sputum_culture <- sputum_culture
bal_culture <- gsub("tuberculosis", "M.tuberculosis",
                    characteristic_data_frame$bal_culture)
characteristic_data_frame$bal_culture <- bal_culture

title <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title) %>% unlist
title_split <- strsplit(title,"_")
PatientID <- sapply(title_split, function(x) x[3])
characteristic_data_frame$PatientID <- PatientID
Tissue <- sapply(title_split, function(x) x[5])
Tissue[grep("Neut", Tissue)] <- "Neutrophils"
Tissue[grep("Mono", Tissue)] <- "Monocytes"
characteristic_data_frame$Tissue <- Tissue

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data information #####
row_data <- map_gene_symbol(GSE19443_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)
##### Create Metadata #####
GSE19443_experimentData <- methods::new("MIAME",
                               name = "Damien Chaussabel",
                               lab = "Baylor Institute for Immunology Research",
                               contact = "DChaussabel@benaroyaresearch.org",
                               title = "An interferon-inducible neutrophil-driven blood transcriptional signature in human tuberculosis",
                               abstract = "Whole blood collected in EDTA tubes from patients with active TB disease and healthy controls. Blood was then processed or separated sequentially into neutrophil, monocyte, CD4+ or CD8+ populations and then processed. All patients were sampled prior to the initiation of any antimycobacterial therapy.",
                               url = "10.1038/nature09247",
                               pubMedIds = "20725040",
                               other = list(Platform = "Illumina HumanHT-12 V3.0 expression beadchip (GPL6947)"))
GSE19443_sobject <- SummarizedExperiment(
  assays = list(GSE19443_Non_normalized_counts= as.matrix(GSE19443_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE19443_experimentData))
save_raw_files(GSE19443_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
