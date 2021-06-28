if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE19439"
sequencePlatform <- "GPL6947"
GSE19439_data_list <- readRawData(geo, sequencePlatform)
GSE19439_Non_normalized_data <- GSE19439_Non_pvalue <- GSE19439_data_list$data_Non_normalized

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
for (i in 1:length(TBStatus_temp)) {
  if (TBStatus_temp[i] == "Control (BCG+)" || TBStatus_temp[i] == "Control (BCG-)") {
    TBStatus[i] <- "Control"
  } else if (TBStatus_temp[i] == "Latent") {
    TBStatus[i] <- "LTBI"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$DiabetesStatus <- "Negative"
characteristic_data_frame$HIVStatus <- "Negative"
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
# Replace n.a. or empty string with NA
TST_temp <- characteristic_data_frame$TST
TST <- as.numeric(ifelse(TST_temp == "n.a", NA, TST_temp))
characteristic_data_frame$TST <- TST
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
sputum_culture <- gsub("tuberculosis", "M.tuberculosis",
                       characteristic_data_frame$sputum_culture)
characteristic_data_frame$sputum_culture <- sputum_culture

bal_culture <- gsub("tuberculosis", "M.tuberculosis",
                    characteristic_data_frame$bal_culture)
characteristic_data_frame$bal_culture <- bal_culture

GFT_GIT <- rep(NA, nrow(characteristic_data_frame))
index_latent <- grep("LTBI", characteristic_data_frame$TBStatus)
GFT_GIT[index_latent] <- "Positive"
characteristic_data_frame$GFT_GIT <- GFT_GIT
characteristic_data_frame$Tissue <- "Whole Blood"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE19435_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)
##### Create Metadata #####
GSE19439_experimentData <- methods::new("MIAME",
                                        name = "Damien Chaussabel",
                                        lab = "Baylor Institute for Immunology Research",
                                        contact = "DChaussabel@benaroyaresearch.org",
                                        title = "An interferon-inducible neutrophil-driven blood transcriptional signature in human tuberculosis",
                                        abstract = "Whole blood collected in EDTA tubes from patients with active TB disease and healthy controls. Blood was then processed or separated sequentially into neutrophil, monocyte, CD4+ or CD8+ populations and then processed. All patients were sampled prior to the initiation of any antimycobacterial therapy.",
                                        url = "10.1038/nature09247",
                                        pubMedIds = "20725040",
                                        other = list(Platform = "Illumina HumanHT-12 V3.0 expression beadchip (GPL6947)"))
GSE19439_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE19439_Non_normalized_data = as.matrix(GSE19439_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE19439_experimentData));GSE19439_sobject
save_raw_files(GSE19439_sobject, path = "data-raw/", geo = geo)
##### Create normalized curated assay #####
GSE19439_normed <- GSE19439_data_list$data_normalized
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE19439_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
