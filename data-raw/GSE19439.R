if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
urls <- GEOquery::getGEOSuppFiles("GSE19439", fetch_files = FALSE)
url_sub <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_sub, temp)
utils::untar(temp,exdir = tempd)
files <- list.files(tempd, pattern = "txt.*")
GSE19439_Non_normalized_list <- lapply(files, function(x)
  read.delim(paste0(tempd, "/", x), header = TRUE,
             col.names = c("ID_REF", gsub("_.*", "", x), paste0(gsub("_.*", "", x), ".Pval")),
             stringsAsFactors = FALSE))
GSE19439_Non_normalized_list_noPvalue <- lapply(GSE19439_Non_normalized_list, function(x)
  x[, -grep(".Pval", colnames(x))])
# Merge list to a matrix based on the ID_REF
GSE19439_Non_normalized <- Reduce(function(x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE19439_Non_normalized_list_noPvalue, function(x) {x}))
row.names(GSE19439_Non_normalized) <- GSE19439_Non_normalized$ID_REF
GSE19439_Non_normalized_counts <- GSE19439_Non_pvalue <- GSE19439_Non_normalized[-1]

##### Create Column data #####
gse <- GEOquery::getGEO("GSE19439", GSEMatrix = FALSE)
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, "[[", x))
characteristic_data_frame <- sub("(.*?): ", "", characteristic_table) %>%
  S4Vectors::DataFrame()
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
TBStatus <- ifelse(TBStatus == "Latent", "LTBI", TBStatus)

characteristic_data_frame$TBStatus <- TBStatus
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))
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
# Annotation from vendor's information
gpl6947 <- GEOquery::getGEO("GPL6947", GSEMatrix = FALSE)
GPL6947_dat <- gpl6947@dataTable@table %>% data.frame()

# Annotation from Bioconductor package
PROBES <- row.names(GSE19439_Non_pvalue)
OUT <- AnnotationDbi::select(illuminaHumanv3.db::illuminaHumanv3.db, PROBES, "SYMBOL")
OUT[is.na(OUT)] <- NA
# Map ProbeID to Gene Symbol
OUT_collapse <- OUT %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::summarise(SYMBOL = paste(SYMBOL, collapse="///"),
                   times = length(unlist(strsplit(SYMBOL, "///"))))
GSE19439_Non_pvalue$ID_REF <- row.names(GSE19439_Non_pvalue)
GSE19439_final <- GSE19439_Non_pvalue %>%
  dplyr::left_join(OUT_collapse, by=c("ID_REF" = "PROBEID")) %>%
  dplyr::left_join(GPL6947_dat, by = c("ID_REF" = "ID"))
# Create row data
row_data <- GSE19439_final %>%
  dplyr::select(-grep("GSM", colnames(GSE19439_final))) %>%
  S4Vectors::DataFrame()
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
  assays = list(GSE19439_Non_normalized_counts= as.matrix(GSE19439_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE19439_experimentData));GSE19439_sobject
save_raw_files(GSE19439_sobject, path = "data-raw/", geo = "GSE19439")
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
