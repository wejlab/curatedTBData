if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

#### Read in raw data ####
geo <- "GSE152218"
sequencePlatform <- "GPL16791"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
tempd <- tempdir()
utils::download.file(urls$url, temp)
data_counts <- read.delim(temp)
colnames(data_counts)[1] <- "ID_REF"
colnames(data_counts) <- gsub("X", "", colnames(data_counts))
data_counts1 <- data_counts |> 
    tibble::column_to_rownames("ID_REF")
title_to_name <- lapply(1:length(GEOquery::GSMList(gse)), function(x) {
    gsm_info <- GEOquery::GSMList(gse)[[x]]@header
    title1 <- gsm_info$title
    data.frame(title = gsub("-", ".", title1), gsm_name = gsm_info$geo_accession)
}) |> 
    dplyr::bind_rows()
index_col <- match(title_to_name$title, colnames(data_counts1))
data_counts1 <- data_counts1[, index_col]
colnames(data_counts1) <- title_to_name$gsm_name

#### Create Column Data ####
characteristic_data_frame <- readRawColData2(gse)

colnames(characteristic_data_frame) <- c("Age", "TBStatus", "Gender", "BMI")
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$GeographicalRegion <- "India"
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
characteristic_data_frame$HIVStatus <- "Negative"
characteristic_data_frame$BMI <- as.numeric(characteristic_data_frame$BMI)
TBStatus <- ifelse(characteristic_data_frame$TBStatus == "LTBI", "LTBI", "PTB")
characteristic_data_frame$TBStatus <- TBStatus
Gender <- ifelse(characteristic_data_frame$Gender == "male", "Male", "Female")
characteristic_data_frame$Gender <- Gender
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create raw data: convert ensembl to gene symbol #####
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables")
library(annotables)
data_counts1_new <- data_counts1 |> 
    tibble::rownames_to_column("ENSEMBL") |> 
    dplyr::inner_join(grch38, by = c("ENSEMBL" = "ensgene")) |> 
    dplyr::filter(symbol != "") |> 
    dplyr::select(colnames(data_counts1), "symbol")
# Merge duplicated gene names
data_counts1_new_combine <- stats::aggregate(. ~ symbol , data_counts1_new, median)  
data_counts1_new_combine <- data_counts1_new_combine |> 
    tibble::column_to_rownames("symbol")

#### Create Row Data ####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(data_counts1_new_combine),
                                     SYMBOL_NEW = row.names(data_counts1_new_combine))

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "William Evan Johnson",
                      lab = "Boston University",
                      contact = "wej@bu.edu",
                      title = "Tuberculosis in Malnourished Individuals",
                      abstract = "Whole blood gene expression profiling from well and malnourished Indian individuals with TB and severely malnourished household contacts with latent TB infection (LTBI). Severe malnutrition was defined as body mass index (BMI) <16. kg/m2 in adults and based on weight-for-height Z scores in children <18 years. Gene expression was measured using RNA-sequencing.",
                      url = "10.3389/fimmu.2022.1011166",
                      pubMedIds = "36248906",
                      other=list(Platform = "Illumina HiSeq 2500 (Homo sapiens)
"))

sobject <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(data_counts1_new_combine)),
    colData = new_col_info,
    rowData = new_row_data,
    metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(data_counts1_new_combine, paste0("data-raw/", geo, "_assay_curated.RDS"))

