if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE50834"
sequencePlatform <- "GPL10558"
GSE50834_data_list <- readRawData(geo, sequencePlatform)
GSE50834_Non_pvalues <- GSE50834_data_list$data_Non_normalized

##### Match colnames to sample ID #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(seq_len(length(names(GEOquery::GSMList(gse)))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$description[2])
description_id <- gsub("-", ".", description_id_raw)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)
indx <- base::match(ID_table$DescriptionID, colnames(GSE50834_Non_pvalues))
GSE50834_Non_pvalues <- GSE50834_Non_pvalues[,indx]
colnames(GSE50834_Non_pvalues) <- ID_table$SampleID
GSE50834_Non_normalized_data <- GSE50834_Non_pvalues

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Gender", "Tissue", "TBStatus")
characteristic_data_frame$GeographicalRegion <- "South Africa"
characteristic_data_frame$HIVStatus <- "Positive"
characteristic_data_frame$Barcode <- ID_table$DescriptionID
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus <- ifelse(TBStatus_temp == "HIV/TB", "PTB", "Control")
characteristic_data_frame$TBStatus <- TBStatus
Treatment <- rep("3TC+D4T+EFV", nrow(characteristic_data_frame))
Treatment[2:3] <- "3TC+D4T+EFV/3TC+D4T+NVP (Unidentified)"
characteristic_data_frame$Treatment <- Treatment
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE50834_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

###### Create Metadata #####
GSE50834_experimentData <- methods::new("MIAME",
                                        name = "Louise C Showe",
                                        lab = "The Wistar Institute",
                                        contact = "lshowe@wistar.org",
                                        title = "Identification of a 251 gene expression signature that can accurately detect M. tuberculosis in patients with and without HIV co-infection.",
                                        abstract = "Diagnosis of TB, especially in the presence of an HIV co-infection, can be challenging when using conventional diagnostic methods. In this study, we analyzed global gene expression data from PBMC of patients that were either mono-infected with HIV or co-infected with HIV and TB in order to identify a TB-specific gene signature.",
                                        url = "10.1371/journal.pone.0089925",
                                        pubMedIds = "24587128",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE50834_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE50834_Non_normalized_data = as.matrix(GSE50834_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE50834_experimentData));GSE50834_sobject
save_raw_files(GSE50834_sobject, path = "data-raw/", geo = geo)

##### Create normalized curated assay #####
GSE50834_normed <- GSE50834_data_list$data_normalized
colnames(GSE50834_normed) <- names(GEOquery::GSMList(gse))
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE50834_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)


