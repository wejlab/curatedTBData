if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE37250"
sequencePlatform <- "GPL10558"
GSE37250_Non_pvalues <- readRawData(geo, sequencePlatform)
# Remove SYMBOL SEARCH_KEY ILMN_GENE CHROMOSOME DEFINITION SYNONYMS
GSE37250_Non_pvalues <- GSE37250_Non_pvalues[, grep("X.*", colnames(GSE37250_Non_pvalues))]

##### Match colnames to sample ID #####
# Obtain raw data information from GEO
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@dataTable@columns$Column[3])
# Convert 5667664060_K.Detection Pval to 5667664060_K
description_id <- gsub("\\..*", "", description_id_raw)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)), DescriptionID = description_id)

# Demo: convert X6247215025_L.AVG_Signal to 6247215025_L, should use \\. when mathcing .
colnames(GSE37250_Non_pvalues) <- gsub(".*X|\\..*", "", colnames(GSE37250_Non_pvalues))
indx <- base::match(ID_table$DescriptionID, colnames(GSE37250_Non_pvalues))
GSE37250_Non_pvalues <- GSE37250_Non_pvalues[,indx]
colnames(GSE37250_Non_pvalues) <- ID_table$SampleID
GSE37250_Non_normalized_counts <- GSE37250_Non_pvalues

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "HIVStatus", "GeographicalRegion",
                                         "Tissue")
characteristic_data_frame$Barcode <- ID_table$DescriptionID
characteristic_data_frame$Tissue <- "Whole Blood"
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
for(i in 1:length(TBStatus)){
  if(TBStatus_temp[i] == "active tuberculosis"){
    TBStatus[i] <- "PTB"
  }
  if(TBStatus_temp[i] == "latent TB infection"){
    TBStatus[i] <- "LTBI"
  }
  if(TBStatus_temp[i] == "other disease"){
    TBStatus[i] <- "OD"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
HIVStatus <- characteristic_data_frame$HIVStatus
characteristic_data_frame$HIVStatus <- ifelse(HIVStatus == "HIV negative",
                                              "Negative", "Positive")
index_PTB <- grep("PTB",characteristic_data_frame$TBStatus)
characteristic_data_frame$sputum_culture <- NA
characteristic_data_frame$sputum_culture[index_PTB] <- "M.tuberculosis"
# Add TST information for LTBI
TST <- rep(NA, nrow(characteristic_data_frame))
index_latent_positive <- Reduce(intersect, list(which(characteristic_data_frame$TBStatus == "LTBI"),
                                                which(characteristic_data_frame$HIVStatus == "Positive")))
index_latent_negative <- Reduce(intersect, list(which(characteristic_data_frame$TBStatus == "LTBI"),
                                                which(characteristic_data_frame$HIVStatus == "Negative")))
TST[index_latent_positive] <- ">= 5mm"
TST[index_latent_negative] <- ">= 10mm"
characteristic_data_frame$TST <- TST

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE37250_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

###### Create Metadata #####
GSE37250_experimentData <- methods::new("MIAME",
                                        name = "Victoria Wright",
                                        lab = "Wright Fleming Institute",
                                        contact = "v.wright@imperial.ac.uk",
                                        title = "Detection of tuberculosis in HIV-infected and -uninfected African adults using whole blood RNA expression signatures: a case-control study.",
                                        abstract = "Adults were recruited from Cape Town, South Africa (n=300) and Karonga, Malawi (n=237) who were either HIV+ or HIV - with either active TB, LTBI or OD.",
                                        url = "10.1371/journal.pmed.1001538",
                                        pubMedIds = "24167453",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE37250_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE37250_Non_normalized_counts= as.matrix(GSE37250_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE37250_experimentData))
save_raw_files(GSE37250_sobject, path = "data-raw/", geo = geo)

# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
