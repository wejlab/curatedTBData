if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE37250"
sequencePlatform <- "GPL10558"
GSE37250_Non_pvalue <- readRawData(geo, sequencePlatform)
# Remove SYMBOL SEARCH_KEY ILMN_GENE CHROMOSOME DEFINITION SYNONYMS
GSE37250_Non_pvalue <- GSE37250_Non_pvalue[, grep("X.*", colnames(GSE37250_Non_pvalue))]

##### Match colnames to sample ID #####
# Obtain raw data information from GEO
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@dataTable@columns$Column[3])
# Convert 5667664060_K.Detection Pval to 5667664060_K
description_id <- gsub("\\..*", "", description_id_raw)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)), DescriptionID = description_id)

# Demo: convert X6247215025_L.AVG_Signal to 6247215025_L, should use \\. when mathcing .
colnames(GSE37250_Non_pvalue) <- gsub(".*X|\\..*", "", colnames(GSE37250_Non_pvalue))
indx <- base::match(ID_table$DescriptionID, colnames(GSE37250_Non_pvalue))
GSE37250_Non_pvalue <- GSE37250_Non_pvalue[,indx]
colnames(GSE37250_Non_pvalue) <- ID_table$SampleID
GSE39939_Non_normalized_counts <- GSE37250_Non_pvalue

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
  if(TBStatus_temp[i]=="latent TB infection"){
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


