if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}

source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE69581"
sequencePlatform <- "GPL10558"
GSE69581_data_list <- readRawData(geo, sequencePlatform)
GSE69581_Non_pvalue <- GSE69581_data_list$data_Non_normalized

##### Match colnames to sample ID #####
# Obtain raw data information from GEO
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)

description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@dataTable@columns$Column[3])
# Example: convert PValue_9421742028_J to 9421742028_J
description_id <- sapply(1:length(description_id_raw),
                         function(x) gsub("^(.*?)_", "", description_id_raw[x]))
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)

# Mapping descriptionID to SampleID
# Example: onvert X9421742043_G to 9421742043_G
colnames(GSE69581_Non_pvalue) <- gsub("^(.*?)X", "", colnames(GSE69581_Non_pvalue))
indx <- base::match(ID_table$DescriptionID, colnames(GSE69581_Non_pvalue))
GSE69581_Non_pvalue <- GSE69581_Non_pvalue[,indx]
colnames(GSE69581_Non_pvalue) <- ID_table$SampleID
GSE69581_Non_normalized_data <- GSE69581_Non_pvalue

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Tissue", "HIVStatus", "Treatment",
                                         "TBStatus", "MeasurementTime", "PatientID")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
for(i in 1:length(TBStatus)){
  if(TBStatus_temp[i] == "Active"){
    TBStatus[i] <- "PTB"
  } else if (TBStatus_temp[i] == "Latent") {
    TBStatus[i] <- "LTBI"
  }
}
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Barcode <- ID_table$DescriptionID
characteristic_data_frame$GeographicalRegion <- "South Africa"
characteristic_data_frame$HIVStatus <- "Positive"
characteristic_data_frame$Treatment <- "Untreated"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE69581_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

###### Create Metadata #####
GSE69581_experimentData <- new("MIAME",
                               name = "Hanif Esmail",
                               lab = "University of Cape Town",
                               contact = "hanif.esmail@ndcls.ox.ac.uk",
                               title = "Complement Pathway Gene Activation and Rising Circulating Immune Complexes Characterize Early Disease in HIV-Associated Tuberculosis",
                               abstract = "Total RNA extracted from whole blood of 35 HIV infected participants with evidence of latent TB (10 with subclinical pathology on PET/CT and 25 without) and 15 HIV infected participants with active TB",
                               url = "10.1073/pnas.1711853115",
                               pubMedIds = "29339504",
                               other=list(platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))

GSE69581_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE69581_Non_normalized_data = as.matrix(GSE69581_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE69581_experimentData));GSE69581_sobject
save_raw_files(GSE69581_sobject, path = "data-raw/", geo = geo)

##### Create normalized curated assay #####
GSE69581_normed <- GSE69581_data_list$data_normalized[,indx]
colnames(GSE69581_normed) <- ID_table$SampleID
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE69581_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
