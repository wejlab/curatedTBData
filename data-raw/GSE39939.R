if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE39939"
sequencePlatform <- "GPL10558"
GSE39939_data_list <- readRawData(geo, sequencePlatform)
GSE39939_Non_pvalue <- GSE39939_data_list$data_Non_normalized

##### Match colnames to sample ID #####
# Obtain raw data information from GEO
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)

description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@dataTable@columns$Column[3])
# Example: Convert 6116725094_J.Detection to 6116725094_J
description_id <- gsub("\\..*", "", description_id_raw)

ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)

# Mapping descriptionID to SampleID
colnames(GSE39939_Non_pvalue) <- gsub(".*X|\\..*", "", colnames(GSE39939_Non_pvalue))
indx <- base::match(ID_table$DescriptionID, colnames(GSE39939_Non_pvalue))
GSE39939_Non_pvalue <- GSE39939_Non_pvalue[,indx]
colnames(GSE39939_Non_pvalue) <- ID_table$SampleID
GSE39939_Non_normalized_data <- GSE39939_Non_pvalue

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus","HIVStatus","GeographicalRegion")
characteristic_data_frame$Barcode <- ID_table$DescriptionID
characteristic_data_frame$PneumoniaStatus <- "Negative"
characteristic_data_frame$HIVStatus <- ifelse(characteristic_data_frame$HIVStatus == "HIV positive",
                                              "Positive", "Negative")
characteristic_data_frame$Tissue <- "Whole Blood"

col_info <- create_standard_coldata(characteristic_data_frame)

relabel_TB <- function(dat_new) {
  dat_new$CultureStatus <- sub(".*\\((.*)\\).*", "\\1", dat_new$TBStatus)

  dat_new$CultureStatus <- ifelse(dat_new$CultureStatus==dat_new$TBStatus,
                                  NA, dat_new$CultureStatus)

  TBStatus <- TBStatus_temp <- dat_new$TBStatus

  for (i in 1:length(TBStatus)){
    if (TBStatus[i] == unique(TBStatus_temp)[1] || TBStatus[i] == unique(TBStatus_temp)[4]){
      TBStatus[i] <- "PTB"
    }
    if (TBStatus[i] == unique(TBStatus_temp)[2] || TBStatus[i] == unique(TBStatus_temp)[5]){
      TBStatus[i] <- "OD"
    }
    if (TBStatus[i] == unique(TBStatus_temp)[3]){
      TBStatus[i] <- "LTBI"
    }
  }
  dat_new$TBStatus <- TBStatus
  return(S4Vecotrs::DataFrame(dat_new))
}

new_col_info <- relabel_TB(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE39939_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)
###### Create Metadata #####
GSE39939_experimentData <- methods::new("MIAME",
                               name = "Victoria Wright",
                               lab = "Wright Fleming Institute",
                               contact = "v.wright@imperial.ac.uk",
                               title = "Diagnosis of childhood tuberculosis and host RNA expression in Africa",
                               abstract = "Children were recruited from 2 hospitals in Coast Province, Kenya (n=157) who were either HIV+ or HIV - with either active TB (culture confirmed), active TB (culture negative), LTBI or OD.",
                               url = "10.1056/NEJMoa1303657",
                               pubMedIds = "24785206",
                               other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE39939_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE39939_Non_normalized_data= as.matrix(GSE39939_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE39939_experimentData))
save_raw_files(GSE39939_sobject, path = "data-raw/", geo = geo)

##### Create normalized curated assay #####
GSE39939_normed <- GSE39939_data_list$data_normalized[,indx]
colnames(GSE39939_normed) <- ID_table$SampleID
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE39939_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))

# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
