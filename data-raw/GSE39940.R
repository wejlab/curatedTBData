if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}

source("data-raw/UtilityFunctionForCuration.R")
##### Read in Non-normalized data #####
geo <- "GSE39940"
sequencePlatform <- "GPL10558"
GSE39940_Non_pvalue <- readRawData(geo, sequencePlatform)

##### Match colnames to sample ID #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1[4])
# Demo: Convert 'barcode: 6247215037_A' to 6247215037_A
description_id <- sapply(1:length(description_id_raw),
                         function(x) gsub("(.*?): ", "", description_id_raw[x]))
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)
colnames(GSE39940_Non_pvalue) <- sapply(1:ncol(GSE39940_Non_pvalue),
                                        function(x) gsub(".*X|\\..*", "",
                                                         colnames(GSE39940_Non_pvalue)[x]))
indx <- base::match(ID_table$DescriptionID, colnames(GSE39940_Non_pvalue))
GSE39940_Non_pvalue <- GSE39940_Non_pvalue[,indx]
colnames(GSE39940_Non_pvalue) <- ID_table$SampleID
GSE39940_Non_normalized_counts <- GSE39940_Non_pvalue

##### Create Column Data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c('TBStatus','HIVStatus','GeographicalRegion','Barcode')
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$PneumoniaStatus <- "Negative"
characteristic_data_frame$HIVStatus <- ifelse(characteristic_data_frame$HIVStatus == "HIV positive",
                                              "Positive", "Negative")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)

for (i in 1:length(TBStatus)) {
  if (TBStatus[i] == unique(TBStatus_temp)[1] ) {
    TBStatus[i] <- "PTB"
  }
  if (TBStatus[i] == unique(TBStatus_temp)[2] ) {
    TBStatus[i] <- "LTBI"
  }
  if (TBStatus[i] == unique(TBStatus_temp)[3]) {
    TBStatus[i] <- "OD"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)
###### Create Row data #####
row_data <- map_gene_symbol(GSE39940_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

###### Create Metadata #####
# Same with GSE39939
GSE39940_experimentData <- methods::new("MIAME",
                                        name = "Victoria Wright",
                                        lab = "Wright Fleming Institute",
                                        contact = "v.wright@imperial.ac.uk",
                                        title = "Diagnosis of childhood tuberculosis and host RNA expression in Africa",
                                        abstract = "Children were recruited from 2 hospitals in Coast Province, Kenya (n=157) who were either HIV+ or HIV - with either active TB (culture confirmed), active TB (culture negative), LTBI or OD.",
                                        url = "10.1056/NEJMoa1303657",
                                        pubMedIds = "24785206",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE39940_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE39940_Non_normalized_counts= as.matrix(GSE39940_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE39940_experimentData))
save_raw_files(GSE39940_sobject, path = "data-raw/", geo = geo)


