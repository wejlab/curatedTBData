if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}

source("data-raw/UtilityFunctionForCuration.R")
##### Read in Non-normalized data #####
urls <- GEOquery::getGEOSuppFiles("GSE39940", fetch_files = FALSE)
url_non_normalized <- as.character(urls$url[2])

temp <- tempfile()
utils::download.file(url_non_normalized, temp)

GSE39940_Non_normalized <- read.delim(gzfile(temp), row.names = 1, header = TRUE)
# Remove P.val
GSE39940_Non_pvalue <- GSE39940_Non_normalized[,-grep('.Pval',
                                                      colnames(GSE39940_Non_normalized))]

##### Match colnames to sample ID #####
gse <- GEOquery::getGEO("GSE39940", GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1[4])
# Demo: Convert 'barcode: 6247215037_A' to 6247215037_A
description_id <- sapply(1:length(description_id_raw),
                         function(x) gsub('(.*?): ', "", description_id_raw[x]))
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
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, '[[',x))
## Demo convert 'geographical region: Kenya' to 'Kenya'
characteristic_data_frame <- gsub('(.*?): ','',characteristic_table) %>%
  S4Vectors::DataFrame()
colnames(characteristic_data_frame) <- c('TBStatus','HIVStatus','GeographicalRegion','Barcode')
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$PneumoniaStatus <- "Negative"
characteristic_data_frame$HIVStatus <- ifelse(characteristic_data_frame$HIVStatus == "HIV positive",
                                              "Positive", "Negative")
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))

new_col_info <- S4Vectors::DataFrame(create_standard_coldata(characteristic_data_frame))
TBStatus <- TBStatus_temp <- as.character(new_col_info$TBStatus)
unique(TBStatus_temp)

for (i in 1:length(TBStatus)){
  if (TBStatus[i] == unique(TBStatus_temp)[1] ){
    TBStatus[i] = 'PTB'
  }
  if (TBStatus[i] == unique(TBStatus_temp)[2] ){
    TBStatus[i] = 'LTBI'
  }
  if (TBStatus[i] == unique(TBStatus_temp)[3]){
    TBStatus[i] = 'OD'
  }
}

new_col_info$TBStatus <- TBStatus

###### Create Row data #####
GPL10558_HumanHT <- GEOquery::getGEO("GPL10558", GSEMatrix = F)
GPL10558_HumanHT_dat <- GPL10558_HumanHT@dataTable@table

# Annotation
PROBES <- row.names(GSE39940_Non_pvalue)
OUT <- AnnotationDbi::select(illuminaHumanv4.db::illuminaHumanv4.db, PROBES, "SYMBOL")
OUT[is.na(OUT)] <- NA
# Map ProbeID to Gene Symbol
OUT_collapse <- OUT %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::summarise(SYMBOL = paste(SYMBOL, collapse="///"),
                   times = length(unlist(strsplit(SYMBOL, "///"))))

GSE39940_Non_pvalue$ID_REF <- row.names(GSE39940_Non_pvalue)
GSE39940_final <- GSE39940_Non_pvalue %>%
  dplyr::left_join(OUT_collapse, by=c("ID_REF" = "PROBEID")) %>%
  dplyr::left_join(GPL10558_HumanHT_dat, by = c("ID_REF" = "ID"))
# Create Row data annotation
row_data <- GSE39940_final %>%
  dplyr::select(-grep("GSM",colnames(GSE39940_final))) %>%
  S4Vectors::DataFrame()
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
save_raw_files(GSE39940_sobject, path = "data-raw/", geo = "GSE39940")


