if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in Non-normalized data #####
urls <- GEOquery::getGEOSuppFiles("GSE39939", fetch_files = FALSE)
url_non_normalized <- as.character(urls$url[2])
temp <- tempfile()
utils::download.file(url_non_normalized,temp)
GSE39939_Non_normalized <- read.delim(gzfile(temp), row.names = 1, header = TRUE)

# Remove .Pval from the column
GSE39939_Non_pvalue <- GSE39939_Non_normalized[,-grep('.Pval', colnames(GSE39939_Non_normalized))]
dim(GSE39939_Non_pvalue)

##### Process Column names of the raw data. Map them to sample name #####
# Obtain raw data information from GEO
gse <- GEOquery::getGEO("GSE39939", GSEMatrix = FALSE)

description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@dataTable@columns$Column[3])
# Example: Convert 6116725094_J.Detection Pval to 6116725094_J
description_id <- sapply(1:length(description_id_raw),
                         function(x) gsub("\\..*", "", description_id_raw[x]))
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)), DescriptionID = description_id)

# Mapping descriptionID to SampleID
colnames(GSE39939_Non_pvalue) <- sapply(1:ncol(GSE39939_Non_pvalue),
                                        function(x) gsub(".*X|\\..*", "", colnames(GSE39939_Non_pvalue)[x]))
indx <- base::match(ID_table$DescriptionID, colnames(GSE39939_Non_pvalue))
GSE39939_Non_pvalue <- GSE39939_Non_pvalue[,indx]
colnames(GSE39939_Non_pvalue) <- ID_table$SampleID
GSE39939_Non_normalized_counts <- GSE39939_Non_pvalue

##### Create column data #####
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, '[[',x))
## Demo convert geographical region: Kenya to Kenya
characteristic_data_frame <- sub('(.*?): ','',characteristic_table) %>% as_tibble()
colnames(characteristic_data_frame) <- c("TBStatus","HIVStatus","GeographicalRegion")

characteristic_data_frame <- dplyr::mutate(characteristic_data_frame,
                                           Barcode=ID_table$DescriptionID) %>% DataFrame()
characteristic_data_frame$PneumoniaStatus <- "Negative"
characteristic_data_frame$HIVStatus <- ifelse(characteristic_data_frame$HIVStatus == "HIV positive",
                                              "Positive", "Negative")
characteristic_data_frame$Tissue <- "Whole Blood"
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))

col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- relabel_TB(col_info)

###### Create Row data #####
GPL10558_HumanHT <- GEOquery::getGEO("GPL10558", GSEMatrix = F)
GPL10558_HumanHT_dat <- GPL10558_HumanHT@dataTable@table

# Annotation
PROBES <- row.names(GSE39939_Non_pvalue)
OUT <- AnnotationDbi::select(illuminaHumanv4.db::illuminaHumanv4.db, PROBES, "SYMBOL")
OUT[is.na(OUT)] <- NA

# Map ProbeID to Gene Symbol
OUT_collapse <- OUT %>% dplyr::group_by(PROBEID) %>%
  dplyr::summarise(SYMBOL = paste(SYMBOL, collapse = "///"),
                   times = length(unlist(strsplit(SYMBOL, "///"))))

GSE39939_Non_pvalue$ID_REF <- row.names(GSE39939_Non_pvalue)
GSE39939_final <- GSE39939_Non_pvalue %>%
                  dplyr::eft_join(OUT_collapse, by = c("ID_REF" = "PROBEID")) %>%
                  dplyr::left_join(GPL10558_HumanHT_dat, by = c("ID_REF" = "ID"))

# Modify Row data annotation
row_data <- GSE39939_final %>%
  dplyr::select(-grep("GSM",colnames(GSE39939_final))) %>% DataFrame()
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
  assays = list(GSE39939_Non_normalized_counts= as.matrix(GSE39939_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE39939_experimentData))
save_raw_files(GSE39939_sobject, path = "data-raw/", geo = "GSE39939")

# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
