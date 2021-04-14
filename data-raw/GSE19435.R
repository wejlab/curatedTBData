if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
urls <- GEOquery::getGEOSuppFiles("GSE19435", fetch_files = FALSE)
url_sub <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_sub, temp)
utils::untar(temp,exdir = tempd)
files <- list.files(tempd, pattern = "txt.*")
GSE19435_Non_normalized_list <- lapply(files, function(x)
  read.delim(paste0(tempd,"/",x), header = TRUE,
             col.names = c("ID_REF", gsub("_.*", "", x), paste0(gsub("_.*", "", x), ".Pval")),
             stringsAsFactors = FALSE))
GSE19435_Non_normalized_list_noPvalue <- lapply(GSE19435_Non_normalized_list, function(x)
  x[, -grep('.Pval', colnames(x))])
# Merge list to a matrix based on the ID_REF
GSE19435_Non_normalized <- Reduce(function(x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE19435_Non_normalized_list_noPvalue, function(x) {x}))
row.names(GSE19435_Non_normalized) <- GSE19435_Non_normalized$ID_REF
GSE19435_Non_normalized_counts <- GSE19435_Non_pvalue <- GSE19435_Non_normalized[-1]

##### Create Column data #####
gse <- GEOquery::getGEO("GSE19435", GSEMatrix = FALSE)
data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
  sapply(data_characteristic, "[[", x))
characteristic_data_frame <- sub("(.*?): ", "", characteristic_table) %>%
  S4Vectors::DataFrame()
colnames(characteristic_data_frame) <- c("Age", "Gender", "Ethnicity", "TBStatus",
                                         "MeasurementTime", "SputumSmearStatus",
                                         "isolate_sensitivity", "BcgVaccinated")
row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))
characteristic_data_frame$GeographicalRegion <- "London"
characteristic_data_frame$Age <- as.numeric(gsub("years.*", "", characteristic_data_frame$Age))
SputumSmearStatus <- characteristic_data_frame$SputumSmearStatus
SputumSmearStatus[grep("n.a", SputumSmearStatus)] <- NA
characteristic_data_frame$SputumSmearStatus <- SputumSmearStatus
isolate_sensitivity <- characteristic_data_frame$isolate_sensitivity
isolate_sensitivity[grep("n.a", isolate_sensitivity)] <- NA
characteristic_data_frame$isolate_sensitivity <- isolate_sensitivity
characteristic_data_frame$DiabetesStatus <- "Negative"
characteristic_data_frame$HIVStatus <- "Negative"
title <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title) %>% unlist
title_split <- strsplit(title,"_")
PatientID <- sapply(title_split, function(x) x[3])
characteristic_data_frame$PatientID <- PatientID
Treatment <- ifelse(characteristic_data_frame$TBStatus == "PTB", "anti-mycobacterial treatment", NA)
characteristic_data_frame$Treatment <- Treatment
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
# Annotation from vendor's information
gpl6947 <- GEOquery::getGEO("GPL6947", GSEMatrix = FALSE)
GPL6947_dat <- gpl6947@dataTable@table %>% data.frame()

# Annotation from Bioconductor package
PROBES <- row.names(GSE19435_Non_pvalue)
OUT <- AnnotationDbi::select(illuminaHumanv3.db::illuminaHumanv3.db, PROBES, "SYMBOL")
OUT[is.na(OUT)] <- NA
# Map ProbeID to Gene Symbol
OUT_collapse <- OUT %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::summarise(SYMBOL = paste(SYMBOL, collapse="///"),
                   times = length(unlist(strsplit(SYMBOL, "///"))))
GSE19435_Non_pvalue$ID_REF <- row.names(GSE19435_Non_pvalue)
GSE19435_final <- GSE19435_Non_pvalue %>%
  dplyr::left_join(OUT_collapse, by=c("ID_REF" = "PROBEID")) %>%
  dplyr::left_join(GPL6947_dat, by = c("ID_REF" = "ID"))
# Create row data
row_data <- GSE19435_final %>%
  dplyr::select(-grep("GSM", colnames(GSE19435_final))) %>%
  S4Vectors::DataFrame()
new_row_data <- match_gene_symbol(row_data)
##### Create Metadata #####
GSE19435_experimentData <- methods::new("MIAME",
                                        name = "Damien Chaussabel",
                                        lab = "Baylor Institute for Immunology Research",
                                        contact = "DChaussabel@benaroyaresearch.org",
                                        title = "An interferon-inducible neutrophil-driven blood transcriptional signature in human tuberculosis",
                                        abstract = "Whole blood collected in EDTA tubes from patients with active TB disease and healthy controls. Blood was then processed or separated sequentially into neutrophil, monocyte, CD4+ or CD8+ populations and then processed. All patients were sampled prior to the initiation of any antimycobacterial therapy.",
                                        url = "10.1038/nature09247",
                                        pubMedIds = "20725040",
                                        other = list(Platform = "Illumina HumanHT-12 V3.0 expression beadchip (GPL6947)"))
GSE19435_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE19435_Non_normalized_counts= as.matrix(GSE19435_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE19435_experimentData));GSE19435_sobject
save_raw_files(GSE19435_sobject, path = "data-raw/", geo = "GSE19435")
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
