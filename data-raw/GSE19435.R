if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE19435"
sequencePlatform <- "GPL6947"
GSE19435_Non_normalized_list_noPvalue <- readRawData(geo, sequencePlatform)
# Merge list to a matrix based on the ID_REF
GSE19435_Non_normalized <- Reduce(function (x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE19435_Non_normalized_list_noPvalue, function(x) {x}))
row.names(GSE19435_Non_normalized) <- GSE19435_Non_normalized$ID_REF
GSE19435_Non_normalized_data <- GSE19435_Non_pvalue <- GSE19435_Non_normalized[-1]

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)

colnames(characteristic_data_frame) <- c("Age", "Gender", "Ethnicity", "TBStatus",
                                         "MeasurementTime", "SputumSmearStatus",
                                         "isolate_sensitivity", "BcgVaccinated")
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
title <- sapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title)
title_split <- strsplit(title,"_")
characteristic_data_frame$PatientID <- sapply(title_split, function(x) x[3])
Treatment <- ifelse(characteristic_data_frame$TBStatus == "PTB", "anti-mycobacterial treatment", NA)
characteristic_data_frame$Treatment <- Treatment
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE19435_Non_pvalue, sequencePlatform)
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
  assays = list(GSE19435_Non_normalized_data = as.matrix(GSE19435_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE19435_experimentData));GSE19435_sobject
save_raw_files(GSE19435_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
