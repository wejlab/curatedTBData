if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE62525"
sequencePlatform <- "GPL16951"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_cel <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_cel,temp)
utils::untar(temp,exdir = tempd)

gprFiles <- list.files(path = tempd, pattern = '*.gpr', full.names = TRUE)
# Label cy5
GSE62525_RG <- limma::read.maimages(gprFiles, source = "genepix",
                                    columns = list(R = "F635 Median", G = "B635 Median"))
GSE62525_Non_normalized <- log(GSE62525_RG$R/GSE62525_RG$G, 2)
col_name1 <- gsub(".*/", "", colnames(GSE62525_RG$R))
colnames(GSE62525_Non_normalized) <- gsub("_.*", "", col_name1)
row.names(GSE62525_Non_normalized) <- GSE62525_RG$genes$ID

##### Create curated assay ####
GSE62525_RG_background <- limma::backgroundCorrect(GSE62525_RG, method = "normexp")
GSE62525_RG_normalized <- limma::normalizeBetweenArrays(GSE62525_RG_background,
                                                        method = "quantile")
NoSymbol <- is.na(GSE62525_RG_normalized$genes$Name)
GSE62525_RG_normalized_filt <- GSE62525_RG_normalized[!NoSymbol, ]
GSE62525_normalized_data <- data.frame(GSE62525_RG_normalized_filt$M)
colnames(GSE62525_normalized_data) <- gsub("_.*", "", col_name1)
GSE62525_normalized_data$SYMBOL <- GSE62525_RG_normalized_filt$genes$Name
curatedExprs <- stats::aggregate(. ~ SYMBOL, data = GSE62525_normalized_data,
                           FUN = median, na.action = na.pass)
row.names(curatedExprs) <- curatedExprs$SYMBOL
curatedExprs <- curatedExprs[,-1]
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus")
characteristic_data_frame$Tissue <- "PBMC"
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)

TBStatus[grep("healthy",TBStatus_temp)] <- "Control"
TBStatus[grep("latent",TBStatus_temp)] <- "LTBI"
TBStatus[grep("active TB",TBStatus_temp)] <- "PTB"
characteristic_data_frame$TBStatus <- TBStatus
Title <- data_characteristic <- sapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title)
title_split <- strsplit(Title, "_")
PatientID <- unlist(lapply(title_split, function(x) x[3]))
characteristic_data_frame$PatientID <- PatientID
characteristic_data_frame$GeographicalRegion <- "Taiwan"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create row data ####
gpl16951 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
gpl16951_annotate <- gpl16951@dataTable@table
gpl16951_annotate$ID <- as.character(gpl16951_annotate$ID)
GSE62525_probeNames <- data.frame(row.names(GSE62525_Non_normalized))
colnames(GSE62525_probeNames) <- c("ID_REF")

new_row_data <- GSE62525_probeNames %>%
  dplyr::left_join(gpl16951_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$Gene_symbol
##### Create Metadata #####
GSE62525_experimentData <- methods::new("MIAME",
                                        name = "Julia Tzu-Ya Weng",
                                        lab = "Yuan Ze University",
                                        contact = "julweng@gmail.com",
                                        title = "Gene expression profiling identifies candidate biomarkers for active and latent tuberculosis",
                                        abstract = "Tuberculosis (TB) is a serious infectious disease, but current methods of detection require improvement in sensitivity, efficiency or specificity.
We conducted a microarray experiment, comparing the gene expression profiles in peripheral blood mononuclear cells among individuals with active TB, latent infection, and healthy conditions in a Taiwanese population. These differentially expressed genes may be potential biomarkers that can differentiate between active TB and latent infection.",
                                        url = "10.1186/s12859-015-0848-x",
                                        pubMedIds = "26818387",
                                        other = list(Platform = "Phalanx Human OneArray Ver. 6 Release 1 (GPL1695)"))
GSE62525_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE62525_Non_normalized_data = as.matrix(GSE62525_Non_normalized)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE62525_Non_normalized));GSE62525_sobject
save_raw_files(GSE62525_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)





