if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

geo <- "GSE41055"
sequencePlatform <- "GPL5175"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_cel <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
utils::download.file(url_cel,temp)
utils::untar(temp,exdir = tempd)

celFiles <- list.files(path = tempd, pattern = "*.CEL", full.names=TRUE)
GSE41055_data.oligo <- oligo::read.celfiles(celFiles)

GSE41055_normalized_rma <- Biobase::exprs(oligo::rma(GSE41055_data.oligo)) # A matrix
colnames(GSE41055_normalized_rma) <- gsub('_.*','',colnames(GSE41055_normalized_rma))
GSE41055_normalized_data <- GSE41055_normalized_rma

##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("Gender", "Age", "sputum_culture",
                                         "TST", "QFT_GIT", "ChestRadiographs")
characteristic_data_frame$GeographicalRegion <- "Venezuela"
characteristic_data_frame$Tissue <- "Whole Blood"
TBStatus <- TBStatus_temp <- sapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$title) %>% as.character()
unique(TBStatus_temp)
TBStatus[grep("healthy control",TBStatus_temp)] <- "Control"
TBStatus[grep("latent",TBStatus_temp)] <- "LTBI"
TBStatus[grep("active",TBStatus_temp)] <- "PTB"
characteristic_data_frame$TBStatus <- TBStatus
sputum_culture <- ifelse(characteristic_data_frame$sputum_culture == "not performed",
                         NA, characteristic_data_frame$sputum_culture)
characteristic_data_frame$sputum_culture <- sputum_culture
characteristic_data_frame$Age <- as.numeric(characteristic_data_frame$Age)
characteristic_data_frame$TST <- as.numeric(characteristic_data_frame$TST)
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
GPL5175 <- GEOquery::getGEO("GPL5175",GSEMatrix = FALSE)
GPL5175_annotate <- GPL5175@dataTable@table
GPL5175_annotate$ID <- as.character(GPL5175_annotate$ID)
GSE41055_probeNames <- data.frame(row.names(GSE41055_normalized_data))
colnames(GSE41055_probeNames) <- c("ID_REF")

new_row_data <- GSE41055_probeNames %>%
  dplyr::left_join(GPL5175_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
SYMBOL_ori <- as.character(new_row_data$gene_assignment)
SYMBOL_ori <- ifelse(SYMBOL_ori == "---", "NA", SYMBOL_ori)
index <- which(SYMBOL_ori != "NA")
SYMBOL <- strsplit(SYMBOL_ori[index], "///")

SYMBOL_new <- lapply(1:length(SYMBOL),
                     function(i) lapply(strsplit(SYMBOL[[i]], "//"), "[[", 2) %>%
                       unlist() %>% unique())
SYMBOL_collapse <- sapply(SYMBOL_new, function(x) gsub(" ", "", paste(x, collapse = "///")))
SYMBOL_ori[index] <- SYMBOL_collapse

new_row_data$SYMBOL_NEW <- SYMBOL_ori

##### Create Metadata #####
GSE41055_experimentData <- methods::new("MIAME",
                                        name = "Aldert Zomer",
                                        lab = "Utrecht University",
                                        contact = "A.L.Zomer@uu.nl",
                                        title = "A predictive signature gene set for discriminating active from latent tuberculosis in Warao Amerindian children.",
                                        abstract = "27 samples in total where analyzed, active TB infection (TB, n=9), Latent TB infection (LTBI, n=9) and healthy controls (HC, n=9) . Gene expression values were log2-transformed and differentially expressed genes were identified based on log2 fold changes (M-values).",
                                        url = "10.1186/1471-2164-14-74",
                                        pubMedIds = "23375113",
                                        other = list(Platform = "Affymetrix Human Exon 1.0 ST Array (GPL5175)"))
GSE41055_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE41055_Normalized_data = as.matrix(GSE41055_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE41055_experimentData));GSE41055_sobject
save_raw_files(GSE41055_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE41055_normalized_data,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))

