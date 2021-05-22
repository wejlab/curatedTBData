if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE28623"
sequencePlatform <- "GPL4133"
temp <- tempfile()
tempd <- tempdir()
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
# Large data file 843.6MB
utils::download.file(as.character(urls$url[1]), temp)
utils::untar(temp, exdir = tempd)
filesPath <- list.files(tempd, pattern = "GSM.*", full.names = TRUE)
GSE28623_raw <- limma::read.maimages(filesPath, source="agilent", green.only = TRUE)
GSE28623_Non_normalized_data <- GSE28623_raw$E
row.names(GSE28623_Non_normalized_data) <- GSE28623_raw$genes$ProbeName
col_name1 <- gsub(".*/", "", colnames(GSE28623_Non_normalized_data))
colnames(GSE28623_Non_normalized_data) <- gsub("_.*", "", col_name1)
GSE28623_Non_pvalue <- GSE28623_Non_normalized_data

##### Create curated assay ####
curatedExprs <- norm_probeToGenes_Agilent(GSE28623_raw, FUN = median)
colnames(curatedExprs) <- colnames(GSE28623_Non_normalized_data)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))


##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "Gender")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
for(i in 1:length(TBStatus)){
  if(TBStatus_temp[i] == "healthy non-infected donors"){
    TBStatus[i] = "Control"
  }
  if(TBStatus_temp[i] == "healthy donors latently infected with M. tuberculosis"){
    TBStatus[i] = "LTBI"
  }
  if(TBStatus_temp[i] == "tuberculosis patients"){
    TBStatus[i] = 'PTB'
  }
}
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$GeographicalRegion <- "The Gambia"
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
gpl4133 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
gpl4133_annotate <- gpl4133@dataTable@table %>%
  dplyr::distinct(NAME, .keep_all = TRUE)

GSE28623_probeNames <- data.frame(row.names(GSE28623_Non_pvalue))
colnames(GSE28623_probeNames) <- c("ID_REF")

new_row_data <- GSE28623_probeNames %>%
  dplyr::left_join(gpl4133_annotate, by = c("ID_REF" = "SPOT_ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$GENE_SYMBOL

##### Create Metadata #####
GSE28623_experimentData <- methods::new("MIAME",
                                        name = "Jeroen Maertzdorf",
                                        lab = "MPIIB",
                                        contact = "maertzdorf@mpiib-berlin.mpg.de",
                                        title = "Functional correlations of pathogenesis-driven gene expression signatures in tuberculosis.",
                                        abstract = "A total of 46 sputum smear and chest x-ray positive TB patients (TB), 25 latently infected healthy donors (LTBI, TB case contacts with Mantoux test induration >= 10mm) and 37 uninfected donors (NID, Mantoux test induration = 0mm) were included from a cohort in The Gambia.",
                                        url = "10.1371/journal.pone.0026938",
                                        pubMedIds = "22046420",
                                        other = list(Platform = "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Feature Number version) (GPL4133)"))
GSE28623_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE28623_Non_normalized_data = as.matrix(GSE28623_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE28623_experimentData));GSE28623_sobject
save_raw_files(GSE28623_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)



