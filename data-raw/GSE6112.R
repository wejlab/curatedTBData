if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
# library(xml2)
# library(tidyverse)
##### Read in Information #####
# Having trouble in matching probeID to gene symbol, using the normalized data for GSE6112
geo <- "GSE6112"
sequencePlatform <- "GPL4475"
# urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
# url_xml <- as.character(urls$url[1])
# temp <- tempfile()
# tempd <- tempdir()
# utils::download.file(url_xml, temp)
# utils::untar(temp, exdir = tempd)
# xmlFiles <- list.files(path = tempd, pattern = '*.xml', full.names=TRUE)
# GSE6112_xml <- lapply(xmlFiles, function(x) xml2::read_xml(x))
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
GSE6112_exprs_list <- lapply(1:length(GEOquery::GSMList(gse)), function(x) {
  x <- GEOquery::GSMList(gse)[[x]]@dataTable@table
  x %>% dplyr::select(ID_REF, Intensity1, Intensity2)
})

GSE6112_normalized_data <- Reduce(function (x, y)
  merge(x, y, by = "ID_REF", all = FALSE),
  lapply(GSE6112_exprs_list, function(x) {x}))
row.names(GSE6112_normalized_data) <- GSE6112_normalized_data$ID_REF
GSE6112_normalized_data <- GSE6112_normalized_data[, -1]
colnames(GSE6112_normalized_data) <- paste0(rep(names(GEOquery::GSMList(gse)), each=2),
                                            rep(c("_ch1","_ch2"), length(names(GEOquery::GSMList(gse)))))
##### Create Column data #####
data_characteristic_1 <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
data_characteristic_2 <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch2)
characteristic_table_1 <- sapply(1:length(data_characteristic_1[[1]]), function(x)
  sapply(data_characteristic_1, "[[", x))
head(characteristic_table_1)

characteristic_table_2 <- sapply(1:length(data_characteristic_2[[1]]), function(x)
  sapply(data_characteristic_2, "[[", x))
characteristic_data_frame <- rbind(characteristic_table_1,characteristic_table_2) %>%
  DataFrame()

row.names(characteristic_data_frame) <- c(paste0(names(GEOquery::GSMList(gse)), "_ch1"),
                                          paste0(names(GEOquery::GSMList(gse)), "_ch2"))
colnames(characteristic_data_frame) <- "TBStatus"
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus[grep('LTBI',TBStatus_temp)] <- "LTBI"
TBStatus[grep('Healthy_contact',TBStatus_temp)] <- "LTBI"
TBStatus[grep('_TB_',TBStatus_temp)] <- "PTB"
TBStatus[grep('Pool',TBStatus_temp)] <- "Control"
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Tissue <- "PBMC"
indx <- match(colnames(GSE6112_normalized_data), row.names(characteristic_data_frame))
characteristic_data_frame <- characteristic_data_frame[indx,]
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)
##### Create Row Data #####
# Annotation for the probe
GPL4475 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
GPL4475_annotate <- GPL4475@dataTable@table
GPL4475_annotate$ID <- as.character(GPL4475_annotate$ID)
GSE6112_probeNames <- data.frame(row.names(GSE6112_normalized_data))
colnames(GSE6112_probeNames) <- c("ID_REF")

new_row_data <- GSE6112_probeNames %>%
  dplyr::left_join(GPL4475_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$Primary.Sequence.Name

##### Create Metadata #####
GSE6112_experimentData <- methods::new("MIAME",
                                      name = "Marc Jacobsen",
                                      lab = "Max-Planck-Institute for Infection Biology",
                                      contact = "jacobsen@mpiib-berlin.mpg.de",
                                      title = "Candidate biomarkers for discrimination between infection and disease caused by Mycobacterium tuberculosis",
                                      abstract = "Infection with Mycobacterium tuberculosis is controlled by an efficacious immune response in about 90% of infected individuals who do not develop disease. Although essential mediators of protection, e.g., interferon-gamma, have been identified, these factors are insufficient to predict the outcome of M. tuberculosis infection. As a first step to determine additional biomarkers, we compared gene expression profiles of peripheral blood mononuclear cells from tuberculosis patients and M. tuberculosis-infected healthy donors by microarray analysis. Differentially expressed candidate genes were predominantly derived from monocytes and comprised molecules involved in the antimicrobial defense, inflammation, chemotaxis, and intracellular trafficking. We verified differential expression for alpha-defensin 1, alpha-defensin 4, lactoferrin, Fcgamma receptor 1A (cluster of differentiation 64 [CD64]), bactericidal permeability-increasing protein, and formyl peptide receptor 1 by quantitative polymerase chain reaction analysis. Moreover, we identified increased protein expression of CD64 on monocytes from tuberculosis patients. Candidate biomarkers were then assessed for optimal study group discrimination. Using a linear discriminant analysis, a minimal group of genes comprising lactoferrin, CD64, and the Ras-associated GTPase 33A was sufficient for classification of (1) tuberculosis patients, (2) M. tuberculosis-infected healthy donors, and (3) noninfected healthy donors.",
                                      url = "10.1007/s00109-007-0157-6",
                                      pubMedIds = "17318616",
                                      other = list(Platform = "MPIIB custom human (AMADID 011412) (GPL4475)"))
GSE6112_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE6112_normalized_data = as.matrix(GSE6112_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE6112_experimentData));GSE6112_sobject
save_raw_files(GSE6112_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE6112_normalized_data,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
