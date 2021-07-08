if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE25534"
sequencePlatform <- "GPL1708"
urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
url_raw <- as.character(urls$url[1])
temp <- tempfile()
tempd <- tempdir()
# Large file 644.7 MB
utils::download.file(url_raw, temp)
utils::untar(temp, exdir = tempd)
files <- list.files(tempd, pattern = "GSM.*", full.names = TRUE)
GSE25534_raw <- limma::read.maimages(files, source = "agilent")
GSE25534_backgroundCorrect <- limma::backgroundCorrect(GSE25534_raw, method = "normexp")
# Cy3 is G (green) Cy5 is R (red)
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
sample_name <- names(GEOquery::GSMList(gse))
sample_name_G <- paste0(sample_name, "_ch1")
sample_name_R <- paste0(sample_name, "_ch2")
GSE25534_Non_normalized_data <- cbind(GSE25534_backgroundCorrect$G, GSE25534_backgroundCorrect$R)
row.names(GSE25534_Non_normalized_data) <- GSE25534_raw$genes$ProbeUID
colnames(GSE25534_Non_normalized_data) <- c(sample_name_G,sample_name_R)

##### Create curated assay ####
xr <- new("EListRaw", list(E = GSE25534_Non_normalized_data, genes = GSE25534_raw$genes))
curatedExprs <- norm_probeToGenes_Agilent(xr, FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))

##### Create Column data #####
data_characteristic_G <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
characteristic_table_G <- sapply(1:length(data_characteristic_G[[1]]), function(x)
  sapply(data_characteristic_G, "[[", x))
characteristic_data_frame_G <- sub("(.*?): ", "", characteristic_table_G)

data_characteristic_R <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
  GEOquery::GSMList(gse)[[x]]@header$characteristics_ch2)
characteristic_table_R <- sapply(1:length(data_characteristic_R[[1]]), function(x)
  sapply(data_characteristic_R, "[[", x))
characteristic_data_frame_R <- sub("(.*?): ","",characteristic_table_R)

characteristic_data_frame <- rbind(characteristic_data_frame_G,
                                   characteristic_data_frame_R) %>% S4Vectors::DataFrame()
colnames(characteristic_data_frame) <- c("TBStatus", "Tissue")
row.names(characteristic_data_frame) <- c(sample_name_G, sample_name_R)
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)

for(i in 1:length(TBStatus)){
  if(TBStatus_temp[i] == "healthy non-infected donors"){
    TBStatus[i] = "Control"
  }
  if(TBStatus_temp[i] == "healthy donors latently infected with M. tuberculosis"){
    TBStatus[i] = "LTBI"
  }
  if(TBStatus_temp[i] == "tuberculosis patients"){
    TBStatus[i] = "PTB"
  }
}

characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$Tissue <- "Whole Blood"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
gpl1708 <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
gpl1708_annotate <- gpl1708@dataTable@table
gpl1708_annotate$ID <- as.character(gpl1708_annotate$ID)
GSE25534_probeNames <- data.frame(row.names(GSE25534_Non_normalized_data))
colnames(GSE25534_probeNames) <- c("ID_REF")

new_row_data <- GSE25534_probeNames %>%
  dplyr::left_join(gpl1708_annotate, by = c("ID_REF" = "ID")) %>%
  S4Vectors::DataFrame()
new_row_data$SYMBOL_NEW <- new_row_data$GENE_SYMBOL

##### Create Metadata #####
GSE25534_experimentData <- methods::new("MIAME",
                                        name = "Jeroen Maertzdorf",
                                        lab = "MPIIB",
                                        contact = "maertzdorf@mpiib-berlin.mpg.de",
                                        title = "Human gene expression profiles of susceptibility and resistance in tuberculosis.",
                                        abstract = "Tuberculosis (TB) still poses a profound burden on global health, owing to significant morbidity and mortality worldwide. Although a fully functional immune system is essential for the control of Mycobacterium tuberculosis infection, the underlying mechanisms and reasons for failure in part of the infected population remain enigmatic. Here, whole-blood microarray gene expression analyses were performed in TB patients and in latently as well as uninfected healthy controls to define biomarkers predictive of susceptibility and resistance. Fc gamma receptor 1B (FCGRIB)was identified as the most differentially expressed gene, and, in combination with four other markers, produced a high degree of accuracy in discriminating TB patients and latently infected donors. We determined differentially expressed genes unique for active disease and identified profiles that correlated with susceptibility and resistance to TB. Elevated expression of innate immune-related genes in active TB and higher expression of particular gene clusters involved in apoptosis and natural killer cell activity in latently infected donors are likely to be the major distinctive factors determining failure or success in controlling M. tuberculosis infection. The gene expression profiles defined in this study provide valuable clues for better understanding of progression from latent infection to active disease and pave the way for defining predictive correlates of protection in TB.",
                                        url = "10.1038/gene.2010.51",
                                        pubMedIds = "20861863",
                                        other = list(Platform = "Agilent-012391 Whole Human Genome Oligo Microarray G4112A (GPL1708)"))
GSE25534_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE25534_Non_normalized_data = as.matrix(GSE25534_Non_normalized_data)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE25534_experimentData));GSE25534_sobject
save_raw_files(GSE25534_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
