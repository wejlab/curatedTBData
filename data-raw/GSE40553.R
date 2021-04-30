if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE40553"
sequencePlatform <- "GPL10558"
GSE40553_Non_pvalues_SA <- readRawData(geo, sequencePlatform, urlIndex = 4)
GSE40553_Non_pvalues_UK <- readRawData(geo, sequencePlatform, urlIndex = 5)

##### Match colnames to sample ID #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
description_id_UK_ori <- description_id_raw[1:ncol(GSE40553_Non_pvalues_UK)]
description_id_UK <- gsub("-", ".", description_id_UK_ori)
decription_id_SA_ori <- gsub(".*[[]([^]]+)[]].*", "\\1",
                             description_id_raw[(ncol(GSE40553_Non_pvalues_UK)+1):length(description_id_raw)])
decription_id_SA <- gsub("/", ".", decription_id_SA_ori)

description_id <- c(description_id_UK, decription_id_SA)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id)
# Merge two datasets
GSE40553_Non_pvalues_ori <- merge(GSE40553_Non_pvalues_SA, GSE40553_Non_pvalues_UK,
                                  by = "row.names", all = TRUE)
row.names(GSE40553_Non_pvalues_ori) <- GSE40553_Non_pvalues_ori$Row.names
GSE40553_Non_pvalues <- GSE40553_Non_pvalues_ori %>% dplyr::select(-c("Row.names"))
colnames(GSE40553_Non_pvalues) <- gsub("X", "", colnames(GSE40553_Non_pvalues))
indx <- base::match(ID_table$DescriptionID, colnames(GSE40553_Non_pvalues))
GSE40553_Non_pvalues <- GSE40553_Non_pvalues[,indx]
colnames(GSE40553_Non_pvalues) <- ID_table$SampleID
GSE40553_Non_normalized_counts <- GSE40553_Non_pvalues

##### Create column data #####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "MeasurementTime", "Tissue",
                                         "GeographicalRegion")
TBStatus <- TBStatus_temp <- as.character(characteristic_data_frame$TBStatus)
unique(TBStatus_temp)
TBStatus <- ifelse(TBStatus_temp == "LTB", "LTBI",TBStatus_temp)
characteristic_data_frame$TBStatus <- TBStatus
characteristic_data_frame$GeographicalRegion <- gsub(" patient.*", "" ,
                                                     characteristic_data_frame$GeographicalRegion)
characteristic_data_frame$PatientID <- gsub("\\..*", "", ID_table$DescriptionID)
MeasurementTime <- paste(characteristic_data_frame$MeasurementTime, "Month(s)")
characteristic_data_frame$MeasurementTime <- MeasurementTime
characteristic_data_frame$Tissue <- "Whole Blood"
characteristic_data_frame$HIVStatus <- "Negative"
Treatment <- rep(0, nrow(characteristic_data_frame))
indexLTBI <- grep("LTBI", characteristic_data_frame$TBStatus)
indexPTB <- grep("PTB", characteristic_data_frame$TBStatus)
indexBaseline <- which(characteristic_data_frame$MeasurementTime == "0 Month(s)")
Treatment[grep("LTBI", characteristic_data_frame$TBStatus)] <- "Untreated"
Treatment[intersect(indexPTB, indexBaseline)] <- "Pre-treatment"
Treatment <- ifelse(Treatment == 0, "anti-TB treatment", Treatment)
characteristic_data_frame$Treatment <- Treatment
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

###### Create Row data #####
row_data <- map_gene_symbol(GSE40553_Non_pvalues, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

###### Create Metadata #####
GSE40553_experimentData <- methods::new("MIAME",
                                        name = "Chole Bloom",
                                        lab = "MRC National Institute for Medical Research",
                                        contact = "cbloom@nimr.mrc.ac.uk",
                                        title = "Detectable changes in the blood transcriptome are present after two weeks of antituberculosis therapy.",
                                        abstract = "Analysis of UK blood transcriptional profiles before treatment to indentify changes that occur during (2 weeks, 2 months), at the end of treatment (6 months) and after treatment (12 months)
Analysis of South African blood transcriptional profiles before treatment to indentify changes that occur during (2 weeks, 2 months), at the end of treatment (6 months) and after treatment (12 months)",
                                        url = "10.1371/journal.pone.0046191",
                                        pubMedIds = "23056259",
                                        other = list(Platform = "Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)"))
GSE40553_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE40553_Non_normalized_counts= as.matrix(GSE40553_Non_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE40553_experimentData))
save_raw_files(GSE40553_sobject, path = "data-raw/", geo = geo)

# Remove files in temporary directory
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)


