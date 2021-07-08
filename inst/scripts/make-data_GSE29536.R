if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
geo <- "GSE29536"
sequencePlatform <- "GPL6102"
GSE29536_data_list <- readRawData(geo, sequencePlatform)

GSE29536_Non_normalized_data <- GSE29536_Non_pvalue <- GSE29536_data_list$data_Non_normalized
colnames(GSE29536_Non_normalized_data) <-
  colnames(GSE29536_Non_pvalue) <- gsub("\\..*", "", colnames(GSE29536_Non_pvalue))
##### Create Column data #####
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("disease_status", "gxbdataset_biirinternal")
index_TB <- grep("TUBERCULOSIS", characteristic_data_frame$gxbdataset_biirinternal)
characteristic_data_frame_sub <- characteristic_data_frame[index_TB,]
TBStatus <- ifelse(characteristic_data_frame_sub$disease_status == "PTB", "PTB", "Control")
BcgVaccinated <- rep(NA, nrow(characteristic_data_frame_sub))
index_bcg_pos <- which(characteristic_data_frame_sub$disease_status == "Control BCG+")
index_bcg_neg <- which(characteristic_data_frame_sub$disease_status == "Control BCG-")
BcgVaccinated[index_bcg_pos] <- "Yes"
BcgVaccinated[index_bcg_neg] <- "No"
col_info1 <- data.frame(TBStatus = TBStatus, BcgVaccinated = BcgVaccinated,
                       Tissue = "Whole Blood", GeographicalRegion = "UK",
                       HIVStatus = "Negative")
row.names(col_info1) <- row.names(characteristic_data_frame_sub)
col_info <- create_standard_coldata(col_info1)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
row_data <- map_gene_symbol(GSE29536_Non_pvalue, sequencePlatform)
new_row_data <- match_gene_symbol(row_data)

##### Create Metadata #####
GSE29536_experimentData <- methods::new("MIAME",
                                        name = "Damien Chaussabel",
                                        lab = "Baylor Institute for Immunology Research",
                                        contact = "DChaussabel@benaroyaresearch.org",
                                        title = "Whole Blood Transcriptional Modules generated on Illumina Hu-6 V2 Beadchips.",
                                        abstract = "This dataset was used to establish whole blood transcriptional modules (n=260) that represent groups of coordinately expressed transcripts that exhibit altered abundance within individual datasets or across multiple datasets. This modular framework was generated to reduce the dimensionality of whole blood microarray data processed on the Illumina Beadchip platform yielding data-driven transcriptional modules with biologic meaning.",
                                        url = "10.1371/journal.pone.0074893; 10.1016/j.immuni.2012.12.008; 10.1371/journal.pone.0034390",
                                        pubMedIds = "24069364;23601689;22496797",
                                        other = list(Platform = "Illumina human-6 v2.0 expression beadchip (GPL6102)"))
GSE29536_Non_normalized_data_sub <- GSE29536_Non_normalized_data[, row.names(new_col_info)]
GSE29536_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSE29536_Non_normalized_data = as.matrix(GSE29536_Non_normalized_data_sub)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSE29536_experimentData));GSE29536_sobject
save_raw_files(GSE29536_sobject, path = "data-raw/", geo = geo)
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

##### Create normalized curated assay #####
GSE29536_normed <- GSE29536_data_list$data_normalized
colnames(GSE29536_normed) <- gsub("\\..*", "", colnames(GSE29536_normed))
GSE29536_normed <- GSE29536_normed[, row.names(new_col_info)]
curatedExprs <- probesetsToGenes(row_data = new_row_data,
                                 data_normalized = GSE29536_normed,
                                 FUN = median)
saveRDS(curatedExprs, paste0("data-raw/", geo, "_assay_curated.RDS"))
