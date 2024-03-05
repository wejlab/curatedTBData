# Make meta-data
if (!require("magrittr", character.only = TRUE)) {
    BiocManager::install("magrittr")
    require("magrittr", character.only = TRUE)
}
dataNamesAll <- list.files("data-raw/", pattern = "*.RDS")
dataNames <- gsub("\\.RDS", "", dataNamesAll)
make_metadata <- function(dataName) {
    # dataCategory <- c("assay_raw", "assay_curated", "column_data", "row_data",
    #                   "meta_data", "assay_reprocess")
    dataNames_seq <- strsplit(dataName, "_") |>
        unlist()
    GSEName <- dataNames_seq[1]
    if (length(grep("assay_raw", dataName)) == 1) {
        Description <- paste("Raw transcriptome data derived from GEO accession:", GSEName)
        RDataClass <- "matrix"
        Tags <- paste(GSEName, "assay_raw", sep = ":")
    }
    if (length(grep("assay_curated", dataName)) == 1) {
        Description <- paste("Curated transcriptome data derived from GEO accession:", GSEName)
        RDataClass <- "matrix"
        Tags <- paste(GSEName, "assay_curated", sep = ":")
    }
    if (length(grep("column_data", dataName)) == 1) {
        Description <- paste("Clinical information for samples derived from GEO accession:", GSEName)
        RDataClass <- "DFrame"
        Tags <- paste(GSEName, "column_data", sep = ":")
    }
    if (length(grep("meta_data", dataName) == 1)) {
        Description <- paste("Meta data information for samples derived from GEO accession:", GSEName)
        RDataClass <- "MIAME"
        Tags <- paste(GSEName, "meta_data", sep = ":")
    }
    if (length(grep("row_data", dataName) == 1)) {
        Description <- paste("Probe set information for samples for derived from GEO accession:", GSEName)
        RDataClass <- "DFrame"
        Tags <- paste(GSEName, "row_data", sep = ":")
    }
    if (length(grep("assay_reprocess", dataName)) == 1) {
        Description <- paste("Reprocessed RNA-seq data derived from GEO accession:", GSEName)
        RDataClass <- "matrix"
        Tags <- paste(GSEName, paste0(dataNames_seq[-1], collapse = "_"), sep = ":")
    }
    ####################################################
    BiocVersion <- "3.18"
    Genome <- as.character(NA)
    SourceType <- rep(base::as.character("tar.gz"))
    SourceUrl <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSEName)
    DataProvider <- base::as.character("The National Center for Biotechnology Information")
    ####################################################
    SourceVersion  <- base::as.character(NA)
    ####################################################
    Species <- base::as.character("Homo sapiens")
    ####################################################
    TaxonomyId <- base::as.character("9606")
    ####################################################
    Coordinate_1_based <- base::as.logical(NA)
    ####################################################
    Maintainer <- base::as.character("Xutao Wang <xutaow@bu.edu>")
    ####################################################
    DispatchClass <- base::as.character("Rda")
    ####################################################
    RDataPath <- sprintf("curatedTBData/data/%s.rda", dataName)
    ####################################################
    metadata1 <- data.frame(Title = dataName, Description, BiocVersion, Genome, SourceType,
                            SourceUrl, SourceVersion, Species, TaxonomyId,
                            Coordinate_1_based, DataProvider, Maintainer, RDataClass,
                            DispatchClass, RDataPath, Tags)
    return(metadata1)
}
metadata_v1 <- lapply(dataNames, make_metadata) |>
    dplyr::bind_rows()
utils::write.csv(metadata_v1, "inst/extdata/metadata_v1.csv", row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata("curatedTBData/")
