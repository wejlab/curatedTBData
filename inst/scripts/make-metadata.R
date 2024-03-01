# Make meta-data
if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
####################################################
data("DataSummary")
TitleAll <- DataSummary$GEOAccession[grep("GSE", DataSummary$GEOAccession)]
# Check the 4th character
indexGEO <- which(suppressWarnings(!is.na(as.numeric(substring(TitleAll, 4, 4)))))
TitleGEOMicroarray <- DataSummary[indexGEO, ] %>%
  dplyr::filter(GeneralType != "Illumina RNA-seq") %>%
  dplyr::select(GEOAccession) %>% unlist(use.names = FALSE)
TitleGEORNASeq <- DataSummary[indexGEO, ] %>%
  dplyr::filter(GeneralType == "Illumina RNA-seq") %>%
  dplyr::select(GEOAccession) %>% unlist(use.names = FALSE)
####################################################
#' @param GSEName The name of the GEO accession.
#' @param isGEO Boolean. Indicate whether the data is downloaded from the Gene Expression Omnibus.
#' @param DataType Whether it is a RNA sequencing or microarray data.
#' @param containReprocess Boolean. Indicate whether the RNA-seq data contains reprocessing count matrix
#' @param reprocessVersion String. The genome reference version specifically for reprocessing RNA-seq.

createMetaData <- function(GSEName, isGEO = TRUE, dataType = c("RNA-seq", "Microarray"),
                           containReprocess = FALSE, reprocessVersion = NULL) {
  dataCategory <- c("assay_raw", "assay_curated", "column_data", "row_data", "meta_data")
  dataType <- match.arg(dataType)
  Title <- paste0(GSEName, sep = "_", dataCategory)
  if (dataType == "RNA-seq") {
    if (containReprocess) {
      if (!is.null(reprocessVersion)) {
        assay_reprocess <- paste0("assay_reprocess_", reprocessVersion)
      } else {
        stop("check your reprocessVersion paratmeter")
      }
      dataCategory <- c(dataCategory, assay_reprocess)
      Title <- paste0(GSEName, sep = "_", dataCategory)
    }
  }
  Description <- rep(0, length(Title))
  RDataClass <- rep(0, length(Title))
  Tags <- rep(0, length(Title)) # Add tags from review
  for (i in seq_len(length(Title))) {
    if (isGEO) {
      if(length(grep("assay_raw", Title[i])) == 1) {
        Description[i] <- paste("Raw transcriptome data derived from GEO accession:", GSEName)
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_raw", sep = ":")
      }
      if (length(grep("assay_curated", Title[i])) == 1) {
        Description[i] <- paste("Curated transcriptome data derived from GEO accession:", GSEName)
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_curated", sep = ":")
      }
      if (length(grep("column_data", Title[i])) == 1) {
        Description[i] <- paste("Clinical information for samples derived from GEO accession:", GSEName)
        RDataClass[i] <- "DFrame"
        Tags[i] <- paste(GSEName, "column_data", sep = ":")
      }
      if (length(grep("meta_data", Title[i]) == 1)) {
        Description[i] <- paste("Meta data information for samples derived from GEO accession:", GSEName)
        RDataClass[i] <- "MIAME"
        Tags[i] <- paste(GSEName, "meta_data", sep = ":")
      }
      if (length(grep("row_data", Title[i]) == 1)) {
        Description[i] <- paste("Probe set information for samples for derived from GEO accession:", GSEName)
        RDataClass[i] <- "DFrame"
        Tags[i] <- paste(GSEName, "row_data", sep = ":")
      }
      if (length(grep("assay_reprocess", Title[i])) == 1) {
        Description[i] <- paste("Reprocessed RNA-seq data derived from GEO accession:", GSEName)
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_reprocess", sep = ":")
      }
    } else {
      GSE0 <- gsub("GSE","", GSEName)
      if (length(grep("assay_raw", Title[i])) == 1) {
        Description[i] <- paste("Raw transcriptome data derived from", GSE0, "et.al")
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_raw", sep = ":")
      }
      if (length(grep("assay_curated", Title[i])) == 1) {
        Description[i] <- paste("Curated transcriptome data derived from", GSE0, "et.al")
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_curated", sep = ":")
      }
      if (length(grep("column_data", Title[i])) == 1) {
        Description[i] <- paste("Clinical information for samples derived from", GSE0, "et.al")
        RDataClass[i] <- "DFrame"
        Tags[i] <- paste(GSEName, "column_data", sep = ":")
      }
      if (length(grep("meta_data", Title[i]) == 1)) {
        Description[i] <- paste("Meta data information for samples derived from", GSE0, "et.al")
        RDataClass[i] <- "MIAME"
        Tags[i] <- paste(GSEName, "meta_data", sep = ":")
      }
      if (length(grep("row_data", Title[i]) == 1)) {
        Description[i] <- paste("Probe set information for samples for derived from", GSE0, "et.al")
        RDataClass[i] <- "DFrame"
        Tags[i] <- paste(GSEName, "row_data", sep = ":")
      }
      if (length(grep("assay_reprocess", Title[i])) == 1) {
        Description[i] <- paste("Reprocessed RNA-seq data derived from", GSE0, "et.al")
        RDataClass[i] <- "matrix"
        Tags[i] <- paste(GSEName, "assay_reprocess", sep = ":")
      }
    }
  }
  ####################################################
  # BiocVersion <-
  #   BiocManager::version() %>%
  #   base::as.character()
  BiocVersion <- "3.14" # Edit based on review
  ####################################################
  Genome <- as.character(NA)
  ####################################################
  if (isGEO) {
    SourceType <- rep(base::as.character("tar.gz"))
    SourceUrl <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSEName)
    DataProvider <- base::as.character("The National Center for Biotechnology Information")
  } else {
    SourceType <- "To be added manually"
    SourceUrl <- "To be added manually"
    DataProvider <- "To be added manually"
  }
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
  RDataPath <- sprintf("curatedTBData/data/%s.rda", Title)
  ####################################################
  # Tags <- base::as.character("curatedTBData")
  ####################################################
  metadata1 <- data.frame(Title, Description, BiocVersion, Genome, SourceType,
                          SourceUrl, SourceVersion, Species, TaxonomyId,
                          Coordinate_1_based, DataProvider, Maintainer, RDataClass,
                          DispatchClass, RDataPath, Tags)
  return(metadata1)
}

metadataMicroarrayList <- lapply(TitleGEOMicroarray, function(x)
  createMetaData(GSEName = x, isGEO = TRUE, dataType = "Microarray"))
metadataMicroarray <- do.call(rbind, metadataMicroarrayList)

reprocssVersion <- rep(0, length(TitleGEORNASeq))
reprocssVersion <- ifelse(TitleGEORNASeq %in% paste0("GSE10799", c(1:4)), "hg38", "hg19")
metadataRNASeqList <- lapply(seq_len(length(TitleGEORNASeq)), function(i)
  createMetaData(GSEName = TitleGEORNASeq[i], isGEO = TRUE, dataType = "RNA-seq",
                 containReprocess = TRUE, reprocessVersion = reprocssVersion[i]))
metadataRNASeq <- do.call(rbind, metadataRNASeqList)

GSEBruno <- TitleAll[-indexGEO][1]
metadataBruno <- createMetaData(GSEName = GSEBruno, isGEO = FALSE, dataType = "Microarray")
metadataBruno$SourceType <- "TSV"
metadataBruno$SourceUrl <- "https://pubmed.ncbi.nlm.nih.gov/28515464/"
metadataBruno$DataProvider <- "Bruno B Andrade"

GSETornheim <- TitleAll[-indexGEO][2]
metadataTornheim <- createMetaData(GSEName = GSETornheim, isGEO = FALSE, dataType = "RNA-seq",
                                containReprocess = FALSE)
metadataTornheim$SourceType <- "FASTQ"
metadataTornheim$SourceUrl <- "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP229386&o=acc_s%3Aa"
metadataTornheim$DataProvider <- "The National Center for Biotechnology Information"

metadata <- rbind(metadataMicroarray, metadataRNASeq, metadataBruno, metadataTornheim)
utils::write.csv(metadata, "inst/extdata/metadata.csv", row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata("curatedTBData/")
