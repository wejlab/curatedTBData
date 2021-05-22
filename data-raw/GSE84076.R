if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
#### Read in raw data ####
geo <- "GSE84076"
urls <-  GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
temp <- tempfile()
download.file(as.character(urls$url[1]), temp)
GSE84076_raw <- read.delim(temp)
row_info <- c("gene_id","gene_short_name","locus")

GSE84076_normalized_counts <- GSE84076_raw %>%
  dplyr::select(-row_info) %>% as.matrix()
row.names(GSE84076_normalized_counts) <- GSE84076_raw$locus
GSE84076_normalized_curated <- GSE84076_raw %>%
  dplyr::select(-c("gene_short_name","locus")) %>%
  dplyr::group_by(gene_id) %>% dplyr::summarise_all(mean)
GSE84076_normalized_curated_avg <- as.matrix(GSE84076_normalized_curated[,-1])
row.names(GSE84076_normalized_curated_avg) <- GSE84076_normalized_curated$gene_id
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
description_id_raw <- sapply(1:length(names(GEOquery::GSMList(gse))),
                             function(x) GEOquery::GSMList(gse)[[x]]@header$title)
ID_table <- data.frame(SampleID = names(GEOquery::GSMList(gse)),
                       DescriptionID = description_id_raw)

indx <- match(ID_table$DescriptionID, colnames(GSE84076_normalized_counts))
GSE84076_normalized_counts <- GSE84076_normalized_counts[,indx]
colnames(GSE84076_normalized_counts) <- ID_table$SampleID

GSE84076_normalized_curated_avg <- GSE84076_normalized_curated_avg[,indx]
colnames(GSE84076_normalized_curated_avg) <- ID_table$SampleID

#### Create Column Data ####
characteristic_data_frame <- readRawColData(gse)
colnames(characteristic_data_frame) <- c("TBStatus", "Tissue")
characteristic_data_frame$Tissue <- "Whole Blood"

TBStatus_list <- strsplit(characteristic_data_frame$TBStatus, "-")
TBStatus <- TBStatus_temp <- unlist(lapply(TBStatus_list, function(x) x[1]),
                                    use.names = FALSE)
unique(TBStatus_temp)
for (i in 1:length(TBStatus_temp)) {
  if (TBStatus_temp[i] == "Active Tuberculosis") {
    TBStatus[i] = "PTB"
  } else if (TBStatus_temp[i] == "Treated Active Tuberculosis") {
    TBStatus[i] = "PTB"
  } else if (TBStatus_temp[i] == "Latent Tuberculosis ") {
    TBStatus[i] = "LTBI"
  } else if (TBStatus_temp[i] == "Control ") {
    TBStatus[i] = "Control"
  }
}
characteristic_data_frame$TBStatus <- TBStatus
BcgVaccinated <- BcgVaccinated_temp <- unlist(lapply(TBStatus_list, function(x) x[3]),
                                              use.names = FALSE)
index_filter <- which(!is.na(BcgVaccinated))
for (i in index_filter){
  if (BcgVaccinated_temp[i] == " vaccinated") {
    BcgVaccinated[i] = "Yes"
  }
  if (BcgVaccinated_temp[i] == " Unvaccinated") {
    BcgVaccinated[i] = "No"
  }
}
characteristic_data_frame$BcgVaccinated <- BcgVaccinated
TreatmentStatus <- rep(NA, nrow(characteristic_data_frame))
TreatmentStatus[which(TBStatus_temp == "Treated Active Tuberculosis")] <- "Treated"
TreatmentStatus[TBStatus_temp == "Active Tuberculosis"] <- "Treatment-naive"
characteristic_data_frame$TreatmentStatus <- TreatmentStatus
characteristic_data_frame$GeographicalRegion <- "Brazil"
characteristic_data_frame$HIVStatus <- "Negative"
col_info <- create_standard_coldata(characteristic_data_frame)
new_col_info <- S4Vectors::DataFrame(col_info)

#### Create Row Data ####
new_row_data <- GSE84076_raw %>%
  dplyr::select(row_info) %>%
  S4Vectors::DataFrame()
colnames(new_row_data) <- c("SYMBOL_NEW", "gene_short_name", "ID_REF")

##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Aravind Tallam",
                      lab = "TWINCORE",
                      contact = "aravind.tallam@twincore.de",
                      title = "Transcriptomic Biomarkers for Tuberculosis: Evaluation of DOCK9. EPHA4, and NPC2 mRNA Expression in Peripheral Blood.",
                      abstract = "Lately, much effort has been made to find mRNA biomarkers for tuberculosis (TB) disease/infection with microarray-based approaches. In a pilot investigation, through RNA sequencing technology, we observed a prominent modulation of DOCK9, EPHA4, and NPC2 mRNA abundance in the blood of TB patients. To corroborate these findings, independent validations were performed in cohorts from different areas. Gene expression levels in blood were evaluated by quantitative real-time PCR (Brazil, n = 129) or reanalysis of public microarray data (UK: n = 96; South Africa: n = 51; Germany: n = 26; and UK/France: n = 63). In the Brazilian cohort, significant modulation of all target-genes was observed comparing TB vs. healthy recent close TB contacts (rCt). With a 92% specificity, NPC2 mRNA high expression (NPC2high) showed the highest sensitivity (85%, 95% CI 65%-96%; area under the ROC curve [AUROC] = 0.88), followed by EPHA4 (53%, 95% CI 33%-73%, AUROC = 0.73) and DOCK9 (19%, 95% CI 7%-40%; AUROC = 0.66). All the other reanalyzed cohorts corroborated the potential of NPC2high as a biomarker for TB (sensitivity: 82-100%; specificity: 94-97%). An NPC2high profile was also observed in 60% (29/48) of the tuberculin skin test positive rCt, and additional follow-up evaluation revealed changes in the expression levels of NPC2 during the different stages of Mycobacterium tuberculosis infection, suggesting that further studies are needed to evaluate modulation of this gene during latent TB and/or progression to active disease. Considering its high specificity, our data indicate, for the first time, that NPC2high might serve as an accurate single-gene biomarker for TB.",
                      url = "10.3389/fmicb.2016.01586",
                      pubMedIds = "27826286",
                      other=list(Platform = "Illumina HiSeq 2500 (Homo sapiens) (GPL16791)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Normalized_counts = as.matrix(GSE84076_normalized_counts)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)
saveRDS(GSE84076_normalized_curated_avg, paste0("data-raw/", geo, "_assay_curated.RDS"))

#### Add reprocessed RNA-seq counts ####
GenomeVersion <- "hg19"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))
