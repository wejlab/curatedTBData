if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")
#### Only reprocessed RNA-seq counts ####
geo <- "GSETornheim"
GenomeVersion <- "hg38"
assay_reprocess <- read.delim(paste0("~/Desktop/practice/ReprocessRNASeqCounts/",
                                     GenomeVersion,"/", geo,"_reprocess_",
                                     GenomeVersion,".txt"))
# saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
#                                       GenomeVersion,".RDS"))
##### Read in Column data #####
# Download from website: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP229386&o=acc_s%3Aa&s=SRR10424375,SRR10424376,SRR10424377,SRR10424378,SRR10424379,SRR10424380,SRR10424381,SRR10424382,SRR10424383,SRR10424384,SRR10424385,SRR10424386,SRR10424387,SRR10424388,SRR10424389,SRR10424390,SRR10424391,SRR10424392,SRR10424393,SRR10424394,SRR10424395,SRR10424396,SRR10424397,SRR10424398,SRR10424399,SRR10424400,SRR10424401,SRR10424402,SRR10424403,SRR10424404,SRR10424405,SRR10424406,SRR10424407,SRR10424408,SRR10424409,SRR10424410,SRR10424411,SRR10424412,SRR10424413,SRR10424414,SRR10424415,SRR10424416,SRR10424417,SRR10424418,SRR10424419,SRR10424420,SRR10424421,SRR10424422,SRR10424423,SRR10424424,SRR10424425,SRR10424426,SRR10424427,SRR10424428,SRR10424429,SRR10424430,SRR10424431,SRR10424432,SRR10424433,SRR10424434,SRR10424435,SRR10424436,SRR10424437,SRR10424438,SRR10424439,SRR10424440,SRR10424441,SRR10424442,SRR10424443,SRR10424444,SRR10424445,SRR10424446,SRR10424447,SRR10424448,SRR10424449,SRR10424450,SRR10424451,SRR10424452,SRR10424453,SRR10424454,SRR10424455,SRR10424456,SRR10424457,SRR10424458,SRR10424459,SRR10424461,SRR10424462,SRR10424463,SRR10424464,SRR10424465,SRR10424466,SRR10424467,SRR10424468,SRR10424469,SRR10424470,SRR10424471,SRR10424472,SRR10424474,SRR10424460,SRR10424473
India_metadata_raw <- read.csv("data-raw/SraRunTable.txt")
row.names(India_metadata_raw) <- India_metadata_raw$Run

sample_baseline <- India_metadata_raw[,c("Run","Isolate","SiteofTB","visitmonth")] %>%
  dplyr::arrange(visitmonth) %>%
  data.frame() %>%
  dplyr::group_by(Isolate) %>%
  dplyr::mutate(first = dplyr::first(Run))
aa <- India_metadata_raw[unique(sample_baseline$first),]
aa$SiteofTB %>% table()
# 32 controls and 16 cases, same with paper
aa$controltype %>% table()
# 13 convertors, same with paper
India_metadata <- India_metadata_raw[,c("AGE","controltype","ConvertMonth",
                                        "Isolate",
                                        "PAXGene_Processing_Batch" ,
                                        "sex", "SiteofTB","Tissue","visitmonth")]
TB_temp <- TBStatus <- India_metadata$SiteofTB
for (i in 1:length(TB_temp)){
  if (TB_temp[i] == "CONTROL") {
    TBStatus[i] = "Control"
  }
}
India_metadata$SiteofTB <- TBStatus
colnames(India_metadata)[colnames(India_metadata) == "PAXGene_Processing_Batch"] <- "Batch"
colnames(India_metadata)[colnames(India_metadata) == "SiteofTB"] <- "TBStatus"
colnames(India_metadata)[colnames(India_metadata) == "sex"] <- "Gender"
India_metadata$Gender <- ifelse(India_metadata$Gender == "female", "Female", "Male")
colnames(India_metadata)[colnames(India_metadata) == "isolate"] <- "PatientID"
colnames(India_metadata)[colnames(India_metadata) == "visitmonth"] <- "MeasurementTime"
India_metadata$MeasurementTime <- paste(India_metadata$MeasurementTime, "Month(s)")
colnames(India_metadata)[colnames(India_metadata) == "AGE"] <- "Age"
India_metadata$Age <- as.numeric(India_metadata$Age)
colnames(India_metadata)[colnames(India_metadata) == "controltype"] <- "Progression"
progress_info <- India_metadata$Progression
progress_info[grep("nonconverter", India_metadata$Progression)] <- "Negative"
progress_info[grep("converter", progress_info)] <- "Positive"
India_metadata$Progression <- progress_info
colnames(India_metadata)[colnames(India_metadata) == "tissue"] <- "Tissue"
India_metadata$Tissue <- "Whole Blood"

colnames(India_metadata)[colnames(India_metadata) == "ConvertMonth"] <- "TimeToTB"
timeToTB <- India_metadata$TimeToTB
timeToTB[!is.na(timeToTB)] <- paste(timeToTB[!is.na(timeToTB)], "Month(s)")
India_metadata$TimeToTB <- timeToTB
India_metadata$GeographicalRegion <- "India"
India_metadata$HIVStatus <- "Negative"

col_info <- create_standard_coldata(India_metadata)
new_col_info <- S4Vectors::DataFrame(col_info)

##### Create Row Data #####
new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(assay_reprocess),
                                     SYMBOL_NEW = row.names(assay_reprocess))
##### Create Metadata #####
##### Create Metadata #####
experimentData <- new("MIAME",
                      name = "Jeffrey A. Tornheim",
                      lab = "Johns Hopkins University School of Medicine",
                      contact = "tornheim@jhu.edu",
                      title = "Transcriptomic Profiles of Confirmed Pediatric Tuberculosis Patients and Household Contacts Identifies Active Tuberculosis, Infection, and Treatment Response Among Indian Children",
                      abstract = "Background. Gene expression profiling is emerging as a tool for tuberculosis diagnosis and treatment response monitoring, but limited data specific to Indian children and incident tuberculosis infection (TBI) exist.
Methods. Sixteen pediatric Indian tuberculosis cases were age- and sex-matched to 32 tuberculosis-exposed controls (13 de- veloped incident TBI without subsequent active tuberculosis). Longitudinal samples were collected for ribonucleic acid sequencing. Differential expression analysis generated gene lists that identify tuberculosis diagnosis and tuberculosis treatment response. Data were compared with published gene lists. Population-specific risk score thresholds were calculated.
Results. Seventy-one genes identified tuberculosis diagnosis and 25 treatment response. Within-group expression was partially explained by age, sex, and incident TBI. Transient changes in gene expression were identified after both infection and treatment. Application of 27 published gene lists to our data found variable performance for tuberculosis diagnosis (sensitivity 0.38–1.00, specificity 0.48–0.93) and treatment response (sensitivity 0.70–0.80, specificity 0.40–0.80). Our gene lists found similarly variable performance when applied to published datasets for diagnosis (sensitivity 0.56–0.85, specificity 0.50–0.85) and treatment response (sensitivity 0.49– 0.86, specificity 0.50–0.84).
Conclusions. Gene expression profiles among Indian children with confirmed tuberculosis were distinct from adult-derived gene lists, highlighting the importance of including distinct populations in differential gene expression models.",
                      url = "doi.org/10.1093/infdis/jiz639",
                      pubMedIds = "31796955",
                      other=list(Platform = "Illumina HiSeq 2500 (Homo sapiens) (GPL16791)"))
sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(Non_normalized_counts = as.matrix(assay_reprocess)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(experimentData));sobject
save_raw_files(sobject, path = "data-raw/", geo = geo)

##### Create curated matched assay #####
counts <- assay_reprocess
counts[counts < 5] <- 5
counts <- log(counts, base=2)
NormFactor <- edgeR::calcNormFactors(counts, method = "TMM")
ScaleFactors <- colSums(counts) * NormFactor
counts_normalized <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))

saveRDS(counts_normalized, paste0("data-raw/", geo, "_assay_curated.RDS"))


