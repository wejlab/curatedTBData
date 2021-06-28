if (!require("magrittr", character.only = TRUE)) {
  BiocManager::install("magrittr")
  require("magrittr", character.only = TRUE)
}
source("data-raw/UtilityFunctionForCuration.R")

##### Read in raw data #####
sequencePlatform <- "GPL10558"
GSEBruno <- read.delim("data-raw/expressionset_normalized_used_for_limma_edots.tsv",
                       stringsAsFactors = FALSE, row.names = 1)
GSEBruno_col_info <- read.delim("data-raw/sample_annotation_edots.tsv",
                                 stringsAsFactors = FALSE)

TB_status <- TBStatus_temp <- GSEBruno_col_info$class
unique(TBStatus_temp)
DMStatus <- rep(0, nrow(GSEBruno_col_info))
for (i in 1:length(TBStatus_temp)) {
  if (TBStatus_temp[i] == "DM_TB" | TBStatus_temp[i] == "NonDM_TB") {
    TB_status[i] = "PTB"
  }

  if (TBStatus_temp[i] == "NonDM_NonTB") {
    TB_status[i] = "Control"
  }
  if (TBStatus_temp[i] == "DM_NonTB") {
    TB_status[i] = "OD"
  }
}

for (i in 1:length(DMStatus)) {
  if (TBStatus_temp[i] == "NonDM_NonTB" | TBStatus_temp[i] == "NonDM_TB") {
    DMStatus[i] = "Negative"
  }
  else {
    DMStatus[i] = "Positive"
  }
}
GSEBruno_col_info$TBStatus <- TB_status
GSEBruno_col_info$DiabetesStatus <- DMStatus
row.names(GSEBruno_col_info) <- GSEBruno_col_info$sample
# Matching sample names to annotation file
index <- match(colnames(GSEBruno), GSEBruno_col_info$sample)
col_info <- GSEBruno_col_info[index, ]
col_info <- col_info[, -c(1:2)]

col_info <- create_standard_coldata(col_info)
new_col_info <- S4Vectors::DataFrame(col_info)

new_row_data <- S4Vectors::DataFrame(ID_REF = row.names(GSEBruno),
                                 SYMBOL_NEW = row.names(GSEBruno))
GSEBruno_experimentData <- new("MIAME",
                               name = "Cesar A. Prada-Medina",
                               lab = "Brazail",
                               contact = "hnakaya@usp.br",
                               title = "Systems Immunology of Diabetes- Tuberculosis Comorbidity
Reveals Signatures of Disease Complications",
                               abstract = "Comorbid diabetes mellitus (DM) increases tuberculosis (TB) risk and adverse outcomes but the pathological interactions between DM and TB remain incompletely understood. We performed an integrative analysis of whole blood gene expression and plasma analytes, comparing South Indian
TB patients with and without DM to diabetic and non-diabetic controls without TB. Luminex assay
of plasma cytokines and growth factors delineated a distinct biosignature in comorbid TBDM in this cohort. Transcriptional profiling revealed elements in common with published TB signatures from cohorts that excluded DM. Neutrophil count correlated with the molecular degree of perturbation, especially in TBDM patients. Body mass index and HDL cholesterol were negatively correlated with molecular degree of perturbation. Diabetic complication pathways including several pathways linked
to epigenetic reprogramming were activated in TBDM above levels observed with DM alone. Our data provide a rationale for trials of host-directed therapies in TBDM, targeting neutrophilic inflammation and diabetic complication pathways to address the greater morbidity and mortality associated with this increasingly prevalent dual burden of communicable and non-communicable diseases.",
                               url = "10.1038/s41598-017-01767-4",
                               pubMedIds = "28515464",
                               other=list(Platform = "Illumina Human HT-12 v4 Expression BeadChips (GPL10558)"))
GSEBruno_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(GSEBruno_normalized_data= as.matrix(GSEBruno)),
  colData = new_col_info,
  rowData = new_row_data,
  metadata = list(GSEBruno_experimentData));GSEBruno_sobject

save_raw_files(GSEBruno_sobject, path = "data-raw/", geo = "GSEBruno")
saveRDS(GSEBruno, paste0("data-raw/", "GSEBruno", "_assay_curated.RDS"))

