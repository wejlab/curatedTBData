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
                                        "Repeat_QGIT_Positive_Month",
                                        "Repeat_TST_Positive_Month",
                                        "sex", "SiteofTB","Tissue","visitmonth")]
TB_temp <- TBStatus <- India_metadata$SiteofTB
for (i in 1:length(TB_temp)){
  if (TB_temp[i] == "CONTROL") {
    TBStatus[i] = "Control"
  }
}
India_metadata$SiteofTB <- TBStatus
colnames(India_metadata)[colnames(India_metadata)=="SiteofTB"] <- "TBStatus"
colnames(India_metadata)[colnames(India_metadata)=="sex"] <- "Gender"
India_metadata$Gender <- ifelse(India_metadata$Gender == "female", "Female", "Male")



# assay_reprocess_final <- matchSRRtoSampleID(gse, assay_reprocess)
saveRDS(assay_reprocess_final, paste0("data-raw/", geo, "_assay_reprocess_",
                                      GenomeVersion,".RDS"))
