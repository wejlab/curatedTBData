context("test-MultiAsssay/TBSig")

# Create a Summarized Experiment Object
library(TBSignatureProfiler)

set.seed(123)

test_count <- rbind(matrix(rnorm(10*length(TBsignatures$Walter_51)), length(TBsignatures$Walter_51), 10,
                             dimnames = list(paste0("probe", 1:length(TBsignatures$Walter_51)),
                                             paste0("sample", 1:10))),
                      matrix(rnorm(10*length(TBsignatures$Walter_PNA_119)), length(TBsignatures$Walter_PNA_119), 10,
                             dimnames = list(paste0("probe", (length(TBsignatures$Walter_51)+1):(length(TBsignatures$Walter_51)+length(TBsignatures$Walter_PNA_119))),
                                             paste0("sample", 1:10))))
test_row_data <- data.frame(ID_REF = paste0("probe", 1:(length(TBsignatures$Walter_51)+length(TBsignatures$Walter_PNA_119))),
                            SYMBOL_NEW = c(TBsignatures$Walter_51,TBsignatures$Walter_PNA_119))
test_col_data <- data.frame(TBStatus=c(rep("PTB",3),rep("Latent",3),rep("OD",4)))
row.names(test_col_data) <- colnames(test_count)

test_sobject <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts= as.matrix(test_count)), rowData = DataFrame(test_row_data), colData = DataFrame(test_col_data))

# test_TBsig_sobject <- get_sobject_TBSig(test_mobject,"PTB", "Latent")
test_ssgsea_TBsig_sobject <- TBSignatureProfiler::runTBsigProfiler(input = test_TBsig_sobject,
                                                                   useAssay = "counts",
                                                                   signatures = TBsignatures,
                                                                   algorithm = "ssGSEA",
                                                                   combineSigAndAlgorithm = TRUE,
                                                                   parallel.sz = 1)
# In TBSignatureProfiler::runTBsigProfiler(input = test_TBsig_sobject,  :
# No identifiers in the gene sets could be matched to the identifiers
#in the expression data for the following signatures:  Blankley_5, Lee_4, Maertzdorf_4, Roe_OD_4, Sloot_HIV_2, Suliman_RISK_4, Thompson_FAIL_13, Thompson_RES_5


# Test for plot
test_ssgsea_auc <- get_pvalue_auc(test_ssgsea_TBsig_sobject, annotationColName = "Disease",
                                  signatureColNames = c("Anderson_42", "Anderson_OD_51", "Berry_393"))



