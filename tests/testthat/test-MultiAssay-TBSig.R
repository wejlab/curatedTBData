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






