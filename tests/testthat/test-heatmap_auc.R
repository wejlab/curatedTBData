mySig <- rep(c("Anderson_42", "Anderson_OD_51", "Berry_393",
               "Berry_OD_86", "Blankley_5"), 2)
myStudy <- rep(c("GSE39939", "GSE19442"), each = 5)
combine_dat_exp <- data.frame(Signature = mySig,
                              AUC = runif(10, 0.5, 1),
                              Study = myStudy)
GSE_sig_exp <- data.frame(TBSignature = mySig[1:4],
                          Study = c("GSE39939", "GSE39940",
                                    "GSE19442", "GSE19443"))
test_that("Column names: \"Signature\", \"Study\", \"AUC\" cannot be
          missing for \"combine_dat\".", {
    for (i in seq_len(ncol(combine_dat_exp))) {
        df <- combine_dat_exp[, -i, drop = FALSE]
        expect_error(heatmap_auc(combine_dat = df))
    }
})

test_that("Column names: \"TBSignature\", \"Study\" cannot be missing
          for \"GSE_sig\"", {
    for (i in seq_len(ncol(GSE_sig_exp))) {
        expect_error(heatmap_auc(combine_dat = combine_dat_exp,
                                 GSE_sig = GSE_sig_exp[, -i, drop = FALSE]))
    }
})
test_that("\"GSE_sig\" argument can be missing", {
    expect_message(heatmap_auc(combine_dat = combine_dat_exp, facet = TRUE))
})

test_that("Return type is a list of ggplot", {
    p <- heatmap_auc(combine_dat = combine_dat_exp, GSE_sig = GSE_sig_exp)
    expect_type(p, "list")
    expect_s3_class(p, "ggplot")
})
