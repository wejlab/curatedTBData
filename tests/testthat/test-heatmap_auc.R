combine_dat_exp <- data.frame(Signature = rep(c("Anderson_42", "Anderson_OD_51",
                              "Berry_393", "Berry_OD_86", "Blankley_5"), 2),
                              AUC = stats::runif(10, 0.5, 1),
                              GSE = rep(c("GSE39939", "GSE19442"), each = 5))
test_that("Column with names \"Signature\", \"Study\", \"AUC\" cannot be missing.", {
  expect_error(heatmap_auc(combine_dat_exp))
})
