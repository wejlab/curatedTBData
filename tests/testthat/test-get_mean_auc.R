returned_resources <- curatedTBData(c("GSE107104", "GSE19435"),
                                    dry.run = FALSE, curated.only = TRUE) |>
    suppressWarnings()
mysignatures <- list(Sweeney_OD_3 = c("DUSP3", "GBP5", "KLF2"),
                     Suliman_4 = c("ANKRD22", "C1QC", "OSBPL10", "TRAV27"))
sig_names <- names(mysignatures)
re1 <- lapply(returned_resources, function(x)
    subset_curatedTBData(x, annotationColName = "TBStatus",
                         annotationCondition = c("Control", "PTB")))
re2 <- lapply(re1, function(x)
    TBSignatureProfiler::runTBsigProfiler(input = x,
                                          useAssay = 1,
                                          signatures = mysignatures,
                                          algorithm = "ssGSEA",
                                          update_genes = FALSE))
df <- combine_auc(re2, annotationColName = "TBStatus",
                   signatureColNames = sig_names)

test_that("Use \"empirical\" as default method for
          Bootstrap Confidence Interval", {
    expect_message(get_mean_auc(df, column_name_variable = "Signature",
                                column_name_value = "AUC"))
})

test_that("Bootstrap Confidence Interval only support
          \"percentile\" and \"empirical\"", {
    expect_error(get_mean_auc(df, column_name_variable = "Signature",
                              column_name_value = "AUC", method = "fakeMethod"))
})
test_that("Return class is data.frame", {
    out <- get_mean_auc(df, column_name_variable = "Signature",
                        column_name_value = "AUC", method = "percentile")
    expect_s3_class(out, "data.frame")
})
