returned_resource <- curatedTBData("GSE74092", dryrun = FALSE,
                                    curated.only = TRUE) %>%
  base::suppressWarnings()
test_that("return type is NULL if at least one of the conditions is not found", {
  returned_subset <- subset_curatedTBData(returned_resource[[1]],
                                             annotationColName = "TBStatus",
                                             annotationCondition = c("Control","LTBI"),
                                             assayName = "assay_curated")
  expect_null(returned_subset)
})

test_that("return type is NULL if annotationColName is not found in the metadata", {
  returned_subset <- subset_curatedTBData(returned_resource[[1]],
                                          annotationColName = "FakeColumnNames",
                                          annotationCondition = c("Control","LTBI"),
                                          assayName = "assay_curated")
  expect_null(returned_subset)
})

test_that("retrun type is SummarizedExperiment", {
  returned_subset <- subset_curatedTBData(returned_resource[[1]],
                                          annotationColName = "TBStatus",
                                          annotationCondition = c("Control","PTB"),
                                          assayName = "assay_curated")
  expect_s4_class(returned_subset, "SummarizedExperiment")
})
