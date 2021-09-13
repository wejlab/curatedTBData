returned_resource <- curatedTBData("GSE74092", dryrun = FALSE,
                                    curated.only = TRUE) %>%
  base::suppressMessages()
mobject <- returned_resource[[1]]
assay_curated <- returned_resource$GSE74092[["assay_curated"]]
col_info <- SummarizedExperiment::colData(mobject)
sobject <- SummarizedExperiment::SummarizedExperiment(list(counts = assay_curated),
                                                      colData = col_info)

test_that("stop when \"annotationColName\" is not found in the metadata", {
  expect_error(subset_curatedTBData(mobject,
                                    annotationColName = "FakeColumnNames",
                                    annotationCondition = c("Control", "PTB"),
                                    assayName = "assay_curated"))
})

test_that("return type is NULL if at least one of the conditions is not found from the column data", {
  expect_null(subset_curatedTBData(mobject,
                                   annotationColName = "TBStatus",
                                   annotationCondition = c("Control", "LTBI"),
                                   assayName = "assay_curated"))
  expect_null(subset_curatedTBData(sobject,
                                   annotationColName = "TBStatus",
                                   annotationCondition = c("Control", "LTBI"),
                                   assayName = "counts"))
})

test_that("stop when \"assayName\" is not found from the input", {
  expect_error(subset_curatedTBData(mobject,
                                    annotationColName = "TBStatus",
                                    annotationCondition = c("Control", "PTB"),
                                    assayName = "EXAMPLE_ASSAY"))
  expect_error(subset_curatedTBData(sobject,
                                    annotationColName = "TBStatus",
                                    annotationCondition = c("Control", "PTB"),
                                    assayName = "EXAMPLE_ASSAY"))
})
test_that("retrun type is SummarizedExperiment", {
  expect_s4_class(subset_curatedTBData(mobject,
                                       annotationColName = "TBStatus",
                                       annotationCondition = c("Control", "PTB"),
                                       assayName = "assay_curated"),
                  "SummarizedExperiment")
  expect_s4_class(subset_curatedTBData(sobject,
                                       annotationColName = "TBStatus",
                                       annotationCondition = c("Control", "PTB"),
                                       assayName = "counts"),
                  "SummarizedExperiment")
})
