returned_resources <- curatedTBData(c("GSE19435", "GSE19439"), dryrun = FALSE,
                                    curated.only = TRUE) %>%
  base::suppressWarnings()

test_that("Argument \"experment_name\" cannot be missing", {
  expect_error(combineObjects(returned_resources))
})

test_that("return type is SummarizedExperiment", {
  expect_s4_class(combineObjects(returned_resources, experiment_name = "assay_curated"),
                  "SummarizedExperiment")
})

test_that("Input list does not have unique name for each object,
          Input list contains only one element", {
  returned_resources_noname <- returned_resources
  names(returned_resources_noname) <- NULL
  expect_error(combineObjects(returned_resources_noname,
                              experiment_name = "assay_curated"))
})

test_that("return type is SummarizedExperiment", {
  expect_s4_class(combineObjects(returned_resources,
                                 experiment_name = "assay_curated"),
                  "SummarizedExperiment")
})
