returned_resources <- curatedTBData(c("GSE107104", "GSE19435", "GSE19443"),
                                    dryrun = FALSE, curated.only = TRUE) %>%
  base::suppressWarnings()

test_that("Argument \"experment_name\" cannot be missing", {
  expect_error(combineObjects(returned_resources))
})

test_that("return type is SummarizedExperiment", {
  expect_s4_class(combineObjects(returned_resources, experiment_name = "assay_curated"),
                  "SummarizedExperiment")
})

test_that("Input list does not have unique name for each element.
          Input list contains only one element", {
  returned_resources_noname <- returned_resources
  names(returned_resources_noname) <- NULL
  expect_error(combineObjects(returned_resources_noname,
                              experiment_name = "assay_curated"))
  expect_error(combineObjects(returned_resources[1],
                              experiment_name = "assay_curated"))
})

test_that("return type is SummarizedExperiment", {
  expect_s4_class(combineObjects(returned_resources,
                                 experiment_name = "assay_curated"),
                  "SummarizedExperiment")
})
