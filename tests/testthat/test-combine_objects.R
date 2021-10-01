returned_resources <- curatedTBData(c("GSE107104", "GSE19435", "GSE19443"),
                                    dry.run = FALSE, curated.only = TRUE) |>
  suppressWarnings()

test_that("Argument \"experment_name\" cannot be missing", {
  expect_error(combine_objects(returned_resources))
})

test_that("return type is SummarizedExperiment", {
  re1 <- combine_objects(returned_resources, experiment_name = "assay_curated")
  re2 <- combine_objects(returned_resources, experiment_name = "assay_curated",
                         update_genes = FALSE)
  expect_s4_class(re1, "SummarizedExperiment")
  expect_s4_class(re2, "SummarizedExperiment")
})

test_that("Input list does not have unique name for each element.
          Input list contains only one element", {
  returned_resources_noname <- returned_resources
  names(returned_resources_noname) <- NULL
  expect_error(combine_objects(returned_resources_noname,
                               experiment_name = "assay_curated"))
  expect_error(combine_objects(returned_resources[1],
                               experiment_name = "assay_curated"))
})

