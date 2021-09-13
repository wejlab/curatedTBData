test_that("Argument \"study_name\" cannot be missing.", {
  expect_error(curatedTBData(dryrun = FALSE, curated.only = FALSE))
  expect_error(curatedTBData(dryrun = TRUE, curated.only = FALSE))
  expect_error(curatedTBData(dryrun = FALSE, curated.only = TRUE))
  expect_error(curatedTBData(dryrun = TRUE, curated.only = TRUE))
})

test_that("Resources not available", {
  expect_error(curatedTBData(study_name = "Example", dryrun = TRUE,
                             curated.only = TRUE))
})

test_that("All resources from metadata exist in the ExperimentHub service", {
  metadata_file_path <-
    base::system.file("extdata/metadata.csv", package = "curatedTBData")
  if (!base::file.exists(metadata_file_path)) {
    metadata_file_path <-
      base::as.character("inst/extdata/metadata.csv")
  }
  returned_resources_names <-
    curatedTBData("", dryrun = TRUE, curated.only = FALSE) %>%
    base::sort()
  metadata <- utils::read.csv(metadata_file_path)
  expect_equal(base::sort(metadata$Title), returned_resources_names)
  # return type is character when dryrun is `TRUE`
  expect_type(returned_resources_names, "character")
})

test_that("Retrun type is a list when dryrun is `FALSE`,
          all elements in the returned list is MultiAssayExperiment", {
  returned_resource <- curatedTBData("GSE74092", dryrun = FALSE,
                                      curated.only = TRUE) %>%
    base::suppressWarnings()
  expect_type(returned_resource, "list")
  length_list <- base::length(returned_resource)
  for (i in base::seq_len(length_list)) {
    expect_s4_class(returned_resource[[i]], "MultiAssayExperiment")
  }
})
