test_that("Argument \"study_name\" cannot be missing.", {
  expect_error(curatedTBData(dry.run = FALSE, curated.only = FALSE))
  expect_error(curatedTBData(dry.run = TRUE, curated.only = FALSE))
  expect_error(curatedTBData(dry.run = FALSE, curated.only = TRUE))
  expect_error(curatedTBData(dry.run = TRUE, curated.only = TRUE))
})

test_that("Resources not available", {
  expect_error(curatedTBData(study_name = "Example", dry.run = TRUE,
                             curated.only = TRUE))
})

test_that("All resources from metadata exist in the ExperimentHub service", {
  metadata_file_path <- system.file("extdata/metadata.csv",
                                    package = "curatedTBData")
  if (!file.exists(metadata_file_path)) {
    metadata_file_path <- as.character("inst/extdata/metadata.csv")
  }
  returned_resources_names <- curatedTBData("", dry.run = TRUE,
                                            curated.only = FALSE) |>
    sort()
  metadata <- utils::read.csv(metadata_file_path)
  expect_equal(sort(metadata$Title), returned_resources_names)
  # return type is character when dry.run is `TRUE`
  expect_type(returned_resources_names, "character")
})

test_that("Retrun type is a list when dry.run is `FALSE`,
          all elements in the returned list is MultiAssayExperiment", {
  returned_resource <- curatedTBData("GSE74092", dry.run = FALSE,
                                      curated.only = TRUE) |>
    base::suppressWarnings()
  expect_type(returned_resource, "list")
  length_list <- length(returned_resource)
  for (i in seq_len(length_list)) {
    expect_s4_class(returned_resource[[i]], "MultiAssayExperiment")
  }
})
