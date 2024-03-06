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
  sys_file_paths <- system.file("extdata", package = "curatedTBData")
  metadata_file_paths <- list.files(sys_file_paths, pattern = "metadata",
                                    full.names = TRUE)
  if (length(metadata_file_paths) == 0) { # For development
    metadata_file_paths <- list.files("inst/extdata", pattern = "metadata",
                                      full.names = TRUE)
  }
  returned_resources_names <- curatedTBData("", dry.run = TRUE,
                                            curated.only = FALSE) |>
    sort()
  metadata_tab_list <- lapply(metadata_file_paths, utils::read.csv)
  metadata <- do.call(rbind, metadata_tab_list)
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
