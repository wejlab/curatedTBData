context("data")

DataSummary <- get(utils::data("DataSummary", package = "curatedTBData"))

geo_sobject <- DataSummary %>% dplyr::as_tibble() %>%
  dplyr::filter(stringr::str_detect(.data$GeneralType, "Microarray|RT-PCR") ) %>%
  dplyr::select(.data$`GEO accession`)

geo_mobject <- DataSummary %>% dplyr::as_tibble() %>%
  dplyr::filter(stringr::str_detect(.data$GeneralType, "RNA-seq")) %>%
  dplyr::select(.data$`GEO accession`)

test_that("data are SummarizedExperiment objects", {

  lapply(geo_sobject, function(x) {
      expect_s4_class(get_curatedTBData(x)[[1]],"SummarizedExperiment")})

})

test_that("data are MultiAssayExperiment objects", {

  lapply(geo_mobject, function(x) {
           expect_s4_class(get_curatedTBData(x)[[1]],"MultiAssayExperiment")})

})
