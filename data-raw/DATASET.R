## code to prepare `DATASET` dataset goes here

# usethis::use_data("DATASET")
# usethis::use_data(NameOfTheDataset, compress = "xz", overwrite = TRUE)
library(purrr)
library(devtools)
dataNamesAll <- list.files("data-raw/", pattern = "*.RDS")
dataPathFull <- list.files("data-raw", pattern = "*.RDS", full.names = TRUE)
dataAllList <- lapply(dataPathFull, function(x) readRDS(x))
names(dataAllList) <- gsub("\\..*", "", dataNamesAll)
purrr::walk2(dataAllList, paste0(names(dataAllList)), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

# Save summary table
DataSummary_raw <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                  sheet = "DataSummary")
DataSummary <- DataSummary_raw[-base::nrow(DataSummary_raw),] %>% ## Remove the last row: Total
  base::as.data.frame()
usethis::use_data(DataSummary, compress = "xz", overwrite = TRUE)

# Save training table
SignatureInfoTraining_raw <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                            sheet = "SignatureInfoTraining")
SignatureInfoTraining <- SignatureInfoTraining_raw %>%
  dplyr::filter(.data$Study != "Unable to find training data") %>% ## Remove signatures with unidentified discovery datasets
  base::as.data.frame()
usethis::use_data(SignatureInfoTraining, compress = "xz", overwrite = TRUE)
