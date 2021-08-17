#' Import curated Tuberculosis Data
#' A function to access available curated Tuberculosis data from the Bioconductor's
#' ExperimentHub services
#'
#' @param study_name A character or vector of characters that contain name of the studies
#' `""` will return all available studies.
#' @param dryrun Boolean. Indicate the whether downloading resources from the
#' ExperimentHub services. If `TRUE` (Default), return the names of the available resources
#' to be downloaded. If `FALSE`, start downloading data.
#' @param curated.only Boolean. Indicate whether downloading resources for the curated version.
#' If `TRUE` (Default), only download only download the curated gene expression profile and the clinical annotation information.
#' If `FALSE`, download all the available resources.
#' @return a `list` containing [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class].
#' @export
#' @examples
#' curatedTBData("GSE39939", dryrun = TRUE)
#' curatedTBData(c("GSE39939", "GSE39940"), dryrun = FALSE, curated.only = TRUE)
curatedTBData <- function(study_name, dryrun = TRUE, curated.only = TRUE) {
  # Access to experimenthub
  eh <- base::suppressWarnings(ExperimentHub::ExperimentHub())
  tbData <- AnnotationHub::query(eh, "curatedTBData")
  # List available data
  names_all <- base::unique(gsub("_.*", "", tbData$title))
  # Let study_name equals to all the the studies when "".
  if (study_name == "") {
    study_name <- names_all
  }
  indexMatch <- base::match(study_name, names_all)
  names_sub <- names_all[indexMatch]

  if (base::any(base::is.na(names_sub))) {
    if (base::all(base::is.na(names_sub))) {
      base::stop(base::sprintf(
      "The curatedTBData for the input geo accession(s):%s is/are not available. Check your input.",
                   base::paste0(study_name, collapse = ", ")), call. = FALSE)
    } else {
      indexNA <- base::which(base::is.na(names_sub))
      base::message(
        base::sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.",
                      base::paste0(study_name[indexNA], collapse = ", ")))
    }
    names_sub <- names_sub[!base::is.na(names_sub)]
  }
  if (curated.only) {
    base::message("Download curated version. Set curated.only = FALSE if want to download raw data.")
  }
  if (dryrun) {
    resources <- c()
    for (geo in names_sub) {
      resources <- c(resources, tbData$title[base::grep(geo, tbData$title)])
    }
    # Subset resources if curated.only = TRUE
    # Remove assay_raw, row_data, meta_data
    if (curated.only) {
      resources <- resources[-base::grep(base::paste0(c("assay_raw", "row_data"),
                                          collapse = "|"), resources)]
    }
    return(base::cat(base::sprintf("Attempt to download following resources for %s from the ExperimentHub service:\n%s",
                       base::paste0(names_sub, collapse = ", "),
                       base::paste0(resources, collapse = "\n"))))
  }
  df <- base::data.frame(ah_id = tbData$ah_id, title = tbData$title)
  object_list <- base::lapply(names_sub, function(x, curated.only) {
    base::message(sprintf("Now downloading: %s", x))
    # Find the index of related resources associated with the GSE names
    indexGSE <- base::grep(paste0(x, "_"), tbData$title)
    # Subset df based on the index
    dfSub <- df[indexGSE, ]
    if (curated.only) {
      indexRemove <- base::grep(base::paste0(c("assay_raw", "row_data"),
                                collapse = "|"), dfSub$title)
      dfSub <- dfSub[-indexRemove, ]
    }
    # Data download step
    data_list <- base::lapply(dfSub$ah_id, function(y)
      base::suppressMessages(tbData[[y]]))
    names(data_list) <- base::sub("^[^_]*_", "", dfSub$title)

    if (curated.only) {
      objlist <- data_list[-base::which(base::names(data_list) %in%
                                          c("column_data", "meta_data"))]
      objlist <- base::lapply(objlist, function(x) base::as.matrix(x))
    } else {
      # Create MultiAssayExperiment object using all the resources
      # Create SummarizedExperiment object from assay_raw and row_data
      exprdat <- SummarizedExperiment::SummarizedExperiment(
        base::list(assay_raw = data_list[["assay_raw"]]),
        rowData = data_list[["row_data"]])
      # Select assay_curated and/or assay_reprocess
      objlist <- data_list[-base::which(names(data_list) %in%
                                          c("column_data", "assay_raw",
                                            "row_data", "meta_data"))]
      objlist <- lapply(objlist, function(x) base::as.matrix(x))
      objlist$object_raw <- exprdat
    }
    listMap <- base::lapply(objlist, function(x) {
      base::data.frame(primary = colnames(x), colname = colnames(x),
                       stringsAsFactors = FALSE)
    })
    dfMap <- MultiAssayExperiment::listToMap(listMap)
    myMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(objlist,
                                                               data_list[["column_data"]],
                                                               dfMap)
    MultiAssayExperiment::metadata(myMultiAssay) <- base::list(data_list[["meta_data"]])
    myMultiAssay
  }, curated.only = curated.only)
  message("Finished!")
  names(object_list) <- names_sub
  return(object_list)
}
