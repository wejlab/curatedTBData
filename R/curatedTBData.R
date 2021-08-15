#' Import curated Tuberculosis Data
#' A function to access available curated Tuberculosis data from the Bioconductor's
#' ExperimentHub services
#'
#' @param study_name A string or vector of strings that contain name of the datasets.
#' @param dryrun Boolean. Indicate the whether downloading resources from the
#' ExperimentHub services. If `TRUE` (Default), return the names of the available resources
#' to be downloaded. If `FALSE`, start downloading data.
#' @param curated.only Boolean. Indicate whether downloading resources for the curated version.
#' If `TRUE` (Default), only download only download the curated gene expression profile and the clinical annotation information.
#' If `FALSE`, download all the available resources.
#' @return a `list` containing [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class].
#' @export
#' @examples
#' curatedTBData("GSE39939")
#' curatedTBData(c("GSE39939", "GSE39940"))
curatedTBData <- function(study_name, dryrun = TRUE, curated.only = TRUE) {

  eh <- suppressWarnings(ExperimentHub::ExperimentHub())
  tbData <- AnnotationHub::query(eh, "curatedTBData")
  # List available data
  names_all <- unique(gsub("_.*", "", tbData$title))
  indexMatch <- match(study_name, names_all)
  names_sub <- names_all[indexMatch]

  if (any(is.na(names_sub))) {
    if (all(is.na(names_sub))) {
      stop(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.
                   Check your input.",
                   paste0(study_name, collapse = ", ")), call. = FALSE)
    } else {
      indexNA <- which(is.na(names_sub))
      message(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.",
                      paste0(study_name[indexNA], collapse = ", ")))
    }
    names_sub <- names_sub[!is.na(names_sub)]
  }
  if (curated.only) {
    message("Download curated version. Set curated.only = FALSE if want to download raw data.")
  }
  if (dryrun) {
    resources <- c()
    for (geo in names_sub) {
      resources <- c(resources, tbData$title[grep(geo, tbData$title)])
    }
    # Subset resources if curated.only = TRUE
    # Remove assay_raw, row_data, meta_data
    if (curated.only) {
      resources <- resources[-grep(paste0(c("assay_raw", "row_data"),
                                          collapse = "|"), resources)]
    }
    return(cat(sprintf("Attempt to download following resources for %s from the ExperimentHub service:\n%s",
                       paste0(names_sub, collapse = ", "),
                       paste0(resources, collapse = "\n"))))
  }
  df <- data.frame(ah_id = tbData$ah_id, title = tbData$title)
  object_list <- lapply(names_sub, function(x, curated.only) {
    message(sprintf("Now downloading: %s", x))
    # Find the index of related resources associated with the GSE names
    indexGSE <- grep(paste0(x, "_"), tbData$title)
    # Subset df based on the index
    dfSub <- df[indexGSE, ]
    if (curated.only) {
      indexRemove <- grep(paste0(c("assay_raw", "row_data"),
                                collapse = "|"), dfSub$title)
      dfSub <- dfSub[-indexRemove, ]
    }
    # Data download step
    data_list <- lapply(dfSub$ah_id, function(y) suppressMessages(tbData[[y]]))
    names(data_list) <- sub("^[^_]*_", "", dfSub$title)

    if (curated.only) {
      objlist <- data_list[-which(names(data_list) %in% c("column_data", "meta_data"))]
      objlist <- lapply(objlist, function(x) as.matrix(x))
    } else {
      # Create MultiAssayExperiment Object using all the resources
      # Create SummarizedExperiment object from assay_raw and row_data
      exprdat <- SummarizedExperiment::SummarizedExperiment(data_list[["assay_raw"]],
                                                            rowData = data_list[["row_data"]])
      # Filtered out assay_curated and/or assay_reprocess
      objlist <- data_list[-which(names(data_list) %in% c("column_data", "assay_raw",
                                                          "row_data", "meta_data"))]
      objlist <- lapply(objlist, function(x) as.matrix(x))
      objlist$object_raw <- exprdat
    }
    listMap <- lapply(objlist, function(x) {
      data.frame(primary = colnames(x), colname = colnames(x), stringsAsFactors = FALSE)
    })
    dfMap <- MultiAssayExperiment::listToMap(listMap)
    myMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(objlist,
                                                               data_list[["column_data"]],
                                                               dfMap)
    MultiAssayExperiment::metadata(myMultiAssay) <- list(data_list[["meta_data"]])
    myMultiAssay
  }, curated.only = curated.only)
  message("Finished!")
  names(object_list) <- names_sub
  return(object_list)
}
