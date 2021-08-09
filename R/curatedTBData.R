#' Import curated Tuberculosis Data
#' A function to access available curated Tuberculosis data from the Bioconductor's
#' ExperimentHub services
#'
#' @param geoAccession A string or vector of strings that contains geo accession information.
#' @param dryrun Boolean. Indicate the whether downloading resources from the
#' ExperimenttHub services. If `TRUE`, return the names of the resources.
#' If `FALSE`, start downloading data.
#' @return a `list` containing [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class]
#' @export
#' @examples
#' curatedTBData("GSE39939")
#' curatedTBData(c("GSE39939", "GSE39940"))
curatedTBData <- function(geoAccession, dryrun = FALSW) {
  eh <- suppressWarnings(ExperimentHub::ExperimentHub())
  tbData <- AnnotationHub::query(eh, "curatedTBData")
  # List available data
  geoAccession_all <- unique(gsub("_.*", "", tbData$title))
  indexMatch <- match(geoAccession, geoAccession_all)
  geoAccession_sub <- geoAccession_all[indexMatch]

  if (any(is.na(geoAccession_sub))) {
    if(all(is.na(geoAccession_sub))) {
      stop(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.
                   Check your input.",
                   paste0(geoAccession, collapse = ", ")), call. = FALSE)
    } else {
      indexNA <- which(is.na(geoAccession_sub))
      message(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.",
                      paste0(geoAccession[indexNA], collapse = ", ")))
    }
    geoAccession_sub <- geoAccession_sub[!is.na(geoAccession_sub)]
  }
  if (dryrun) {
    return(sprintf("Attempt to download associated data for %s from the ExperimentHub service",
                   paste0(geoAccession_sub, collapse = ", ")))
  }

  object_list <- lapply(geoAccession_sub, function(x) {
    message(sprintf("Now downloading: %s", x))
    indexGSE <- grep(paste0(x, "_"), tbData$title)
    data_list <- lapply(tbData$ah_id[indexGSE], function(y) suppressMessages(tbData[[y]]))
    names(data_list) <- sub("^[^_]*_", "", tbData$title[indexGSE])

    # Create MultiAssayExperiment Object
    exprdat <- SummarizedExperiment::SummarizedExperiment(data_list[["assay_raw"]],
                                                          rowData = data_list[["row_data"]])
    objlist <- data_list[-which(names(data_list) %in% c("column_data", "assay_raw",
                                                        "row_data", "meta_data"))]
    objlist <- lapply(objlist, function(x) as.matrix(x))
    objlist$object_raw <- exprdat
    listMap <- lapply(objlist, function(x) {
      data.frame(primary = colnames(x), colname = colnames(x), stringsAsFactors = FALSE)
    })
    dfMap <- MultiAssayExperiment::listToMap(listMap)
    myMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(
      objlist, data_list[["column_data"]], dfMap)
    MultiAssayExperiment::metadata(myMultiAssay) <- list(data_list[["meta_data"]])
    myMultiAssay
  })
  names(object_list) <- geoAccession_sub
  return(object_list)
}
