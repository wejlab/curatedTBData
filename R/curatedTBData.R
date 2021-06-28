#' Import curatedTBData
#' A function to load available curatedTBData from the package.
#'
#' @param geoAccession A character/vector that contains geo accession information.
#' If `geoAccession == All, load all avaible studies in the package.
#' @return A MultiAssayExperiment Object correpsonding to the `geoAccession`
#' @export
#' @examples
#' curatedTBData("GSE39939")
#' curatedTBData(c("GSE39939", "GSE39940"))
#' curatedTBData("All")
curatedTBData <- function(geoAccession) {
  if (!is.na(match("All", geoAccession))) {
    objectList <- .curatedTBDataLoadData(geoAccession, loadAll = TRUE)
  } else {
    objectList <- .curatedTBDataLoadData(geoAccession, loadAll = FALSE)
  }
  return(objectList)
}

#' Load selected curatedTBData
#' A function to load selected datasets from the package.
#' @inheritParams curatedTBData
#' @param loadAll Logical. An logical value to indicate whether load all available studies.
#' @return A list of MultiAssayExperiment objects
.curatedTBDataLoadData <- function(geoAccession, loadAll) {
  fileNamesFull <- utils::data(package = "curatedTBData")[["results"]][, "Item"]
  fileNamesSelected <- fileNamesFull[grep("GSE", fileNamesFull)]
  geoAccession_all <- unique(gsub("_.*", "", fileNamesSelected))
  if (loadAll) {
    finalObjectList <- .createMultiAssayExperimentObject(geoAccession_all, fileNamesSelected)
    return(finalObjectList)
  } else {
    indexMatch <- match(geoAccession, geoAccession_all)
    geoAccession_sub <- geoAccession_all[indexMatch]
    if (any(is.na(geoAccession_sub))) {
      if(all(is.na(geoAccession_sub))) {
        stop(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.
                   Please double check your input geo accession.",
                     paste0(geoAccession, collapse = ", ")))
      } else {
        indexNA <- which(is.na(geoAccession_sub))
        message(sprintf("The curatedTBData for the input geo accession(s):%s is/are not available.
                   Please double check.",  paste0(geoAccession[indexNA], collapse = ", ")))
      }
      geoAccession_sub <- geoAccession_sub[!is.na(geoAccession_sub)]
    }
    finalObjectList <- .createMultiAssayExperimentObject(geoAccession_sub, fileNamesSelected)
    return(finalObjectList)
  }
}

#' Create MultiAssayExperiment Object
#' A function to to create MultiAssayExperiment Object given geo accession number
#' from the curatedTBData package.
#' @inheritParams curatedTBData
#' @param fileNameSelected A character for slected file path.
#' @return A list of MultiAssayExperiment objects
.createMultiAssayExperimentObject <- function(geoAccession, fileNamesSelected) {
  param <- BiocParallel::SerialParam(progressbar=TRUE)
  object_list <- BiocParallel::bplapply(geoAccession, function(x) {
    message(sprintf("Now loading: %s", x))
    indexGSE <- grep(paste0(x, "_"), fileNamesSelected)
    data_load <-  utils::data(list = fileNamesSelected[indexGSE],
                              package = "curatedTBData")
    data_list <- lapply(data_load, function(y) get(y))
    names(data_list) <- sub("^[^_]*_", "", data_load)
    # Remove data from environment
    objs <- ls(pos = ".GlobalEnv")
    rm(list = objs[which(objs %in% fileNamesSelected[indexGSE])], pos = ".GlobalEnv")

    # Create MultiAssatExperiment Object
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
  }, BPPARAM = param)
  names(object_list) <- geoAccession
  return(object_list)
}

