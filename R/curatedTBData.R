#' Import curated Tuberculosis Data
#'
#' A function to access available curated tuberculosis transcriptomic data
#'   from the Bioconductor's ExperimentHub services
#'
#' @param study_name A character or vector of characters that contain name
#'   of the studies. When \code{any(study_name == "") == TRUE},
#'   the function will return all available studies.
#' @param dry.run Boolean. Indicate the whether downloading resources from the
#'   ExperimentHub services.
#'   If \code{TRUE} (Default), return the names of the available resources
#'   to be downloaded. If \code{FALSE}, start downloading data.
#' @param curated.only Boolean. Indicate whether downloading resources
#'   for the curated version.
#'   If \code{TRUE} (Default), only download the curated gene expression
#'   profile and the clinical annotation information
#'   If \code{FALSE}, download both raw and curated resources.
#' @return A \code{list} of
#'   \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}
#'   objects.
#' @export
#' @examples
#' curatedTBData("GSE39939", dry.run = TRUE)
#' curatedTBData(c("GSE39939", "GSE39940"), dry.run = FALSE, curated.only = TRUE)
#' @importFrom AnnotationHub query
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom rlang .data
#' @importFrom magrittr %>%
curatedTBData <- function(study_name, dry.run = TRUE, curated.only = TRUE) {
    ## Access to ExperimentHub
    if (missing(study_name)) {
        stop("Argument \"study_name\" is missing, with no default.")
    }
    eh <- ExperimentHub::ExperimentHub()
    tbData <- AnnotationHub::query(eh, "curatedTBData")
    ## List available data
    names_all <- unique(gsub("_.*", "", tbData$title))
    ## Let study_name equals to all the the studies when ''.
    if (any(study_name == "")) {
        study_name <- names_all
    }
    indexMatch <- match(study_name, names_all)
    names_sub <- names_all[indexMatch]

    if (any(is.na(names_sub))) {
        if (all(is.na(names_sub))) {
            sprintf("The input \"study_name\": %s",
                    paste0(study_name, collapse = ", "))  |>
                paste("is/are not available.") |>
                stop(call. = FALSE)
        } else {
            indexNA <- which(is.na(names_sub))
            sprintf("The input \"study_name\": %s",
                    paste0(study_name[indexNA], collapse = ", ")) |>
                base::paste("is/are not available.") |>
                base::stop(call. = FALSE)
        }
        names_sub <- names_sub[!is.na(names_sub)]
    }
    ## Check dry.run status
    if (dry.run) {
        paste("dry.run = TRUE, listing dataset(s) to be downloaded",
                    "\nSet dry.run = FALSE to download dataset(s).") |>
            message()
        resources <- NULL
        for (geo in names_sub) {
            resources <- c(resources,
                           tbData$title[grep(geo, tbData$title)])
        }
        ## Subset resources if curated.only = TRUE.
        ## Remove assay_raw, row_data, meta_data
        if (curated.only) {
            idx <- grep(paste0(c("assay_raw", "row_data"), collapse = "|"),
                        resources)
            resources <- resources[-idx]
        }
        msg <- sprintf("%s from the ExperiementHub:\n%s",
                       paste0(names_sub, collapse = ", "),
                       paste0(resources, collapse = "\n"))
        paste("Will download the following resources for", msg) |>
            message()
        # invisible: The results will not be printed if not assigned
        return(invisible(resources))
    }
    ## check whether downloading the curated version
    if (curated.only) {
        base::paste("curated.only = TRUE",
                    "Download curated version.",
                    "\nSet curated.only = FALSE",
                    "if want to download both raw and curated data.") %>%
            base::message()
    }
    df <- base::data.frame(ah_id = tbData$ah_id, title = tbData$title)
    object_list <- base::lapply(names_sub, function(x, curated.only) {
        base::message(sprintf("Now downloading: %s", x))
        ## Find the index of related resources associated with the GSE names
        indexGSE <- base::grep(paste0(x, "_"), tbData$title)
        ## Subset df based on the index
        dfSub <- df[indexGSE, ]
        if (curated.only) {
            indexRemove <- base::grep(base::paste0(c("assay_raw", "row_data"),
                                                   collapse = "|"),
                                      dfSub$title)
            dfSub <- dfSub[-indexRemove, ]
        }
        # Data download step
        data_list <- base::lapply(dfSub$ah_id, function(y)
            base::suppressMessages(tbData[[y]]))
        names(data_list) <- base::sub("^[^_]*_", "", dfSub$title)
        if (curated.only) {
            objlist <- data_list[-which(base::names(data_list) %in%
                                            c("column_data", "meta_data"))]
            objlist <- base::lapply(objlist, function(x) base::as.matrix(x))
        } else {
            ## Create MultiAssayExperiment object using all the resources
            ## Create SummarizedExperiment object from assay_raw and row_data
            exprdat <- SummarizedExperiment::SummarizedExperiment(
                base::list(assay_raw = data_list[["assay_raw"]]),
                rowData = data_list[["row_data"]])
            # Select assay_curated and/or assay_reprocess
            objlist <- data_list[-which(names(data_list) %in%
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
        myMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(
            objlist, data_list[["column_data"]], dfMap)
        MultiAssayExperiment::metadata(myMultiAssay) <-
            base::list(data_list[["meta_data"]])
        myMultiAssay
    }, curated.only = curated.only)
    message("Finished!")
    names(object_list) <- names_sub
    return(object_list)
}
