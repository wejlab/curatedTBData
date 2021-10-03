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
curatedTBData <- function(study_name, dry.run = TRUE, curated.only = TRUE) {
    ## Access to ExperimentHub
    if (missing(study_name)) {
        stop("Argument \"study_name\" is missing, with no default.")
    }
    eh <- ExperimentHub::ExperimentHub()
    tbData <- AnnotationHub::query(eh, "curatedTBData")
    ## List available data from metadata
    metadata_file_path <- system.file("extdata/metadata.csv",
                                      package = "curatedTBData") |>
        suppressWarnings()
    if (!file.exists(metadata_file_path)) {
        metadata_file_path <- as.character("inst/extdata/metadata.csv")
    }
    metadata_tab <- utils::read.csv(metadata_file_path)
    metadata_title <- metadata_tab$Title
    if (!all(sort(metadata_title) == sort(tbData$title))) {
        paste("Resources in metadata are not consistent with",
              "resources on ExperimentHub.") |>
            warning()
    }
    names_all <- unique(gsub("_.*", "", metadata_title))
    ## Let study_name equals to all the the studies when ''.
    if (any(study_name == "")) {
        study_name <- names_all
    }
    names_sub <- names_all[match(study_name, names_all)]
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
                paste("is/are not available.") |>
                stop(call. = FALSE)
        }
        names_sub <- names_sub[!is.na(names_sub)]
    }
    ## Check dry.run status
    if (dry.run) {
        paste0("dry.run = TRUE, listing dataset(s) to be downloaded",
                    "\nSet dry.run = FALSE to download dataset(s).") |>
            message()
        resources <- lapply(names_sub, function(x)
            metadata_title[grep(x, metadata_title)]) |>
            unlist(use.names = FALSE)
        ## Subset resources if curated.only = TRUE.
        ## Remove assay_raw, row_data
        if (curated.only) {
            idx <- grep(paste0(c("assay_raw", "row_data"), collapse = "|"),
                        resources)
            resources <- resources[-idx]
        }
        msg <- sprintf("%s from the ExperiementHub:\n%s\n",
                       paste0(names_sub, collapse = ", "),
                       paste0(resources, collapse = "\n"))
        paste("Will download the following resources for", msg) |>
            message()
        # invisible: The results will not be printed if not assigned
        return(invisible(resources))
    }
    ## check whether downloading the curated version
    if (curated.only) {
        paste("curated.only = TRUE.", "Download curated version.",
              "\nSet curated.only = FALSE",
              "if want to download both raw and curated data.") |>
            message()
    }
    df <- data.frame(ah_id = tbData$ah_id, title = tbData$title)
    ## Add progress bar for lapply
    lapply_pb <- function(X, FUN, ...) {
        pb <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
        cur <- 0
        env <- environment()
        wrapper <- function(...) {
            i <- get("cur", envir = env)
            assign("cur", i + 1, envir = env)
            utils::setTxtProgressBar(get("pb", envir = env), i + 1)
            FUN(...)
        }
        re <- lapply(X, wrapper, ...)
        close(pb)
        re
    }
    object_list <- lapply(names_sub, function(x, curated.only) {
        message(sprintf("Downloading: %s", x))
        ## Find the index of related resources associated with the GSE names
        indexGSE <- grep(paste0(x, "_"), tbData$title)
        ## Subset df based on the index
        dfSub <- df[indexGSE, ]
        if (curated.only) {
            indexRemove <- grep(paste0(c("assay_raw", "row_data"),
                                       collapse = "|"),
                                dfSub$title)
            dfSub <- dfSub[-indexRemove, ]
        }
        ## Data download step with lapply_pb defined before
        data_list <- lapply_pb(dfSub$ah_id, function(y)
            suppressMessages(tbData[[y]]))
        names(data_list) <- sub("^[^_]*_", "", dfSub$title)
        if (curated.only) {
            objlist <- data_list[-which(names(data_list) %in%
                                            c("column_data", "meta_data"))]
            objlist <- lapply(objlist, function(x) as.matrix(x))
        } else {
            ## Create MultiAssayExperiment object using all the resources
            ## Create SummarizedExperiment object from assay_raw and row_data
            exprdat <- SummarizedExperiment::SummarizedExperiment(
                list(assay_raw = data_list[["assay_raw"]]),
                rowData = data_list[["row_data"]])
            # Select assay_curated and/or assay_reprocess
            objlist <- data_list[-which(names(data_list) %in%
                                            c("column_data", "assay_raw",
                                              "row_data", "meta_data"))]
            objlist <- lapply(objlist, function(x) as.matrix(x))
            objlist$object_raw <- exprdat
        }
        listMap <- lapply(objlist, function(x) {
            data.frame(primary = colnames(x), colname = colnames(x),
                             stringsAsFactors = FALSE)
        })
        myMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(
            objlist, data_list[["column_data"]])
        MultiAssayExperiment::metadata(myMultiAssay) <-
            list(data_list[["meta_data"]])
        myMultiAssay
    }, curated.only = curated.only)
    message("Finished!")
    names(object_list) <- names_sub
    return(object_list)
}
