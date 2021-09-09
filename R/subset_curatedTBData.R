#' @title Subset curatedTBData based on single/multiple conditions
#' @description The function selects desired samples from curatedTBData
#' database based pre-specified conditions
#' @name subset_curatedTBData
#' @param theObject A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} or
#' \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment} object
#' @param annotationColName A character indicates feature of interest in the object's annotation data
#' @param annotationCondition A vector of character indicates conditions want to be selected
#' @param assayName A character indicates the name of the assay from the input object.
#' The default is \code{NULL}. When \code{assayName} is \code{NULL}, the function selects the first
#' assay along list of assays
#' @param useAssay A character indicates the name of the assay when a
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} is selected from the \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment} object.
#' The default is \code{NULL}. When \code{assayName} is \code{NULL}, the function selects the first
#' assay along list of assays
#' @param ... Extra named arguments passed to function
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' object containing subjects with desired annotation conditions
#' @examples
#' obj <-  curatedTBData("GSE74092", dryrun = FALSE, curated.only = TRUE)
#' subset_curatedTBData(obj[[1]], annotationColName = "TBStatus",
#'                      annotationCondition = c("Control","PTB"))
#' @rdname subset_curatedTBData-methods
#' @exportMethod subset_curatedTBData

setGeneric(name = "subset_curatedTBData", function(theObject, ...) {
  standardGeneric("subset_curatedTBData")
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData", signature = "SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL, ...) {
    # Run some diagnostics check assayName whether it is specified by the users
    if (base::is.null(assayName)) {
        if (base::length(SummarizedExperiment::assays(theObject)) >= 1L) {
            base::message("assayName not specified, select the first assay as default.")
            assay_name_exclude <- -1
        } else {
            base::stop("No available assay from the input.")
        }
    } else {
        # Check whether assay exists in the object
        assay_names <- SummarizedExperiment::assayNames(theObject)
        assay_name_index <- base::which(assay_names %in% assayName)
        if (base::length(assay_name_index) == 0L) {
            base::stop(base::sprintf("Assay with name: %s is not found from the input. The available assay(s) are %s.",
                                     assayName, base::paste0(assay_names, collapse = ", ")),
                       call. = FALSE)
        }
        assay_name_exclude <- base::which(assay_names != assayName)
    }
    col_info <- SummarizedExperiment::colData(theObject)
    # Check whether annotationColName exists in the col data
    if (!base::any(base::colnames(col_info) == annotationColName)) {
        base::message(base::sprintf("annotationColName: %s is not found in the clicnial information. NULL is returned.", annotationColName))
        return()
    }
    theObject_filter <- theObject[, col_info[, annotationColName] %in% annotationCondition]
    # Set other assays to be NULL
    SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL

    annotation <- SummarizedExperiment::colData(theObject_filter)[, annotationColName]
    if (base::length(base::unique(annotation)) == base::length(annotationCondition)) {
        return(theObject_filter)
    } else {
        conditionNotFound <- annotationCondition[-match(unique(annotation),
                                                        annotationCondition)]
        base::message(base::sprintf("The condition %s is not found from the input data, NULL is returned.",
                                    base::paste0(conditionNotFound, collapse = ", ")))
    }
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData", signature = "MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL, useAssay = NULL, ...) {
    ## Check whether assay exists in the object
    if (base::is.null(assayName)) {
        if (base::length(base::names(theObject)) >= 1L) {
            base::message("assayName not specified, select the first assay as default.")
            assayName <- 1
        } else {
            stop("No available assay from the input.")
        }
    } else {
        experiment_name_index <- base::which(base::names(theObject) %in% assayName)
        if (base::length(experiment_name_index) == 0) {
            base::stop(base::sprintf("Assay with name: %s is not found from the input.
                                   The available assay(s) are %s.",
                                     assayName, base::paste0(base::names(theObject), collapse = ", ")))
        }
    }
    theObject_sub <- theObject[[assayName]]
    col_data <- SummarizedExperiment::colData(theObject)
    if (base::ncol(theObject_sub) != base::nrow(col_data)) {
        index <- stats::na.omit(base::match(base::colnames(theObject_sub),
                                            base::row.names(col_data)))
        col_data <- col_data[index, ]
    }
    if (methods::is(theObject_sub, "SummarizedExperiment")) {
        ## assay_raw is selected in the full version For those datasets that do not include all samples from the study
        SummarizedExperiment::colData(theObject_sub) <- col_data
        return(subset_curatedTBData(theObject_sub, annotationColName,
                                    annotationCondition, useAssay))
    } else if (methods::is(theObject_sub, "matrix") || methods::is(theObject_sub, "data.frame")) {
        ## assay_curated/assay_reprocess is selected Set attribute to be NULL, ensure that row/column names have NULL attributes
        base::colnames(theObject[[assayName]]) <- base::as.character(base::colnames(theObject[[assayName]]))
        base::row.names(theObject[[assayName]]) <- base::as.character(base::row.names(theObject[[assayName]]))

        sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays = base::list(assay1 = base::as.matrix(theObject[[assayName]])),
                                                                    colData = col_data)
        return(subset_curatedTBData(sobject_TBSig, annotationColName, annotationCondition, "assay1"))
    } else {
        base::stop(base::sprintf("The class of the selected assay class is not recognized.
                           Selected assay has class: %s"),
                   paste0(base::class(theObject_sub), collapse = ", "))
    }
})
