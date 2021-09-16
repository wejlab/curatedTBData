#' @title Subset curatedTBData based on single/multiple conditions
#' @description The function selects desired samples from curatedTBData
#' database based pre-specified conditions
#' @name subset_curatedTBData
#' @param theObject A
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} or
#'   \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment} object.
#' @param annotationColName A character indicates
#'   feature of interest in the object's annotation data.
#' @param annotationCondition A vector of character indicates
#'   conditions want to be selected.
#' @param assayName A character indicates
#'   the name of the assay from the input object. The default is \code{NULL}.
#'   When \code{assayName} is \code{NULL}, the function selects
#'   the first assay along assay list.
#' @param useAssay A character indicates the name of the assay when a
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   is selected from the
#'   \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}
#'   object. The default is \code{NULL}.
#'   When \code{assayName} is \code{NULL}, the function selects
#'   the first assay along assay list.
#' @param ... Extra named arguments passed to function
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing subjects with desired annotation conditions.
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
    ## Check whether annotationColName exists in the column data
    col_info <- SummarizedExperiment::colData(theObject)
    if (!base::any(base::colnames(col_info) == annotationColName)) {
        base::paste("annotationColName:", annotationColName,
                    "is not found in the colData(theObject)",
                    "\nNULL is returned.") %>%
            base::message()
        return()
    }
    ## Check assayName whether it is specified by the users
    if (base::is.null(assayName)) {
        theObject_length <- length(SummarizedExperiment::assays(theObject))
        if (theObject_length >= 1L) {
            base::paste("assayName not specified",
                        "select the first assay as default.") %>%
                base::message()
            assay_name_exclude <- -1
        } else {
            base::stop("No available assay from the input.")
        }
    } else {
        ## Check whether assay exists in the object
        assay_names <- SummarizedExperiment::assayNames(theObject)
        assay_name_index <- base::which(assay_names %in% assayName)
        if (base::length(assay_name_index) == 0L) {
            msg1 <- base::sprintf("assay with name: %s not found",
                                  assayName)
            msg2 <- base::sprintf("\nThe available assay(s) is/are: %s",
                                  base::paste0(assay_names, collapse = ", "))
            base::paste0(msg1, msg2) %>%
                base::stop(call. = FALSE)
        }
        assay_name_exclude <- base::which(assay_names != assayName)
    }
    col_info <- SummarizedExperiment::colData(theObject)
    theObject_filter <- theObject[, col_info[, annotationColName] %in%
                                    annotationCondition]
    ## Set other assays to be NULL
    SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL
    col_info_filter <- SummarizedExperiment::colData(theObject_filter)
    annotation <- col_info_filter[, annotationColName]
    if (length(base::unique(annotation)) == length(annotationCondition)) {
        return(theObject_filter)
    } else {
        conditionNotFound <- annotationCondition[-match(unique(annotation),
                                                        annotationCondition)]
        msg <- base::sprintf("The condition %s is not found from the input.",
                             base::paste0(conditionNotFound, collapse = ", "))
        base::paste(msg, "NULL is returned,") %>%
            base::message()
        return()
    }
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData", signature = "MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL, useAssay = NULL, ...) {
    ## Check whether assay exists in the object
    if (base::is.null(assayName)) {
        if (base::length(base::names(theObject)) >= 1L) {
            base::paste("assayName not specified",
                        "select the first assay as default.") %>%
                base::message()
            assayName <- 1
        } else {
            stop("No available assay from the input.")
        }
    } else {
        experiment_name_index <- which(base::names(theObject) %in% assayName)
        if (base::length(experiment_name_index) == 0L) {
            ## Use names(theObject) when theObject is MultiAssayExperiment
            msg1 <- base::sprintf("assay with name: %s not found",
                                  assayName)
            msg2 <- base::sprintf("\nThe available assay(s) is/are: %s",
                                  paste0(names(theObject), collapse = ", "))
            base::paste0(msg1, msg2) %>%
                base::stop(call. = FALSE)
        }
    }
    theObject_sub <- theObject[[assayName]]
    col_data <- SummarizedExperiment::colData(theObject)
    if (base::ncol(theObject_sub) != base::nrow(col_data)) {
        index <- stats::na.omit(base::match(base::colnames(theObject_sub),
                                            base::row.names(col_data)))
        col_data <- col_data[index, , drop = FALSE]
    }
    if (methods::is(theObject_sub, "SummarizedExperiment")) {
        ## assay_raw is selected in the full version
        ## For those datasets that do not include all samples from the study
        SummarizedExperiment::colData(theObject_sub) <- col_data
        return(subset_curatedTBData(theObject_sub, annotationColName,
                                    annotationCondition, useAssay))
    } else if (methods::is(theObject_sub, "matrix") ||
               methods::is(theObject_sub, "data.frame")) {
        ## assay_curated/assay_reprocess is selected Set attribute to be NULL
        ## ensure that row/column names have NULL attributes
        base::colnames(theObject[[assayName]]) <-
            base::as.character(base::colnames(theObject[[assayName]]))
        base::row.names(theObject[[assayName]]) <-
            base::as.character(base::row.names(theObject[[assayName]]))
        sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
            assays = base::list(assay1 = as.matrix(theObject[[assayName]])),
            colData = col_data)
        return(subset_curatedTBData(sobject_TBSig, annotationColName,
                                    annotationCondition, "assay1"))
    } else {
        base::paste("The class of the selects assay is not recognized",
                    "\n Selected assay has class:",
                    paste0(base::class(theObject_sub), collapse = ", ")) %>%
            base::stop(call. = FALSE)
    }
})
