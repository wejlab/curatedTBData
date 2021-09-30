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
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing subjects with desired annotation conditions.
#' @export
#' @examples
#' obj <-  curatedTBData("GSE74092", dryrun = FALSE, curated.only = TRUE)
#' subset_curatedTBData(obj[[1]], annotationColName = "TBStatus",
#'                      annotationCondition = c("Control","PTB"))
#'
subset_curatedTBData <- function(theObject, annotationColName, annotationCondition,
                                 assayName = NULL) {
    ## Input class's class type
    check_type <- methods::is(theObject, "MultiAssayExperiment") ||
        methods::is(theObject, "MultiAssayExperiment")
    if (!check_type) {
        stop(paste("subset_curatedTBData() only support for",
        "SummarizedExperiment/MultiAssayExperiment object"),
             call. = FALSE)
    }
    ## Check whether annotationColName exists in the column data
    col_data <- SummarizedExperiment::colData(theObject)
    if (!any(colnames(col_data) == annotationColName)) {
        paste("annotationColName:", annotationColName,
              "is not found in the colData(theObject)",
              "\nNULL is returned.") |>
            message()
        return()
    }
    if (methods::is(theObject) == "SummarizedExperiment") {
        theObject_filter <- .subset_curatedTBData(theObject, annotationColName,
                                                  annotationCondition, assayName)
        return(theObject_filter)
    } else {
        if (is.null(assayName)) {
            if (length(base::names(theObject)) >= 1L) {
                paste("assayName not specified",
                      "select the first assay as default.") |>
                    base::message()
                assayName <- 1
            } else {
                stop("No available assay from the input.")
            }
        } else {
            experiment_name_index <- which(names(theObject) %in% assayName)
            if (base::length(experiment_name_index) == 0L) {
                ## Use names(theObject) when theObject is MultiAssayExperiment
                msg1 <- sprintf("assay with name: %s not found", assayName)
                msg2 <- sprintf("\nThe available assay(s) is/are: %s",
                                paste0(names(theObject), collapse = ", "))
                paste0(msg1, msg2) |>
                    stop(call. = FALSE)
            }
        }
        index_filter <- col_data[, annotationColName] %in% annotationCondition
        theObject_reduced <- theObject[, index_filter, assayName]
        col_info <- SummarizedExperiment::colData(theObject_reduced)
        theObject_sub <- theObject_reduced[[assayName]]
        if (methods::is(theObject_sub, "SummarizedExperiment")) {
            ## assay_raw is selected in the full version
            ## For those datasets that do not include all samples from the study
            SummarizedExperiment::colData(theObject_sub) <- col_info
            return(theObject_sub)
        } else if (methods::is(theObject_sub, "matrix") ||
                   methods::is(theObject_sub, "data.frame")) {
            ## assay_curated/assay_reprocess is selected. S
            ## set attribute to be NULL
            ## ensure that row/column names have NULL attributes
            colnames(theObject_sub) <- colnames(theObject_sub) |>
                as.character()
            row.names(theObject_sub) <- row.names(theObject_sub) |>
                as.character()
            sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
                assays = list(assay1 = as.matrix(theObject_sub)),
                colData = col_info)
            return(.subset_curatedTBData(sobject_TBSig, annotationColName,
                                        annotationCondition, "assay1"))
        } else {
            base::paste("The class of the selects assay is not recognized",
                        "\n Selected assay has class:",
                        paste0(base::class(theObject_sub), collapse = ", ")) %>%
                base::stop(call. = FALSE)
        }
    }
}

#' Subset curatedTBData based on single/multiple conditions
#' @name .subset_curatedTBData
#' @inheritParams subset_curatedTBData
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing subjects with desired annotation conditions.
.subset_curatedTBData <- function(theObject, annotationColName, annotationCondition,
                           assayName) {
    ## Check assayName whether it is specified by the users
    if (is.null(assayName)) {
        theObject_length <- length(SummarizedExperiment::assays(theObject))
        if (theObject_length >= 1L) {
            paste("assayName not specified",
                        "select the first assay as default.") |>
                message()
            assay_name_exclude <- -1
        } else {
            base::stop("No available assay from the input.")
        }
    } else {
        ## Check whether assay exists in the object
        assay_names <- SummarizedExperiment::assayNames(theObject)
        assay_name_index <- which(assay_names %in% assayName)
        if (length(assay_name_index) == 0L) {
            msg1 <- sprintf("assay with name: %s not found", assayName)
            msg2 <- sprintf("\nThe available assay(s) is/are: %s",
                            paste0(assay_names, collapse = ", "))
            paste0(msg1, msg2) |>
                stop(call. = FALSE)
        }
        assay_name_exclude <- which(assay_names != assayName)
    }
    col_info <- SummarizedExperiment::colData(theObject)
    theObject_filter <- theObject[, col_info[, annotationColName] %in%
                                      annotationCondition]
    ## Set other assays to be NULL
    SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL
    col_info_filter <- SummarizedExperiment::colData(theObject_filter)
    annotation <- col_info_filter[, annotationColName]
    if (length(unique(annotation)) == length(annotationCondition)) {
        return(theObject_filter)
    } else {
        conditionNotFound <- annotationCondition[-match(unique(annotation),
                                                        annotationCondition)]
        msg <- sprintf("The condition %s is not found from the input.",
                       paste0(conditionNotFound, collapse = ", "))
        paste(msg, "NULL is returned,") |>
            message()
        return()
    }
}
