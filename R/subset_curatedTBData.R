#' @title Subset curatedTBData based on single/multiple conditions
#' @description \code{subset_curatedTBData} selects desired samples from curatedTBData
#' database based pre-specified conditions.
#' @name subset_curatedTBData
#' @param theObject A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' /[MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param annotationColName A character indicates feature of interest in the object's annotation data.
#' @param annotationCondition A vector of character indicates conditions want to be selected.
#' @param assayName A character indicates the name of the assay from the input object.
#' The default is `NULL`. When assayName is `NULL`, the function selects the first
#' assay along the assay list.
#' @param useAssay A character indicates the name of the assay when the
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class] is
#' selected from the [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] object.
#' The default is `NULL`. When assayName is `NULL`, the function selects the first
#' assay along the assay list.
#' @param ... Extra named arguments passed to function.
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
setMethod("subset_curatedTBData",
          signature = "SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL, ...) {
            # Run some diagnostics
            # check assayName whether it is specified by the users
            if (base::is.null(assayName)) {
              if (base::length(SummarizedExperiment::assays(theObject)) >= 1) {
                base::message("assayName not specified, select the first assay as default.")
                assay_name_exclude <- -1
              } else {
                base::stop("No available assay from the input.")
              }
            } else {
              # Check whether assay exists in the object
              assay_names <- SummarizedExperiment::assayNames(theObject)
              assay_name_index <- base::which(assay_names %in% assayName)
              if (base::length(assay_name_index) == 0) {
                base::stop(base::sprintf(
                  "Assay with name: %s is not found from the input. The available assay(s) are %s.",
                   assayName, base::paste0(names(theObject), collapse = ", ")))
              }
              assay_name_exclude <- base::which(assay_names != assayName)
            }
            col_info <- SummarizedExperiment::colData(theObject)
            # Check whether annotationColName exists in the col data
            if (!base::any(base::colnames(col_info) == annotationColName)) {
              base::message(base::sprintf(
                "annotationColName: %s is not found in the clicnial information. NULL is returned.",
                annotationColName))
              return()
            }
            theObject_filter <- theObject[, col_info[, annotationColName] %in% annotationCondition]
            # Set other assays to be NULL
            SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL

            annotation <- SummarizedExperiment::colData(theObject_filter)[, annotationColName]
            if (base::length(base::unique(annotation)) ==
                base::length(annotationCondition)) {
              return(theObject_filter)
            } else {
              base::message(base::sprintf(
                "The condition %s is not found from the input data, NULL is returned.",
                 base::paste0(annotationCondition[-match(unique(annotation), annotationCondition)],
                              collapse = ", ")))
            }
          })

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature = "MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL, useAssay = NULL, ...) {
            # Check whether assay exists in the object
            if (base::is.null(assayName)) {
              if (base::length(base::names(theObject)) >= 1) {
                base::message("assayName not specified, select the first assay as default.")
                assayName <- 1
              } else {
                stop("No available assay from the input.")
              }
            } else {
              experiment_name_index <- base::which(base::names(theObject) %in%
                                                     assayName)
              if (base::length(experiment_name_index) == 0) {
                base::stop(base::sprintf("Assay with name: %s is not found from the input.
                                   The available assay(s) are %s.", assayName,
                                         base::paste0(base::names(theObject),
                                                      collapse = ", ")))
              }
            }
            theObject_sub <- theObject[[assayName]]
            col_data <-  SummarizedExperiment::colData(theObject)
            if (base::ncol(theObject_sub) != base::nrow(col_data)) {
              index <- stats::na.omit(base::match(base::colnames(theObject_sub),
                                                  base::row.names(col_data)))
              col_data <- col_data[index, ]
            }
            if (base::class(theObject_sub)[1] == "SummarizedExperiment") {
              # assay_raw is selected in the full version
              # For those datasets that do not include all samples from the study
              SummarizedExperiment::colData(theObject_sub) <- col_data
              return(subset_curatedTBData(theObject_sub, annotationColName,
                                          annotationCondition, useAssay))
            } else if (base::class(theObject_sub)[1] == "matrix" ||
                       base::class(theObject_sub)[1] == "data.frame") {
              # assay_curated/assay_reprocess is selected
              # Set attribute to be NULL, ensure that row/column names have NULL attributes
              base::colnames(theObject[[assayName]]) <- base::as.character(
                base::colnames(theObject[[assayName]]))
              base::row.names(theObject[[assayName]]) <- base::as.character(
                base::row.names(theObject[[assayName]]))

              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
                assays = base::list(assay1 = base::as.matrix(theObject[[assayName]])),
                colData = col_data)
              return(subset_curatedTBData(sobject_TBSig, annotationColName,
                                          annotationCondition, "assay1"))
            } else {
              base::stop(base::sprintf("The class of the selected assay class is not recognized.
                           Selected assay has class: %s"),
                         paste0(base::class(theObject_sub), collapse = ", "))
            }
          })

#' Check the annotation column name in the colData function
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param annotationCondition A vector indicates conditions want to be chosen.
#' @export
check_annotation <- function(theObject, annotationColName, annotationCondition) {

  col_names <- colnames(SummarizedExperiment::colData(theObject))
  n <- length(annotationCondition)
  if (!is.na(match(annotationColName, col_names))) {

    theObject_sub <- theObject[, SummarizedExperiment::colData(theObject)
                               [, annotationColName] %in% annotationCondition]
    result <- SummarizedExperiment::colData(theObject_sub)[, annotationColName]
    if (length(unique(result)) == n) {
      return(theObject_sub)
    }
  }
}

