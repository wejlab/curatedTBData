#' @title Subsetting curatedTBData based on single/multiple conditions
#' @description \code{subset_curatedTBData} selects desired samples from curatedTBData
#' database based pre-specified conditions.
#' @name subset_curatedTBData
#' @param theObject A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' /[MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param annotationColName A string indicates feature of interest in the object's annotation data.
#' @param annotationCondition A vector of string indicates conditions want to be selected.
#' @param assayName A string indicates the name of the assay from the input object.
#' The default is `NULL`. When assayName is `NULL`, the function selects the first
#' assay along the assay list.
#' @param useAssay A string indicates the name of the assay when the
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

setGeneric(name = "subset_curatedTBData", function(theObject,...){
  standardGeneric("subset_curatedTBData")
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature = "SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName = NULL,...) {
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
              if(base::length(assay_name_index) == 0) {
                base::stop(base::sprintf(
                  "Assay with name: %s is not found from the input. The available assay(s) are %s.",
                   assayName, base::paste0(names(theObject), collapse = ", ")))
              }
              assay_name_exclude <- base::which(assay_names != assayName)
            }
            col_info <- SummarizedExperiment::colData(theObject)
            # Check whether annotationColName exists in the col data
            if (!base::any(base::colnames(col_info) == annotationColName )) {
              base::stop(base::sprintf(
                "annotationColName: %s is not found in the clicnial information."))
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
                   assayName = NULL, useAssay = NULL, ...){
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
              if(base::length(experiment_name_index) == 0) {
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
check_annotation <- function(theObject, annotationColName, annotationCondition){

  col_names <- colnames(SummarizedExperiment::colData(theObject))
  n <- length(annotationCondition)
  if(!is.na(match(annotationColName, col_names))){

    theObject_sub <- theObject[, SummarizedExperiment::colData(theObject)
                               [, annotationColName] %in% annotationCondition]
    result <- SummarizedExperiment::colData(theObject_sub)[, annotationColName]
    if(length(unique(result)) == n){
      return(theObject_sub)
    }
  }
}

#' Merge samples with common genes from selected studies
#' @name combineObjects
#' @param object_list A list contains [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] objects.
#' The assays contains object's assay contain expression data with probes mapped to gene symbol.
#' `names(object_list)` should NOT be `NULL`.
#' @param experiment_name A string/vector of string to choose the name of the assay from
#' [MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @return A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' object contains combined data from the input.
#' @examples
#' geo <-  c("GSE19435", "GSE19439")
#' data_list <-  curatedTBData(geo, dryrun = FALSE, curated.only = TRUE)
#' sobject <- combineObjects(data_list, experiment_name = "assay_curated")
#' @export
combineObjects <- function(object_list, experiment_name = NULL){
  # check the experiment_name argument
  if(base::is.null(experiment_name)) {
    base::stop("Missing experiment name for the MultiAssayExperiment object")
  }
  # check names of the input objects
  if (base::is.null(base::names(object_list))) {
    base::stop("names(object_list) should not be NULL. Add unique names for each object within the list.")
  }
  if (base::length(experiment_name) > 1) {
    # experiment name for the list of object is different
    if (base::length(experiment_name) == base::length(object_list)) {
      dat_exprs_match <- base::mapply(function(x, y) {
        MultiAssayExperiment::experiments(x)[[y]] %>%
          base::as.data.frame()
      }, object_list, experiment_name)
    } else {
      base::stop("The length of input list is not the same as the length of the experiment name vector.")
    }
  } else {
    dat_exprs_match <- base::lapply(object_list, function(x)
      MultiAssayExperiment::experiments(x)[[experiment_name]] %>%
        base::as.data.frame())
  }
  # Combine sample with common genes from a list of objects.
  # Input data type should be data.frame
  dat_exprs_combine <- base::Reduce(function(x, y)
    base::merge(x, y, by = "id", all = FALSE),
    base::lapply(dat_exprs_match, function(x) {
      x$id <- base::row.names(x)
      x
    })
  )
  row_names <- dat_exprs_combine$id
  dat_exprs_count <- dat_exprs_combine %>%
    dplyr::select(-.data$id) %>%
    base::as.data.frame()
  base::row.names(dat_exprs_count) <- row_names

  # Create combined column data information
  Sample1 <- base::lapply(object_list, function(x)
    SummarizedExperiment::colData(x) %>%
      base::row.names())

  Sample <- base::unlist(Sample1, use.names=FALSE)
  col_data <- base::lapply(base::seq_len(base::length(object_list)), function(x) {
    col_data <- SummarizedExperiment::colData(object_list[[x]])
    col_data$Study <- base::names(object_list[x])
    base::as.data.frame(col_data)
  })

  # Combine list into data frame with unequal columns, fill in NA
  # when columns from studies are not found
  rbindx <- function(dfs) {
    ns <- base::unique(base::unlist(base::sapply(dfs, base::colnames)))
    base::do.call(base::rbind, base::lapply(dfs, function(x) {
      for(n in ns[!ns %in% base::colnames(x)]) {x[[n]] <- NA}; x }))
  }
  col_info <- rbindx(col_data)
  base::row.names(col_info) <- Sample
  # Remove samples that does not exist in the count
  index <- stats::na.omit(base::match(base::colnames(dat_exprs_count), Sample))
  col_info <- col_info[index,]
  # Create output in the format of SummarizedExperiment
  result <- SummarizedExperiment::SummarizedExperiment(
    assays = base::list(assay1 = as.matrix(dat_exprs_count)),
    colData = col_info)
  return(result)
}
