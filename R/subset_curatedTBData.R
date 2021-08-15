#' @title Subsetting curatedTBData based on single/multiple conditions
#' @description \code{subset_curatedTBData} selects desired samples from curatedTBData
#' database according to pre-specified conditions.
#' @name subset_curatedTBData
#' @param theObject A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' /[MultiAssayExperiment][MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param annotationCondition A vector indicates conditions want to be selected.
#' @param assayName A character indicates the name of the assay from the input object.
#' @param ... Extra named arguments passed to function.
#' @rdname subset_curatedTBData-methods
#' @exportMethod subset_curatedTBData

setGeneric(name = "subset_curatedTBData", function(theObject,...){
  standardGeneric("subset_curatedTBData")
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature = "SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition,
                   assayName,...){
            if (missing(assayName)) {
              stop("Missing assayName value, please specify the assay name of the input object.")
            }
            # Check whether assay exists in the object
            assay_names <- SummarizedExperiment::assayNames(theObject)
            assay_name_index <- which(assay_names %in% assayName)
            assay_name_exclude <- which(assay_names != assayName)
            if(length(assay_name_index) == 0) {
              stop(sprintf("%s is/are not found from the object. Available assay(s) is/are %s",
                           assayName, paste0(assay_names, collapse = ", ")))
            }

            theObject_filter <- theObject[, SummarizedExperiment::colData(theObject)
                                          [, annotationColName] %in% annotationCondition]
            SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL

            annotation <- SummarizedExperiment::colData(theObject_filter)[, annotationColName]
            if (length(unique(annotation)) == length(annotationCondition)) {
              return(theObject_filter)
            } else {
              message(sprintf("The condition %s is not found from the input data, NULL is returned.",
                              annotationCondition[-match(unique(annotation), annotationCondition)]))
            }
          }
)

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature = "MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition, assayName,...){
            # Check whether assay exists in the object
            if (!is.null(assayName)) {
              experiment_name_index <- which(names(theObject) %in% assayName)
              if(length(experiment_name_index) == 0){
                stop(sprintf("%s is/are not found within the object", assayName))
              }
            }

            # For assayName == "all".
            # Perform whole MultiAssayExperiment selection, output is MultiAsaayExperiment
            if(assayName == "All"){

              n <- length(annotationCondition)
              theObject_filter <- theObject[,SummarizedExperiment::colData(theObject)
                                            [,annotationColName] %in% annotationCondition]
              result <- SummarizedExperiment::colData(theObject_filter)[,annotationColName]
              if(length(unique(result)) == n){
                return(theObject_filter)

              }
            } else if (assayName == "assay_raw"){
              # Perform individual selection, assay_raw is SummarizedExperiment
              # output is reduced SummarizedExperiment

              theObject_sub <- theObject[[assayName]]
              n <- length(annotationCondition)
              col_data <-  SummarizedExperiment::colData(theObject)

              # For those datasets that do not include all samples from the study
              if (ncol(theObject[[assayName]]) != nrow(col_data)){
                index <- stats::na.omit(match(colnames(theObject[[experiment_name]]),
                                         row.names(col_data)))
                col_data <- col_data[index,]
              }

              SummarizedExperiment::colData(theObject_sub) <- col_data
              # subsetting annotationCondition
              sobject_TBSig_filter <- theObject_sub[,SummarizedExperiment::colData(theObject_sub)
                                                    [,annotationColName] %in% annotationCondition]
              result <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
                # check if both status are in the column data
              if(length(unique(result)) == n){
                return(sobject_TBSig_filter)
              }

            } else {
              # Perform individual selection, output is SummarizedExperiment
              # Potentially for TBSignatureProfiler
              # assay_reduce matrix
              col_data <-  SummarizedExperiment::colData(theObject)

              # when not all samples are included in the expression matrix
              # This is the cases with some RNA-seq studies
              if (ncol(theObject[[assayName]]) != nrow(col_data)){
                index <- stats::na.omit(match(colnames(theObject[[assayName]]),
                                       row.names(col_data)))
                col_data <- col_data[index,]
              }

              # Set attribute to be NULL, ensure that row/column names have NULL attributes
              colnames(theObject[[assayName]]) <- as.character(
                                         colnames(theObject[[assayName]]))

              row.names(theObject[[assayName]]) <- as.character(
                                         row.names(theObject[[assayName]]))

              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
                                  assays=list(counts = as.matrix(theObject[[assayName]])),
                                              colData = col_data)

              # subsetting annotationCondition
              sobject_TBSig_filter <- sobject_TBSig[,SummarizedExperiment::colData(sobject_TBSig)
                                                    [,annotationColName] %in% annotationCondition]
              result <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
              # check if both conditions are in the column data
              if(length(unique(result)) == length(annotationCondition)){
                return(sobject_TBSig_filter)
              } else {
                message(sprintf("The condition %s is not found from the input data, NULL is returned.",
                                annotationCondition[-match(unique(result), annotationCondition)]))
              }
            }
        }
)


#' Check the annotation column name in the colData function
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param annotationCondition A vector indicates conditions want to be subsetted.
#'
#' @export
check_annotation <- function(theObject, annotationColName, annotationCondition){

  col_names <- colnames(SummarizedExperiment::colData(theObject))
  n <- length(annotationCondition)
  if(!is.na(match(annotationColName, col_names))){

    theObject_sub <- theObject[, SummarizedExperiment::colData(theObject)
                               [,annotationColName] %in% annotationCondition]
    result <- SummarizedExperiment::colData(theObject_sub)[,annotationColName]
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
  if(is.null(experiment_name)) {
    stop("Missing experiment name for the MultiAssayExperiment object")
  }
  # check names of the input objects
  if (is.null(names(object_list))) {
    stop("names(object_list) should not be NULL. Add unique names for each object within the list.")
  }
  if (length(experiment_name) > 1) {
    # experiment name for the list of object is different
    if (length(experiment_name) == length(object_list)) {
      dat_exprs_match <- mapply(function(x, y) {
        MultiAssayExperiment::experiments(x)[[y]] %>% as.data.frame
      }, object_list, experiment_name)
    } else {
      stop("The length of input list is not the same as the length of the experiment name vector.")
    }

  } else {
    dat_exprs_match <- lapply(object_list, function(x)
      MultiAssayExperiment::experiments(x)[[experiment_name]] %>% as.data.frame)
  }


  # Combine sample with common genes from a list of objects.
  # Input data type should be data.frame
  dat_exprs_combine <- Reduce(function(x, y) merge(x, y, by = "id", all = FALSE),
    lapply(dat_exprs_match, function(x) {
      x$id <- rownames(x)
      x}))
  row_names <- dat_exprs_combine$id
  dat_exprs_count <- dat_exprs_combine %>% dplyr::select(-.data$id) %>% as.data.frame()
  row.names(dat_exprs_count) <- row_names

  # Create combined column data information
  Sample1 <- lapply(object_list, function(x)
    SummarizedExperiment::colData(x) %>% row.names())

  Sample <- unlist(Sample1, use.names=FALSE)
  col_data <- lapply(seq_len(length(object_list)), function(x) {
    col_data <- SummarizedExperiment::colData(object_list[[x]])
    col_data$Study <-names(object_list[x])
    col_data
  })

  # Combine list into data frame with unequal columns
  col_info <- plyr::rbind.fill(lapply(col_data, function(x) as.data.frame(x) ))
  row.names(col_info) <- Sample

  # Remove samples that does not exist in the count
  index <- stats::na.omit(match(colnames(dat_exprs_count), Sample))
  col_info <- col_info[index,]

  # Create output in the format of SummarizedExperiment
  result <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(dat_exprs_count)),
                                                       colData = col_info)
  return(result)
}
