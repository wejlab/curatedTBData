#' @title Subsetting curatedTBData based on single/multiple conditions
#' @description \code{subset_curatedTBData} selects desired samples from curatedTBData
#' database according to pre-specified conditions.
#' @name subset_curatedTBData
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param annotationCondition A vector indicates conditions want to be subsetted.
#' @param UseAssay A character indicates the name of the assay (expression matrix) within the object.
#' Need this argument when the input is a SummarizedExperiment object.
#' @param experiment_name A character indicates the name of the experiment within MultiAssayExperiment object.
#' Expect \code{theObject[[experiment_name]]} to be a matrix. Two special cases are:
#' When experiment_name is "all". Perform whole MultiAssayExperiment subsetting, output is in the form of MultiAsaayExperiment object.
#' When experiment_name is "assay_raw". Perform subsetting on the SummarizedExperiment Object.
#' @param ... Extra named arguments passed to function.
#' @rdname subset_curatedTBData-methods
#' @exportMethod subset_curatedTBData

setGeneric(name="subset_curatedTBData", function(theObject,...){
  standardGeneric("subset_curatedTBData")
})

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature="SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition, UseAssay,...){

            # Check whether assay exists in the object
            assay_names <- SummarizedExperiment::assayNames(theObject)
            assay_name_index <- which(assay_names %in% UseAssay)
            assay_name_exclude <- which(assay_names != UseAssay)

            if(length(assay_name_index)==0){
              stop(paste(UseAssay,"is/are not found within the object"))
            }

            n <- length(annotationCondition)

            theObject_filter <- theObject[,SummarizedExperiment::colData(theObject)
                                          [,annotationColName] %in% annotationCondition]

            # Omit Assays that are not selected in the object
            SummarizedExperiment::assays(theObject_filter)[assay_name_exclude] <- NULL

            annotation <- SummarizedExperiment::colData(theObject_filter)[,annotationColName]
            if(length(unique(annotation)) == n){
              return(theObject_filter)
            }

          }
)

#' @rdname subset_curatedTBData-methods
setMethod("subset_curatedTBData",
          signature="MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition, UseAssay=NULL,
                   experiment_name,...){
            # Check whether experiment exists in the object
            if(experiment_name != "All"){
              experiment_name_index <- which(names(theObject) %in% experiment_name)
              if(length(experiment_name_index) == 0){
                stop(paste(experiment_name,"is not found within the object"))
              }
            }

            # For experiment_name == "all".
            # Perform whole MultiAssayExperiment selection, output is MultiAsaayExperiment
            if(experiment_name == "All"){

              n <- length(annotationCondition)
              theObject_filter <- theObject[,SummarizedExperiment::colData(theObject)
                                            [,annotationColName] %in% annotationCondition]
              TB_status <- SummarizedExperiment::colData(theObject_filter)[,annotationColName]
              if(length(unique(TB_status)) == n){
                return(theObject_filter)

              }
            }
              # Perform individual selection, assay_raw is SummarizedExperiment
              # output is reduced SummarizedExperiment
            else if (experiment_name == "assay_raw"){

                theObject_sub <- theObject[[experiment_name]]
                n <- length(annotationCondition)
                col_data <-  SummarizedExperiment::colData(theObject)

                # For those datasets that do not include all samples from the study
                if (ncol(theObject[[experiment_name]]) != nrow(col_data)){
                  index <- stats::na.omit(match(colnames(theObject[[experiment_name]]),
                                         row.names(col_data)))
                  col_data <- col_data[index,]
                }

                SummarizedExperiment::colData(theObject_sub) <- col_data
                # subsetting annotationCondition
                sobject_TBSig_filter <- theObject_sub[,SummarizedExperiment::colData(theObject_sub)
                                                      [,annotationColName] %in% annotationCondition]
                TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
                # check if both status are in the column data
                if(length(unique(TB_status)) == n){
                  return(sobject_TBSig_filter)
                }

              }
            # Perform individual selection, output is SummarizedExperiment
            # Potentially for TBSignatureProfiler
            # assay_reduce matrix
            else{

              n <- length(annotationCondition)
              col_data <-  SummarizedExperiment::colData(theObject)

              # when not all samples are included in the expression matrix
              # This is the cases with some RNA-seq studies
              if (ncol(theObject[[experiment_name]]) != nrow(col_data)){
                index <- stats::na.omit(match(colnames(theObject[[experiment_name]]),
                                       row.names(col_data)))
                col_data <- col_data[index,]
              }

              # Set atrribute to be NULL, ensure that row/column names have NULL attribute
              colnames(theObject[[experiment_name]]) <- as.character(
                                         colnames(theObject[[experiment_name]]))

              row.names(theObject[[experiment_name]]) <- as.character(
                                         row.names(theObject[[experiment_name]]))

              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
                                  assays=list(counts = as.matrix(theObject[[experiment_name]])),
                                              colData = col_data)

              # subsetting annotationCondition
              sobject_TBSig_filter <- sobject_TBSig[,SummarizedExperiment::colData(sobject_TBSig)
                                                    [,annotationColName] %in% annotationCondition]
              TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
              # check if both conditions are in the column data
              if(length(unique(TB_status)) == n){
                return(sobject_TBSig_filter)
              }

            }

        }
)

#' Combine samples with common genes from selected studies,
#' usually run after \code{\link{MatchProbe}}
#' @name CombineObjects
#' @param object_list A list contains expression data with probes mapped to gene symbol.
#' @param list_name A character/vector contains object name to be selected to merge.
#' @param experiment_name A character/vector to choose the name of the experiment from MultiAssayExperiment Object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' This argument passes to `mod` parameter in \code{\link[sva]{ComBat}}. Default is NULL.
#' @param batch.adjust A logical value indicating whether adjust for the batch effect.
#' Default is TRUE.
#' @return A SummarizedExperiment Object contains combined data from several objects.
#' @examples
#' list_name <-  c("GSE101705","GSE54992","GSE19444")
#' data_list <-  get_curatedTBData(list_name)
#' object_norm <- lapply(data_list, function(x)
#'                                Normalization(x, microarray_method = "quantile",
#'                                RNAseq_method = "TMM", experiment_name = "assay_raw"))
#' object_match <- lapply(object_norm, function(x)
#'                                MatchProbe(x, UseAssay = c("TMM","quantile","RMA"),
#'                                createExperimentName = "assay_MatchProbe"))
#' sobject <- CombineObjects(object_match, list_name,
#'                           experiment_name = "assay_MatchProbe",
#'                           annotationColName = "TBStatus",
#'                          batch.adjust = TRUE)
#' @export
CombineObjects <- function(object_list,list_name=NULL,
                           experiment_name, annotationColName = NULL,batch.adjust = TRUE){

  if(is.null(list_name)){
    list_name <-  names(object_list)
    dat_exprs_match <- lapply(object_list, function(x)
         MultiAssayExperiment::experiments(x)[[experiment_name]] %>% data.frame)
  }
  else {
    dat_exprs_match <- lapply(object_list[list_name], function(x)
         MultiAssayExperiment::experiments(x)[[experiment_name]] %>% data.frame)
  }

  # Combine sample with common genes from selected objects.
  # Input data type should be data.frame
  dat_exprs_combine <- Reduce(
    function(x, y) merge(x, y, by = "id", all = FALSE),
    lapply(dat_exprs_match, function(x) { x$id <- rownames(x); x }))
  row_names <- dat_exprs_combine$id
  dat_exprs_count <- dat_exprs_combine %>% dplyr::select(-.data$id) %>% data.frame()
  row.names(dat_exprs_count) <- row_names

  # Create combined column data information
  Sample1 <- lapply(object_list[list_name], function(x)
    SummarizedExperiment::colData(x) %>% row.names())

  Sample <- unlist(Sample1, use.names=FALSE)
  col_data <- lapply(seq_len(length(object_list[list_name])), function(x) {
    col_data <- SummarizedExperiment::colData(object_list[list_name][[x]])
    col_data$GSE <-names(object_list[list_name][x])
    col_data
  })

  # Combine list into dataframe with unequal columns
  col_info <- plyr::rbind.fill(lapply(col_data,function(x){as.data.frame(x)}))
  row.names(col_info) <- Sample

  # Remove samples that does not exist in the count
  index <- stats::na.omit(match(colnames(dat_exprs_count), Sample))
  col_info <- col_info[index,]

  if (batch.adjust){

    # Batch Correction
    if(is.null(annotationColName)){
      batch1 <- col_info$GSE
      combat_edata1 <- sva::ComBat(dat=as.matrix(dat_exprs_count), batch=batch1,
                                   mod=NULL)
      result <- SummarizedExperiment::SummarizedExperiment(assays = list(Batch_counts = as.matrix(combat_edata1)),
                                                           colData = col_info)
      return(result)
    }
    else{
      mod1 <- stats::model.matrix(~as.factor(col_info[,annotationColName]), data = col_info)
      batch1 <- col_info$GSE
      combat_edata1 <- sva::ComBat(dat=as.matrix(dat_exprs_count), batch=batch1,
                                   mod=mod1)
      result <- SummarizedExperiment::SummarizedExperiment(assays = list(Batch_counts = as.matrix(combat_edata1)),
                                                           colData = col_info)
      return(result)

    }
  }

  else{
  # Create output in the format of SummarizedExperiment
  result <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(dat_exprs_count)),
                                                       colData = col_info)
  return(result)
  }

}
