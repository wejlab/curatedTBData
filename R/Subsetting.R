#' Subsetting objects based on single/multiple conditions
#' @name Subset_curatedTBData
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param annotationCondition A vector indicates conditions want to be subsetted.
#' @param experiment_type A character indicates the name of the experiment within MultiAssayExperiment object.
#' Expect theObject[[experiment_type]] to be a matrix. Two special cases are:
#' When experiment_type is "all". Perform whole MultiAssayExperiment subsetting, output is in the form of MultiAsaayExperiment object.
#' When experiment_type is "assay_raw". Perform subsetting on the SummarizedExperiment Object.
#' @param ... Extra named arguments passed to function.
#' @rdname Subset_curatedTBData-methods
#' @exportMethod Subset_curatedTBData

setGeneric(name="Subset_curatedTBData", function(theObject,...){
  standardGeneric("Subset_curatedTBData")
})

#' @rdname Subset_curatedTBData-methods
setMethod("Subset_curatedTBData",
          signature="SummarizedExperiment",
          function(theObject, annotationColName, annotationCondition, UseAssay,...){

            # Check whether assay exists in the object
            assay_names <- SummarizedExperiment::assayNames(theObject_filter)
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

#' @rdname Subset_curatedTBData-methods
setMethod("Subset_curatedTBData",
          signature="MultiAssayExperiment",
          function(theObject, annotationColName, annotationCondition,
                   experiment_type){
            # Check whether experiment exists in the object
            if(experiment_type != "All"){
              experiment_name_index <- which(names(theObject) %in% experiment_type)
              if(length(experiment_name_index) == 0){
                stop(paste(experiment_type,"is not found within the object"))
              }
            }

            # For experiment_type == "all".
            # Perform whole MultiAssayExperiment selection, output is MultiAsaayExperiment
            if(experiment_type == "All"){

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
            else if (experiment_type == "assay_raw"){

                theObject_sub <- theObject[[experiment_type]]
                n <- length(annotationCondition)
                col_data <-  SummarizedExperiment::colData(theObject)

                # For those datasets that do not include all samples from the study
                if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                  index <- na.omit(match(colnames(theObject[[experiment_type]]),
                                         row.names(col_data)))
                  col_data <- col_data[index,]
                }

                SummarizedExperiment::colData(theObject_sub) <- col_data
                # subsetting annotationCondition
                sobject_TBSig_filter <- theObject_sub[,colData(theObject_sub)
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
              if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                index <- na.omit(match(colnames(theObject[[experiment_type]]),
                                       row.names(col_data)))
                col_data <- col_data[index,]
              }

              # Set atrribute to be NULL, ensure that row/column names have NULL attribute
              colnames(theObject[[experiment_type]]) <- as.character(
                                         colnames(theObject[[experiment_type]]))

              row.names(theObject[[experiment_type]]) <- as.character(
                                         row.names(theObject[[experiment_type]]))

              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(
                                  assays=list(counts = as.matrix(theObject[[experiment_type]])),
                                              colData = col_data)

              # subsetting annotationCondition
              sobject_TBSig_filter <- sobject_TBSig[,colData(sobject_TBSig)
                                                    [,annotationColName] %in% annotationCondition]
              TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
              # check if both conditions are in the column data
              if(length(unique(TB_status)) == n){
                return(sobject_TBSig_filter)
              }

            }

        }

)

#' Combine samples with common genes from selected studies, usually run after `MatchProbe`
#' @name CombineObjects
#' @param object_list A list contains expression data with probes mapped to gene symbol.
#' @param list_name A character/vector contains object name to be selected if want to combine sub-list.
#' @param experiment_type A character/vector to choose the name of the experiment from MultiAssayExperiment Object.
#' @return A SummarizedExperiment Object contains combined data from several objects.
#' @examples
#' list_name <-  c("GSE101705","GSE107104","GSE54992","GSE19444")
#' experiment_type <- "assay_MatchProbe"
#' data_list <-  curatedTBData::get_curatedTBData(list_name)
#'
#' object_norm <- lapply(data_list, function(x)
#'                                  curatedTBData::Normalization(x,
#'                                  microarray_method = "quantile",
#'                                  RNAseq_method = "TMM",
#'                                  experiment_type = "assay_raw")
#' object_match <- lapply(object_norm, function(x)
#'                                  curatedTBData::MatchProbe(x,
#'                                  UseAssay = c("TMM","quantile","RMA"),
#'                                  experiment_type = experiment_type))
#' CombineObjects(object_match, list_name, experiment_type= experiment_type)
#' @export
CombineObjects <- function(object_list,list_name=NULL,
                           experiment_type){

  if(is.null(list_name)){
    list_name <-  names(object_list)
    dat_exprs_match <- lapply(object_list, function(x) experiments(x)[[experiment_type]]
                              %>% data.frame)
  }
  else {
    dat_exprs_match <- lapply(object_list[list_name], function(x) experiments(x)[[experiment_type]]
                              %>% data.frame)
  }

  # Combine sample with common genes from selected objects.
  # Input data type should be data.frame
  dat_exprs_combine <- Reduce(
    function(x, y) merge(x, y, by = "id", all = FALSE),
    lapply(dat_exprs_match, function(x) { x$id <- rownames(x); x }))
  row.names(dat_exprs_combine) <- dat_exprs_combine$id
  dat_exprs_count <- dat_exprs_combine[,-1]
  # Create combined column data information
  Sample1 <- lapply(object_list[list_name], function(x) MultiAssayExperiment::colData(x)
                    %>% row.names())
  Sample <- unlist(Sample1, use.names=FALSE)
  col_data <- lapply(1:length(object_list[list_name]), function(x) {
    col_data <- MultiAssayExperiment::colData(object_list[list_name][[x]])
    col_data$GSE <-names(object_list[list_name][x])
    col_data
  })

  # Combine list into dataframe with unequal columns
  col_info <- plyr::rbind.fill(lapply(col_data,function(x){as.data.frame(x)}))
  row.names(col_info) <- Sample

  # Remove samples that does not exist in the count
  index <- na.omit(match(colnames(dat_exprs_count), Sample))
  col_info <- col_info[index,]

  # Create output in the format of SummarizedExperiment
  result <- SummarizedExperiment::SummarizedExperiment(assays = list(as.matrix(dat_exprs_count)),
                                                       colData = col_info)
  return(result)

}
