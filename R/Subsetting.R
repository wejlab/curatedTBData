#' Subsetting objects based on single/multiple conditions
#' Output should be in the orginal format either summarizedExperiment or MultiAlssayExperiment
#' experiment=NULL in default
#' Output is in the form of SummarizedExperiment format for TBsignatureProfier
#' through indication of which assay you are going to work with.
#' @name SubsetTBStatus
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param diseases A vector indicates conditions want to subset.
#' @param experiment_type A character indicates the name of the experiment within MultiAssayExperiment object.
#' @param ... Extra named arguments passed to function.
#' @rdname SubsetTBStatus-methods
#' @exportMethod SubsetTBStatus

setGeneric(name="SubsetTBStatus", function(theObject,...){
  standardGeneric("SubsetTBStatus")
})

#' @rdname SubsetTBStatus-methods
setMethod("SubsetTBStatus",
          signature="SummarizedExperiment",
          function(theObject,annotationColName, diseases, ...){

            # check eligibility
            if(!all(diseases %in% c("Control", "Latent", "PTB", "OD", "NA"))){
              stop(cat("Invalid disease, only support for",paste("Control", "Latent", "PTB", "OD", "NA", collapse = ","),
                       ". Disease type",
                       diseases[-which(diseases %in% c("Control", "Latent", "PTB", "OD", "NA"))], "is not recognized"
              ))
            }

            n <- length(diseases)
            theObject_filter <- theObject[,colData(theObject)[,annotationColName] %in% diseases]
            annotation <- SummarizedExperiment::colData(theObject_filter)[,annotationColName]
            if(length(unique(annotation)) == n){
              return(theObject_filter)
            }

          }
)

#' @rdname SubsetTBStatus-methods
setMethod("SubsetTBStatus",
          signature="MultiAssayExperiment",

          function(theObject, annotationColName, diseases,
                   experiment_type = c("all", "assay_raw", "assay_reprocess", "assay_reduce")){

            # check eligibility
            if(!all(diseases %in% c("Control", "Latent", "PTB", "OD","NA"))){
              stop(cat("Invalid disease, only support for",paste("Control", "Latent", "PTB", "OD","NA", collapse = ","),
                       ". Disease type",
                   diseases[-which(diseases %in% c("Control", "Latent", "PTB", "OD", "NA"))], "is not recognized"
                   ))
            }

            experiment_type <- match.arg(experiment_type)

            # Perform whole MultiAssayExperiment selection, output is MultiAsaayExperiment
            if(experiment_type == "all"){

              n <- length(diseases)
              theObject_filter <- theObject[,colData(theObject)[,annotationColName] %in% diseases]
              TB_status <- SummarizedExperiment::colData(theObject_filter)[,annotationColName]
              if(length(unique(TB_status)) == n){
                return(theObject_filter)

              }
            }

              # Perform individual selection, output is SummarizedExperiment
              # Potentially for TBSignatureProfiler
              # assay_reduce matrix
              if(experiment_type == "assay_reduce" || experiment_type == "assay_reprocess"){


              n <- length(diseases)
              col_data <-  colData(theObject)

              # when not all samples are included in the expression matrix
              # This is the cases with some RNA-seq data
              if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                index <- na.omit(match(colnames(theObject[[experiment_type]]), row.names(col_data)))
                col_data <- col_data[index,]
              }

              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=as.matrix(theObject[[experiment_type]])), colData = col_data)

              # subsetting diseases
              sobject_TBSig_filter <- sobject_TBSig[,colData(sobject_TBSig)[,annotationColName] %in% diseases]
              TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
              # check if both status are in the column data
              if(length(unique(TB_status)) == n){
                return(sobject_TBSig_filter)
              }

              }

              # Perform individual selection, assay_raw is SummarizedExperiment
              # output is reduced SummarizedExperiment
              if (experiment_type == "assay_raw"){
                theObject_sub <- theObject[[experiment_type]]
                n <- length(diseases)
                col_data <-  colData(theObject)
                if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                  index <- na.omit(match(colnames(theObject[[experiment_type]]), row.names(col_data)))
                  col_data <- col_data[index,]
                }

                colData(theObject_sub) <- col_data
                # subsetting diseases
                sobject_TBSig_filter <- theObject_sub[,colData(theObject_sub)[,annotationColName] %in% diseases]
                TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)[,annotationColName]
                # check if both status are in the column data
                if(length(unique(TB_status)) == n){
                  return(sobject_TBSig_filter)
                }

              }

            }

)


#' Combine samples with common genes from selected studies, usually run after `MatchProbe`
#' @name CombineObjects
#' @param object_list A list contains expression data with mapped gene symbol.
#' @param gse_name A character/vector (GEO accession number) contains object name want to combine.
#' @param experiment_type A character/vector to choose the name of the experiment from MultiAssayExperiment Object.
#' @return A SummarizedExperiment Object contains combined data from several objects.
#'
#' @export
CombineObjects <- function(object_list,gse_name=NULL,
                           experiment_type = c("assay_reduce","assay_reprocess_norm")){

  experiment_type <- match.arg(experiment_type)

  if(is.null(gse_name)){
    gse_name <-  names(object_list)
    dat_exprs_match <- lapply(object_list, function(x) experiments(x)[[experiment_type]] %>% data.frame)
  }
  else {
    dat_exprs_match <- lapply(object_list[gse_name], function(x) experiments(x)[[experiment_type]] %>% data.frame)
  }

  # Combine sample with common genes from selected objects.
  # Input data type should be data.frame
  dat_exprs_combine <- Reduce(
    function(x, y) merge(x, y, by = "id", all = F),
    lapply(dat_exprs_match, function(x) { x$id <- rownames(x); x }))
  row.names(dat_exprs_combine) <- dat_exprs_combine$id
  dat_exprs_count <- dat_exprs_combine[,-1]
  # Create combined column data information
  Sample1 <- lapply(object_list[gse_name], function(x) MultiAssayExperiment::colData(x) %>% row.names())
  Sample <- unlist(Sample1, use.names=FALSE)
  col_data <- lapply(1:length(object_list[gse_name]), function(x) {
    col_data <- MultiAssayExperiment::colData(object_list[gse_name][[x]])
    col_data$GSE <-names(object_list[gse_name][x])
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

#' Remove empty objects from list contains both SummariexExperiment and MultiAssayExpriment objects.
#' @name remove_empty_object
#' @param k A list contains both SummariexExperiment/MultiAssayExpriment objects.
#' @return A list contains non-empty SummariexExperiment/MultiAssayExpriment object.
#'
#' @export
remove_empty_object <- function(k){
  x <- k
  for (i in which(sapply(x, function(x) class(x) == "SummarizedExperiment"))){
    if(nrow(colData(x[[i]]))==0){
      x[[i]] <- NA
    }
  }
  for (j in which(sapply(x, function(x) class(x) == "MultiAssayExperiment"))){
    if(length(experiments(x[[j]]))==0){
      x[[j]] <- NA
    }
  }
  x <- x[!is.na(x)]
  return(x)
}
