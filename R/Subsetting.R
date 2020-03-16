#' Subsetting objects based on single/multiple conditions
setGeneric(name="SubsetSample", function(theObject,...){
  standardGeneric("SubsetSample")
})

setMethod("SubsetSample",
          signature="SummarizedExperiment",
          function(theObject,diseases){
            n <- length(diseases)
            theObject_filter <- theObject[,theObject$TBStatus %in% diseases]
            TB_status <- SummarizedExperiment::colData(theObject_filter)["TBStatus"][,1]
            if(length(unique(TB_status)) == n){
              return(theObject_filter)
            }

    }
)

setMethod("SubsetSample",
          signature="MultiAssayExperiment",

          function(theObject,diseases, experiment_type = NULL){
            if(!is.null(experiment_type)){
              n <- length(diseases)
              #col_info <- colData(theObject)
              #col_data <- data.frame(Sample=row.names(col_info) %>% as.factor(),
              #                       Disease = col_info$TBStatus %>% as.factor())
              #row.names(col_data) <- row.names(col_info)
              col_data <-  colData(theObject)

              # when not all samples are included in the expression matrix
              # This is the cases with some RNA-seq data
              if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                index <- sapply(1:length(colnames(theObject[[experiment_type]])), function (i)
                  which(row.names(col_data) %in% colnames(theObject[[experiment_type]])[i]))

                col_data <- col_data[index,]
              }


              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(theObject[[experiment_type]])), colData = col_data)

              # subsetting disease1 and disease2
              sobject_TBSig_filter <- sobject_TBSig[,sobject_TBSig$TBStatus %in% diseases]
              TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)["TBStatus"][,1]
              # check if both status are in the column data
              if(length(unique(TB_status)) == n){
                return(sobject_TBSig_filter)
              }

            }
            n <- length(diseases)
            theObject_filter <- theObject[,theObject$TBStatus %in% diseases]
            TB_status <- SummarizedExperiment::colData(theObject_filter)["TBStatus"][,1]
            if(length(unique(TB_status)) == n){
              return(theObject_filter)
            }

          }
)

#' Remove objects based on single/multiple conditions
setGeneric(name="RemoveSample", function(theObject,...){
  standardGeneric("RemoveSample")
})

setMethod("RemoveSample",
          signature="SummarizedExperiment",
          function(theObject,ColName, Con){

            n <- length(Con)

            theObject_filter <- theObject[,theObject[,colName] != Con]
            result <- SummarizedExperiment::colData(theObject_filter)[ColName][,1]
            if(length(unique(Con)) == n){
              return(theObject_filter)
            }

          }
)

setMethod("SubsetSample",
          signature="MultiAssayExperiment",

          function(theObject,ColName, Con, experiment_type = NULL){
            if(!is.null(experiment_type)){
              n <- length(Con)
              #col_info <- colData(theObject)
              #col_data <- data.frame(Sample=row.names(col_info) %>% as.factor(),
              #                       Disease = col_info$TBStatus %>% as.factor())
              #row.names(col_data) <- row.names(col_info)
              col_data <-  colData(theObject)

              # when not all samples are included in the expression matrix
              # This is the cases with some RNA-seq data
              if (ncol(theObject[[experiment_type]]) != nrow(col_data)){
                index <- sapply(1:length(colnames(theObject[[experiment_type]])), function (i)
                  which(row.names(col_data) %in% colnames(theObject[[experiment_type]])[i]))

                col_data <- col_data[index,]
              }


              sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(theObject[[experiment_type]])), colData = col_data)

              # subsetting disease1 and disease2
              sobject_TBSig_filter <- sobject_TBSig[,sobject_TBSig[,ColName] != Con]
              result <- SummarizedExperiment::colData(sobject_TBSig_filter)[ColName][,1]
              # check if both status are in the column data
              if(length(unique(result)) == n){
                return(sobject_TBSig_filter)
              }

            }
            n <- length(Con)
            theObject_filter <- theObject[,theObject[,ColName] != Con]
            result <- SummarizedExperiment::colData(theObject_filter)[ColName][,1]
            if(length(unique(result)) == n){
              return(theObject_filter)
            }

          }
)

#' Remove empty objects from list contains both SummariexExperiment and MultiAssayExpriment objects
#' @name remove_empty_object
#' @param k A list contains both SummariexExperiment/MultiAssayExpriment objects
#' @return A list contains non-empty SummariexExperiment/MultiAssayExpriment object
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



#' Combine samples with common genes from selected objects
#' @name CombineObjects
#' @param object_list A list contains expression data with mapped gene symbol
#' @param gse_name A vector contains the name of the objects that you want to combine
#' @return A SummarizedExperiment Object contains combined data from several objects
#'
#'
#' @export
CombineObjects <- function(object_list,gse_name=NULL){
  if(is.null(gse_name)){
    gse_name = names(object_list)
    dat_exprs_match <- lapply(object_list, function(x) experiments(x)[["assay_reduce"]] %>% data.frame)
  }
  else {
    dat_exprs_match <- lapply(object_list[gse_name], function(x) experiments(x)[["assay_reduce"]] %>% data.frame)
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
