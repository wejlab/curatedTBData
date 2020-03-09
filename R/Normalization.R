#' Normalization for microarray and RNA-seq transcriptome data.

setGeneric(name="NormalizeReads", function(theObject,...){
  standardGeneric("NormalizeReads")
})

setMethod("NormalizeReads",
          signature="SummarizedExperiment",

          function(theObject,experiment_type = "assay_raw",method = "quantile"){
            # set counts less than 1 to be 1.
            counts <- assays(theObject)[[1]]
            counts[counts<1] <- 1

            # Normalize between arrays
            norm_counts <- limma::normalizeBetweenArrays (counts,
                                                          method = method)
            assays(theObject)[["NormalizedData"]] <- norm_counts
            return(theObject)
          }
)

setMethod("NormalizeReads",
          signature = "MultiAssayExperiment",
          function(theObject, experiment_type = "assay_raw",method = "TMM"){
            # Get raw counts from assay_raw experiment for a MultiAssayExperiment Object
            if (experiment_type == "assay_raw"){
              counts <- assays(experiments(theObject)[[experiment_type]])[[1]]
              counts[counts<1] <- 1
              NormFactor <- calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              assays(experiments(theObject)[[experiment_type]])[["NormalizedData"]] <- Exp
              return(theObject)

            }
            if (experiment_type == "assay_reprocess"){
              counts <- experiments(theObject)[[experiment_type]]
              counts[counts<1] <- 1
              NormFactor <- calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              # add the normalized data as a new experiment in the MultiAssay

              theObject_norm <- c(theObject, assay_reprocess_norm = Exp)

              return(theObject_norm)
            }

          }
)

#' Combine samples with common genes from selected objects
#' @name CombineObjects
#' @param object_list A list contains expression data with mapped gene symbol
#' @param gse_name A vector contains the name of the objects that you want to combine
#' @return A SummarizedExperiment Object contains combined data from several objects
#'
#'
#' @export
CombineObjects <- function(object_list,gse_name){
  dat_exprs_match <- lapply(object_list[gse_name], function(x) experiments(x)[["assay_reduce"]] %>% data.frame)

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
  col_info <- plyr:::rbind.fill(lapply(col_data,function(x){as.data.frame(x)}))
  row.names(col_info) <- Sample

  result <- SummarizedExperiment::SummarizedExperiment(assays = list(as.matrix(dat_exprs_count)),
                                                       colData = col_info)
  return(result)

  }
