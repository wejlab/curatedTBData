#' Normalization for microarray and RNA-seq.
#' @name norm_reads
#'
#' @param object A SummarizedExperiment/MultiAssayExperiment object
#' @param method
#' @param

#' @return A summarized experiment object for TB signature profiling
#'
#' @examples
#' sobject <- dat("GSE39939_sobject")
#' mobject <- Create_MultiAssay_object(sobject)
#' sobject_TBSig <- get_sobject_TBSig(mobject,"PTB","Latent")
#' @export

setGeneric(name="NormalizeReads", function(theObject,...){
  standardGeneric("NormalizeReads")
})

setMethod("NormalizeReads",
          signature="SummarizedExperiment",
          function(theObject, method = "quantile"){
            norm_counts <- limma::normalizeBetweenArrays (assays(theObject)[[1]],
                                                          method = method)
            assays(theObject)[["NormalizedReads"]] <- norm_counts
            return(theObject)
          }
)

setMethod("NormalizeReads",
          signature = "MultiAssayExperiment",
          function(theObject, experiment_type = "assay_raw",method = "TMM"){
            # Get raw counts from assay_raw experiment for a MultiAssayExperiment Object
            counts <- assay(experiments(theObject)[[experiment_type]])
            NormFactor <- calcNormFactors(counts, method = method)
            ScaleFactors <- colSums(counts) * NormFactor

            Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
            # add the normalized data as a new experiment in the MultiAssay

            norm <- paste0(experiment_type,"_norm")

            theObject_norm <- c(theObject, norm = Exp)

            # rename the normalized experiment
            names(theObject_norm)[which(names(theObject_norm) == "norm")] <- norm

            return(theObject_norm)

          }
)

