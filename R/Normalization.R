#' Normalization for microarray and RNA-seq transcriptome data.

setGeneric(name="NormalizeReads", function(theObject,...){
  standardGeneric("NormalizeReads")
})

setMethod("NormalizeReads",
          signature="SummarizedExperiment",

          function(theObject, method = "quantile"){
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
              NormFactor <- calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              assays(experiments(theObject)[[experiment_type]])[["NormalizedData"]] <- Exp
              return(theObject)

            }
            if (experiment_type == "assay_reprocess"){
              counts <- experiments(theObject)[[experiment_type]]

              NormFactor <- calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              # add the normalized data as a new experiment in the MultiAssay

              theObject_norm <- c(theObject, assay_reprocess_norm = Exp)

              return(theObject_norm)
            }


          }
)

