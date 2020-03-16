#' Normalization for microarray and RNA-seq transcriptome data.

setGeneric(name="Normalization", function(theObject,...){
  standardGeneric("Normalization")
})

# theObject <- GSE39939_sobject

setMethod("Normalization",
          signature="SummarizedExperiment",

          function(theObject,method = "quantile", experiment_type = "assay_raw"){

            # check whether the object has been normalized.
            norm_GSE <- paste(c("GSE54992","GSE36238","GSE31348","GSE73408","GSE41055" ,"GSEXXXXX"),
                              collapse="|")
            if (length(grep(norm_GSE,assayNames(theObject))) == 1){
              assays(theObject)[["NormalizedData"]] <- assays(theObject)[[1]]
              return(theObject)

            }

            # Remove outliers for microarray????
            ################################# Function to be inserted ############

            # set counts less than 10 to be 10.
            counts <- assays(theObject)[[1]]
            counts[counts<10] <- 10

            # log2 transformed data
            counts <- log(counts,base=2)
            # Normalize between arrays
            norm_counts <- limma::normalizeBetweenArrays (counts,
                                                          method = method)
            assays(theObject)[["NormalizedData"]] <- norm_counts
            return(theObject)
          }
)

setMethod("Normalization",
          signature = "MultiAssayExperiment",
          function(theObject, experiment_type = "assay_raw",method = "TMM"){
            # Get raw counts from assay_raw experiment for a MultiAssayExperiment Object
            if (experiment_type == "assay_raw"){
              counts <- assays(experiments(theObject)[[experiment_type]])[[1]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts,base=2)

              NormFactor <- edgeR::calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              assays(experiments(theObject)[[experiment_type]])[["NormalizedData"]] <- Exp
              return(theObject)

            }
            if (experiment_type == "assay_reprocess"){

              counts <- experiments(theObject)[[experiment_type]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts,base=2)

              NormFactor <- calcNormFactors(counts, method = method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              # add the normalized data as a new experiment in the MultiAssay

              theObject_norm <- c(theObject, assay_reprocess_norm = Exp)

              return(theObject_norm)
            }

          }
)

#' Remove outliers using arrayQualityMetrics R package
#'
#'
