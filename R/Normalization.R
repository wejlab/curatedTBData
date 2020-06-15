#' Normalization for microarray and RNA-seq transcriptome data.
#' @name Normalization
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param microarray_method A character string specifying the normalization method to be used. See `limma::normalizeBetweenArrays` for microarray data.
#' The default is "quantile"
#' @param RNAseq_method A character string specifying the normalization method to be used. See `edgeR::calcNormFactors` for RNA sequence data.
#' The default is "TMM"
#' @param experiment_type A character indicates the name of the experiment within MultiAssayExperiment object. Only applicable in multiple assays conditions.
#' Choices for experiment_type are "assay_raw" and "assay_reprocess".
#' If experiment_type is "assay_raw", perform normalization on the assay provided by the authors.
#' If experiment_type is "assay_reprocess", perform normalization on the reprocessed assay.
#' @param ... Extra named arguments passed to function.
#' @rdname Normalization-methods
#' @exportMethod Normalization
setGeneric(name="Normalization", function(theObject,...){
  standardGeneric("Normalization")
})

# Normalization for microarray
#' @rdname Normalization-methods
setMethod("Normalization",
          signature="SummarizedExperiment",

          function(theObject, geo_access = NULL, microarray_method = "quantile", RNAseq_method = "TMM", ...){

            if (microarray_method == "SCAN"){

                ### Do some thing here
            }

            # Use normalizeBetweenArrays from limma package for normalization

            else{
              # Check whether the object is derived from Affymetrix or GSEXXXXX.

              norm_GSE <- paste(c("GSE54992","GSE36238","GSE31348","GSE73408","GSE41055" ,"GSEXXXXX"),
                                collapse="|")
              if (length(grep(norm_GSE,SummarizedExperiment::assayNames(theObject))) == 1){
                SummarizedExperiment::assays(theObject)[["NormalizedData"]] <- SummarizedExperiment::assays(theObject)[[1]]
                return(theObject)

              }

              # Remove outliers for microarray????
              ############## Function to be inserted ############

              # set counts less than 10 to be 10.
              counts <- SummarizedExperiment::assays(theObject)[[1]]
              counts[counts<10] <- 10
              counts <- log(counts,base=2) # log2 transformed data

              # Normalize between arrays
              norm_counts <- limma::normalizeBetweenArrays (counts, method = microarray_method)
              assays(theObject)[["NormalizedData"]] <- norm_counts
              return(theObject)
            }


          }
)

#' @rdname Normalization-methods
setMethod("Normalization",
          signature = "MultiAssayExperiment",
          function(theObject, microarray_method = "quantile", RNAseq_method = "TMM",
                   experiment_type = c("assay_raw","assay_reprocess")){

            # Identify Normalization method
            # microarray_method <-  match.arg(microarray_method)
            # RNAseq_method <- match.arg(RNAseq_method)

            # Get raw counts from assay_raw experiment for a MultiAssayExperiment Object
            if(missing(experiment_type)){experiment_type="assay_raw"}
            experiment_type <- match.arg(experiment_type)
            if (experiment_type == "assay_raw"){
              counts <- SummarizedExperiment::assays(MultiAssayExperiment::experiments(theObject)[[experiment_type]])[[1]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts,base=2)

              NormFactor <- edgeR::calcNormFactors(counts, method = RNAseq_method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              SummarizedExperiment::assays(MultiAssayExperiment::experiments(theObject)[[experiment_type]])[["NormalizedData"]] <- Exp
              return(theObject)

            }
            if (experiment_type == "assay_reprocess"){

              counts <- MultiAssayExperiment::experiments(theObject)[[experiment_type]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts,base=2)

              NormFactor <- edgeR::calcNormFactors(counts, method = RNAseq_method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))
              # add the normalized data as a new experiment in the MultiAssay

              theObject_norm <- c(theObject, assay_reprocess_norm = Exp)

              return(theObject_norm)
            }

          }
)

# Remove outliers using arrayQualityMetrics R package

