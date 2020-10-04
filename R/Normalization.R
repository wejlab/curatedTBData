#' @title Normalization for curatedTBData
#' @description Perform normalization for curatedTBDta.
#' \code{microarray_method} (except "SCAN") argument passes to
#' \code{\link[limma]{normalizeBetweenArrays}}
#' for microarray data. \code{RNAseq_method argument passes to
#' \code{\link[edgeR]{calcNormFactors}}}
#' for RNA sequencing data.
#' @param theObject A \code{SummarizedExperiment} or \code{MultiAssayExperiment} object.
#' @param geo_access GEO accession number for Tuberculosis Study. Used when
#' microarray_method is "SCAN".
#' @param microarray_method A \code{character} string specifying the normalization
#' method to be used.
#' See \code{\link[limma]{normalizeBetweenArrays}} for deatils. The default is "quantile".
#' @param RNAseq_method A \code{character} string specifying the normalization
#' method to be used.
#' See \code{\link[edgeR]{calcNormFactors}} for RNA sequence data. The default is "TMM".
#' @param experiment_name A \code{character} indicates the name of the experiment
#' within MultiAssayExperiment object. Only applicable in multiple assays conditions.
#' Choices for experiment_name are "assay_raw" and "assay_reprocess".
#' The deafult is experiment_name="assay_raw"
#' If experiment_name is "assay_raw", perform normalization on the assay provided by the authors.
#' If experiment_name is "assay_reprocess", perform normalization on the reprocessed assay.
#' @param ... Extra named arguments passed to function.
#' @return A \code{SummarizedExperiment} object with additional normalized assay if
#' input is a \code{SummarizedExperiment} object.
#' A MultiAssayExperiment object with additional normalized assay if input is a
#' \code{MultiAssayExperiment} object.
#' @examples
#' object_list <- get_curatedTBData(geo_access = c("GSE39939","GSE107993"))
#' GSE39939_sobject_norm <- Normalization(object_list$GSE39939, microarray_method = "quantile")
#' GSE107993_sobject_norm <- Normalization(object_list$GSE107993, RNAseq_method = "TMM")
#' @rdname Normalization-methods
#' @exportMethod Normalization
setGeneric(name="Normalization", function(theObject,...){
  standardGeneric("Normalization")
})

#' @rdname Normalization-methods
setMethod("Normalization",
          signature="SummarizedExperiment",

          function(theObject, geo_access = NULL, microarray_method = "quantile",
                   RNAseq_method = NULL, ...){

            # Import Data Summarized table into the function
            # DataSummary <- get(data("DataSummary",package="curatedTBData"))
            geo_access_name <- strsplit(SummarizedExperiment::assayNames(theObject),
                                        "_")[[1]][1]

            # SCAN takes long time to process, we include the preprocessed data
            # in the data package. Include function in case users need
            if (microarray_method == "SCAN"){
              message("Use SCAN normalization method, may take long time")
              param <- BiocParallel::SerialParam(progressbar = TRUE)
              gse <- suppressMessages(GEOquery::getGEO(geo_access, GSEMatrix = FALSE))
              sample_name <- names(GEOquery::GSMList(gse))

              SCAN_list <- BiocParallel::bplapply(sample_name,
                                     function (x) SCAN.UPC::SCAN(x), BPPARAM = param)
              dat_exprs_match <- lapply(SCAN_list, function(x) data.frame(Biobase::exprs(x)))
              dat_exprs_combine <- Reduce(
                function(x, y) merge(x, y, by = "id", all = FALSE),
                lapply(dat_exprs_match, function(x) { x$id <- rownames(x); x }))
              rowname_data <- dat_exprs_combine$id

              dat_final <- dat_exprs_combine %>% dplyr::select(-.data$id) %>%
                                                 as.matrix()
              row.names(dat_final) <- rowname_data
              colname_data <- colnames(dat_final)
              index <- lapply(sample_name, function(x)
                                           grep(x,colname_data)) %>% unlist()

              colnames(dat_final) <- sample_name[index]

              unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

              assay_name <- paste0(geo_access_name, "_",microarray_method)

              SummarizedExperiment::assays(theObject)[[assay_name]] <- dat_final
              return(theObject)

            }

            # Use normalizeBetweenArrays from limma package for normalization

            else{
              # Check whether the object is derived from Affymetrix or GSEXXXXX.

              norm_GSE <- paste(c("GSE31348", "GSE36238", "GSE41055", "GSE54992",
                                  "GSE73408", "GSE107731"), collapse="|")
              if (length(grep(norm_GSE,
                              SummarizedExperiment::assayNames(theObject))) != 0){

                assay_name <- paste0(geo_access_name, "_","RMA")

                SummarizedExperiment::assays(theObject)[[assay_name]] <- SummarizedExperiment::assays(theObject)[[1]]
                return(theObject)

              }

              # set counts less than 10 to be 10.
              counts <- SummarizedExperiment::assays(theObject)[[1]]
              counts[counts<10] <- 10
              counts <- log(counts, base=2) # log2 transformed data

              # Normalize between arrays
              norm_counts <- limma::normalizeBetweenArrays (counts,
                                                     method = microarray_method)

              assay_name <- paste0(geo_access_name, "_",microarray_method)
              SummarizedExperiment::assays(theObject)[[assay_name]] <- norm_counts
              return(theObject)
            }

          }
)

#' @rdname Normalization-methods
setMethod("Normalization",
          signature = "MultiAssayExperiment",
          function(theObject, geo_access = NULL, microarray_method = NULL,
                   RNAseq_method = "TMM",
                   experiment_name = c("assay_raw","assay_reprocess")){

            # Identify Normalization method
            # microarray_method <-  match.arg(microarray_method)
            # RNAseq_method <- match.arg(RNAseq_method)

            # Get raw counts from assay_raw experiment for a MultiAssayExperiment Object
            if(missing(experiment_name)){experiment_name="assay_raw"}
            experiment_name <- match.arg(experiment_name)

            sobject1 <- theObject[["assay_raw"]]
            geo_access_name <- strsplit(SummarizedExperiment::assayNames(sobject1),"_")[[1]][1]

            if (experiment_name == "assay_raw"){
              assay_name <- paste0(geo_access_name,"_raw")
              counts <- SummarizedExperiment::assays(theObject[[experiment_name]])[[assay_name]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts, base=2)

              NormFactor <- edgeR::calcNormFactors(counts, method = RNAseq_method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))

              assay_name <- paste0(geo_access_name,"_",RNAseq_method)

              SummarizedExperiment::assays(theObject[[experiment_name]])[[assay_name]] <- Exp
              return(theObject)

            }
            if (experiment_name == "assay_reprocess"){

              counts <- theObject[[experiment_name]]
              counts[counts<10] <- 10

              # log2 transformed data
              counts <- log(counts,base=2)

              NormFactor <- edgeR::calcNormFactors(counts, method = RNAseq_method)
              ScaleFactors <- colSums(counts) * NormFactor

              Exp <- round(t(t(counts)/ScaleFactors) * mean(ScaleFactors))

              # Include the normalized data as a new experiment in the MultiAssayExperiment Object
              theObject_norm <- c(theObject, assay_reprocess_norm = Exp)

              return(theObject_norm)
            }

          }
)


