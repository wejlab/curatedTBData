#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' S4 method matches probeID to gene symbol by creating MultiAssayExperiement object from either SummarizedExperiment or MultiAssayExperiment Object.
#' @name MatchProbe
#' @param theObject A SummarizedExperiment/MultiAssayExperiment object.
#' @param experiment_type A character indicates the name of the experiment within MultiAssayExperiment object,
#' experiment_type = assay_raw/assay_raw_norm for matching probe to gene symbol using non-normalized data or normalized data respectively.
#' @param ... Extra named arguments passed to function.
#' @rdname MatchProbe-methods
#' @exportMethod MatchProbe

setGeneric(name="MatchProbe", function(theObject,...){
  standardGeneric("MatchProbe")
})

#' @rdname MatchProbe-methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("MatchProbe",
          signature="SummarizedExperiment",
          function(theObject, experiment_type = "assay_raw"){
            if (experiment_type == "assay_raw"){
              experiment_type = 1
            }
            else{
              experiment_type = "NormalizedData"
            }

            sobject_exprs <- SummarizedExperiment::assays(theObject)[[experiment_type]]

            #### row data is NULL, special case for those normalized data with unique gene symbol as row names
            if (ncol(SummarizedExperiment::rowData(theObject)) == 0){

              ### A special case for those normalized data with unique gene symbol as row names GSEXXXX
              ## Matching process has already done
              mobject1 <- new("Mobject", assay_reprocess = assays(theObject)[[experiment_type]], assay_raw = assays(theObject)[[experiment_type]],
                              row_data = data.frame(theObject@elementMetadata), primary = data.frame(theObject@colData),
                              meta_data = theObject@metadata[[1]])
              simpleMultiAssay <- CreateObject(mobject1,assay_type = "assay_reduce")
              return(simpleMultiAssay)

            }

            # For regular cases
            row_data <- SummarizedExperiment::rowData(theObject) %>% data.frame()

            if (!any(colnames(row_data)=="SYMBOL_NEW")){
              stop("RowData of the input Summarized Experiment Object does not have SYMBOL_NEW column, add SYMBOL_NEW column that includes gene symbols.")
            }
            if (!any(colnames(row_data)=="ID_REF")){
              stop("RowData of the input Summarized Experiment Object does not have ID_REF column, add ID_REF column that includes probe IDs.")
            }
            if (!all(row.names(assays(theObject)[[1]])==row_data$ID_REF)){
              stop("Input Summarized Experiment Object row names are not exactly the same with ID_REF from row Data, consider change")
            }
            # Add new column to the expression matrix

            sobject_exprs_new <- sobject_exprs %>% dplyr::as_tibble() %>% dplyr::mutate(SYMBOL=row_data$SYMBOL_NEW) %>% dplyr::filter(SYMBOL!='NA')

            # Expand probe sets for non-specific probes if apllicable
            if(!identical(grep("///",sobject_exprs_new$SYMBOL), integer(0))){
              sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
            }

            # Create expression matrix with gene symbol, take mean values for same genes
            # Where the function spends most of the time, use dtplyr to speed up. Reduce ~1/3 time using dtplyr
            sobject_exprs_new1 <- dtplyr::lazy_dt(sobject_exprs_new)
            sobject_exprs_symbol <- sobject_exprs_new1 %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise_all(mean) %>% as.data.frame()
            row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

            sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')] %>% as.matrix()

            ## Create MultiAssayExperiment object Use methods from CreateSobject
            # assay_reprocess changes to name of assay_reduce in here
            mobject1 <- new("Mobject", assay_reprocess = sobject_exprs_symbol, assay_raw = assays(theObject)[[experiment_type]],
                            row_data = data.frame(theObject@elementMetadata), primary = data.frame(theObject@colData),
                            meta_data = theObject@metadata[[1]])

            simpleMultiAssay <- CreateObject(mobject1,assay_type = "assay_reduce")
            return(simpleMultiAssay)

          }
)

#' @rdname MatchProbe-methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setMethod("MatchProbe",
          signature="MultiAssayExperiment",
          function(theObject, experiment_type = c("assay_raw","assay_raw_norm")){
              if(experiment_type == "assay_raw_norm"){
                assay_name = "NormalizedData"
              }
              if(experiment_type == "assay_raw"){
                assay_name = 1
              }

              sobject_ori <- MultiAssayExperiment::experiments(theObject)[["assay_raw"]]
              sobject_exprs <- MultiAssayExperiment::assays(sobject_ori)[[assay_name]]
              row_data <- SummarizedExperiment::rowData(sobject_ori) %>% data.frame()
              if (!any(colnames(row_data)=="SYMBOL_NEW")){
                stop("RowData of the input Summarized Experiment Object does not have SYMBOL_NEW column, add SYMBOL_NEW column that includes gene symbols.")
              }
              if (!any(colnames(row_data)=="ID_REF")){
                stop("RowData of the input Summarized Experiment Object does not have ID_REF column, add ID_REF column that includes probe IDs.")
              }
              if (!all(row.names(MultiAssayExperiment::assays(sobject_ori)[[1]])==row_data$ID_REF)){
                stop("Input Summarized Experiment Object row names are not exactly the same with ID_REF from row Data, consider change")
              }

              # Create Multi-assay experiment

              ## Add new column to the expression matrix
              sobject_exprs_new <- sobject_exprs %>% dplyr::as_tibble() %>% dplyr::mutate(SYMBOL=row_data$SYMBOL_NEW) %>% dplyr::filter(SYMBOL!='NA')

              ## Expand probe sets for non-specific probes if apllicable
              if(!identical(grep("///",sobject_exprs_new$SYMBOL), integer(0))){
                sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
              }

              # Create expression matrix with gene symbol, take mean values for same genes
              sobject_exprs_new1 <- dtplyr::lazy_dt(sobject_exprs_new)
              sobject_exprs_symbol <- sobject_exprs_new1 %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise_all(mean) %>% as.data.frame()
              row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

              sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')] %>% as.matrix()

              simpleMultiAssay <- c(theObject,assay_reduce=sobject_exprs_symbol)
              return(simpleMultiAssay)


          })

#' Expand probe set for non-specific probes.
#' @name expandProbesets
#'
#' @param dat_example A matrix with one column represents gene symbol with column name `SYMBOL`, and rest of columns are gene expression value from each sample.
#' @return A matrix that split non-uniquely mapped features to one per row.
#' @examples
#' dat_example <- data.frame(SYMBOL=c("OR7E14P///OR7E12P","CEP104///LILRA6",
#' "SNAR-A1///SNAR-A2///SNAR-A12","ALG6"),sample1=rnorm(4),sample2=rnorm(4))
#' expandProbesets(dat_example, sep="///")
#' @export
expandProbesets <- function(sobject_exprs_new, sep="///"){

  # Get index with duplicated symbol
  index <- grep("///",sobject_exprs_new$SYMBOL)
  sobject_exprs_dup <- sobject_exprs_new[index,]

  x_list <- strsplit(as.character(sobject_exprs_dup$SYMBOL), sep)
  symbol_dup <- gsub(" ","",unlist(x_list))
  sobject_exprs_dup$times <- sapply(x_list, length)

  # Expand rows based on their frequency
  sobject_exprs_expand <- as.data.frame(lapply(sobject_exprs_dup, rep, sobject_exprs_dup$times))
  sobject_exprs_expand$SYMBOL <- symbol_dup
  sobject_exprs_expand_result <- sobject_exprs_expand %>% dplyr::select(-times)

  # double-check SYMBOL does not contain NA
  if(any(is.na(sobject_exprs_expand_result$SYMBOL))){
    stop("Expand probe sets contain NA's, please check")
  }

  sobject_exprs_result <- rbind(sobject_exprs_new[-index,],sobject_exprs_expand_result)
  return(sobject_exprs_result)

}

