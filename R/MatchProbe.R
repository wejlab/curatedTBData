#' Match ProbeID to gene symbol by creating MultiAssayExperiement object from SummarizedExperiment/MultiAssayExperiment Object.
#' @name MatchProbe
#'
#' @param sobject A summarized experiment object
#' @return A MultiAssay Experiment Object contains 1. Input Summarized Experiment Object. 2. read count matrix collapse by gene symbol
#' @examples
#' dat("GSE39939_sobject")
#' MatchProbe(GSE39939_sobject)
#' @export
MatchProbe <- function(sobject, experiment_type = "raw"){

  experiment_type <- match.arg(experiment_type)
  if (experiment_type == "raw"){
     experiment_type = 1
  }
  if (experiment_type != "raw"){
    # execute normalized reaads
    experiment_type_sobject = "NormalizedReads"
    experiment_type_mobject1 = "assay_reprocess_norm"
    experiment_type_mobject2 = "assay_raw_norm"
  }
  ## For SummarizedExperiment object
  if (class(sobject) == "SummarizedExperiment"){
    # read in data
    sobject_exprs <- assays(sobject)[[experiment_type]]
    #### row data is NULL, special case for those normalized data with unique gene symbol as row names
    if (ncol(rowData(sobject)) == 0){
      #### row data is NULL, special case for those normalized data with unique gene symbol as row names

      mobject1 <- new("Mobject", assay_reprocess = sobject@assays[[experiemnt_type]], assay_raw = sobject@assays[[experiment_type]],
                      row_data = data.frame(sobject@elementMetadata), primary = data.frame(sobject@colData),
                      meta_data = sobject@metadata[[1]])
      simpleMultiAssay <- CreateObject(mobject1,assay_type = "assay_reduce")
      return(simpleMultiAssay)

    }
    row_data <- rowData(sobject) %>% data.frame()
    if (!any(colnames(row_data)=="SYMBOL_NEW")){
      stop("RowData of the input Summarized Experiment Object does not have SYMBOL_NEW column, add SYMBOL_NEW column that includes gene symbols.")
    }
    if (!any(colnames(row_data)=="ID_REF")){
      stop("RowData of the input Summarized Experiment Object does not have ID_REF column, add ID_REF column that includes probe IDs.")
    }
    if (!all(row.names(assay(sobject))==row_data$ID_REF)){
      stop("Input Summarized Experiment Object row names are not exactly the same with ID_REF from row Data, consider change")
    }

    # Add new column to the expression matrix
    sobject_exprs_new <- sobject_exprs %>% as_tibble() %>% mutate(SYMBOL=row_data$SYMBOL_NEW) %>% filter(SYMBOL!='NA')

    # Expand probe sets for non-specific probes if apllicable
    if(!identical(grep("///",sobject_exprs_new$SYMBOL), integer(0))){
      sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
    }

    # Create expression matrix with gene symbol, take mean values for same genes
    sobject_exprs_symbol <- sobject_exprs_new %>% group_by(SYMBOL) %>% summarise_all(mean) %>% data.frame()
    row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

    sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')] %>% as.matrix()

    ## Create MultiAssayExperiment object Use methods from CreateSobject
    mobject1 <- new("Mobject", assay_reprocess = sobject_exprs_symbol, assay_raw = sobject@assays[[experiment_type]],
                    row_data = data.frame(sobject@elementMetadata), primary = data.frame(sobject@colData),
                    meta_data = sobject@metadata[[1]])

    simpleMultiAssay <- CreateObject(mobject1,assay_type = "assay_reduce")
    return(simpleMultiAssay)
  }

  ## For MultiAssayExperiment object
  if (class(sobject) == "MultiAssayExperiment"){
    sobject_ori <- MultiAssayExperiment::experiments(sobject)[["assay_raw"]]
    sobject_exprs <- assay(sobject_ori)
    row_data <- rowData(sobject_ori) %>% data.frame()
    if (!any(colnames(row_data)=="SYMBOL_NEW")){
      stop("RowData of the input Summarized Experiment Object does not have SYMBOL_NEW column, add SYMBOL_NEW column that includes gene symbols.")
    }
    if (!any(colnames(row_data)=="ID_REF")){
      stop("RowData of the input Summarized Experiment Object does not have ID_REF column, add ID_REF column that includes probe IDs.")
    }
    if (!all(row.names(assay(sobject_ori))==row_data$ID_REF)){
      stop("Input Summarized Experiment Object row names are not exactly the same with ID_REF from row Data, consider change")
    }

    ## Starting create Multi-assay experiment

    # Add new column to the expression matrix
    sobject_exprs_new <- sobject_exprs %>% as_tibble() %>% mutate(SYMBOL=row_data$SYMBOL_NEW) %>% filter(SYMBOL!='NA')

    # Expand probe sets for non-specific probes if apllicable
    if(!identical(grep("///",sobject_exprs_new$SYMBOL), integer(0))){
      sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
    }

    # Create expression matrix with gene symbol, take mean values for same genes
    sobject_exprs_symbol <- sobject_exprs_new %>% group_by(SYMBOL) %>% summarise_all(mean) %>% data.frame()
    row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

    sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')] %>% as.matrix()

    simpleMultiAssay <- c(sobject,assay_reduce=sobject_exprs_symbol)
    return(simpleMultiAssay)
  }

}

#' Expand probe set for non-specific probes
#' @name expandProbesets
#'
#' @param dat_example A matrix with one column represents gene symbol with column name `SYMBOL`, and rest of columns are gene expression value from each sample
#' @return A matrix that split non-uniquely mapped features to one per row
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

