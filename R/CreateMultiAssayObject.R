#' Create MultiAssay Experiment Object from Summarized Experiment Object.
#' @name Create_MultiAssay_object
#'
#' @param sobject A summarized experiment object
#' @return A MultiAssay Experiment Object contains 1. Input Summarized Experiment Object. 2. read count matrix collapse by gene symbol
#' @examples
#' dat("GSE39939_sobject")
#' Create_MultiAssay_object(GSE39939_sobject)
#' @export
Create_MultiAssay_object <- function(sobject){

  # check input type
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (! "MultiAssayExperiment" %in% installed.packages()) BiocManager::install("MultiAssayExperiment")
  if(class(sobject) != "SummarizedExperiment"){
    stop(paste("Invalid input data type. Only supported for SummarizedExperiment objects. Your input:,", typeof(sobject)))
  }
  library(MultiAssayExperiment)
  library(dplyr)
  # Read in Data
  sobject_exprs <- assay(sobject)
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

    sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')]

    ## Create MultiAssayExperiment object
    doubleExp <- list('Probe' = sobject, 'Symbol' = sobject_exprs_symbol)
    symbol_map <- data.frame(primary = row.names(colData(sobject)),
                             colname = row.names(colData(sobject)), stringsAsFactors = FALSE)
    probe_map <- data.frame(primary_info = rep("n.a.",length(row.names(colData(sobject)))),stringsAsFactors = FALSE)
    row.names(probe_map) <- row.names(colData(sobject))

    listmap <- list(symbol_map,symbol_map)
    names(listmap) <- c("Probe", "Symbol")
    dfmap <- listToMap(listmap)

    simpleMultiAssay <- MultiAssayExperiment(experiments=doubleExp,
                                             colData=probe_map,dfmap)

    return(simpleMultiAssay)

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

