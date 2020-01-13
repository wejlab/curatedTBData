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

    ## Create expression matrix with gene symbol, take mean values for same genes
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
