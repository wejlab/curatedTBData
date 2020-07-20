#' S4 method matches probeID to gene symbol by creating MultiAssayExperiement object from either SummarizedExperiment or MultiAssayExperiment Object.
#' @name MatchProbe
#' @param theObject A SummarizedExperiment/MultiAssayExperiment Object.
#' @param UseAssay A character indicats the assay names (partial name) of the SummarizedExperiment Object.
#' @param createExperimentName A character specifying the names of the new experiment matrix.
#' @param only.matrix A logical value. Default is FALSE. TRUE for only output the matched matrix.
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
          function(theObject, UseAssay, createExperimentName = "assay_MatchProbe",
                   only.matrix=FALSE){

            UseAssay <- paste(UseAssay,collapse = "|")
            assay_name_index <- grep(UseAssay,SummarizedExperiment::assayNames(theObject))

            if(length(assay_name_index) > 1){
              assay_name_index <- assay_name_index[1]

              message(paste("More than one assay selected. Use assay:",
                            SummarizedExperiment::assayNames(theObject)[assay_name_index],
                            "for probes matching." ))
            }
            if(length(assay_name_index) == 0){
              stop(paste("No assay with the name:",UseAssay))
            }

            assay_name <- SummarizedExperiment::assayNames(theObject)[assay_name_index]


            geo_access_name <- strsplit(SummarizedExperiment::assayNames(theObject),
                                        "_")[[1]][1]

            sobject_exprs <- SummarizedExperiment::assays(theObject)[[assay_name]]

            #### row data is NULL, special case for those normalized data with unique gene symbol as row names
            if (ncol(SummarizedExperiment::rowData(theObject)) == 0){


              ### A special case for those normalized data with unique gene symbol as row names GSEXXXX
              ## Matching process has already done
              mobject1 <- methods::new("Mobject", assay_reprocess = SummarizedExperiment::assays(theObject)[[assay_name]],
                              assay_raw = SummarizedExperiment::assays(theObject)[[assay_name]],
                              row_data = S4Vectors::DataFrame(theObject@elementMetadata),
                              primary = S4Vectors::DataFrame(theObject@colData),
                              meta_data = theObject@metadata[[1]])
              simpleMultiAssay <- CreateObject(mobject1, createExperimentName = createExperimentName)

              sobject1_final <- simpleMultiAssay[["assay_raw"]]

              # Provide assay names
              names(SummarizedExperiment::assays(sobject1_final)) <- assay_name
              simpleMultiAssay[["assay_raw"]] <- sobject1_final

              return(simpleMultiAssay)

            }

            # For regular cases
            row_data <- SummarizedExperiment::rowData(theObject)

            if (!any(colnames(row_data)=="SYMBOL_NEW")){
              stop("RowData of the input Summarized Experiment Object does not have SYMBOL_NEW column, add SYMBOL_NEW column that includes gene symbols.")
            }
            if (!any(colnames(row_data)=="ID_REF")){
              stop("RowData of the input Summarized Experiment Object does not have ID_REF column, add ID_REF column that includes probe IDs.")
            }
            if (!all(row.names(SummarizedExperiment::assays(theObject)[[1]])==row_data$ID_REF)){
              stop("Input Summarized Experiment Object row names are not exactly the same with ID_REF from row Data, consider change")
            }
            # Add new column to the expression matrix

            sobject_exprs_new <- sobject_exprs %>% dplyr::as_tibble() %>%
              dplyr::mutate(SYMBOL=row_data$SYMBOL_NEW) %>%
              dplyr::filter(.data$SYMBOL != "NA")

            # Expand probe sets for non-specific probes if apllicable
            if(length(grep("///",sobject_exprs_new$SYMBOL)) != 0){
              sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
            } else if(length(grep(";",sobject_exprs_new$SYMBOL)) != 0){
              sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = ";")
            }

            # Create expression matrix with gene symbol, take mean values for same genes
            # Where the function spends most of the time, use dtplyr to speed up. Reduce ~1/3 time using dtplyr

            # sobject_exprs_new1 <- dtplyr::lazy_dt(sobject_exprs_new)
            sobject_exprs_new1 <- sobject_exprs_new

            ##### Think how to reduce processing time for the following code??
            #sobject_exprs_symbol <- sobject_exprs_new1 %>%
            #  dplyr::group_by(.data$SYMBOL) %>%
            #  dplyr::summarise_all(mean) %>% as.data.frame()

            # Try aggregate from stats package. Avoid using sobject_exprs_new1$SYMBOL, slow down the process
            sobject_exprs_symbol <- stats::aggregate(. ~ SYMBOL, data = sobject_exprs_new1,
                                                    FUN = mean)

            row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

            sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol)
                                                 %in% 'SYMBOL')] %>% as.matrix()

            if(only.matrix){return(sobject_exprs_symbol)}

            ## Create MultiAssayExperiment object Use methods from CreateSobject
            # assay_reprocess changes to name of assay_reduce in here
            mobject1 <- new("Mobject", assay_reprocess = sobject_exprs_symbol,
                            assay_raw = SummarizedExperiment::assays(theObject)[[assay_name]],
                            row_data = S4Vectors::DataFrame(theObject@elementMetadata),
                            primary = S4Vectors::DataFrame(theObject@colData),
                            meta_data = theObject@metadata[[1]])

            simpleMultiAssay <- CreateObject(mobject1,createExperimentName = createExperimentName)

            sobject1_final <- simpleMultiAssay[["assay_raw"]]
            names(SummarizedExperiment::assays(sobject1_final)) <- assay_name
            simpleMultiAssay[["assay_raw"]] <- sobject1_final

            return(simpleMultiAssay)

          }
)

#' @rdname MatchProbe-methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setMethod("MatchProbe",
          signature="MultiAssayExperiment",
          function(theObject, UseAssay, createExperimentName = "assay_MatchProbe",
                   only.matrix=FALSE){

              sobject_ori <- MultiAssayExperiment::experiments(theObject)[["assay_raw"]]

              UseAssay <- paste(UseAssay,collapse = "|")
              assay_name_index <- grep(UseAssay,SummarizedExperiment::assayNames(sobject_ori))

              if(length(assay_name_index) > 1){
                assay_name_index <- assay_name_index[1]

                message(paste("More than one assay selected. Use assay:",
                              SummarizedExperiment::assayNames(theObject)[assay_name_index],
                              "for probes matching." ))
              }
              if(length(assay_name_index) == 0){
                stop(paste("No assay with the name:",UseAssay))
              }

              assay_name <- SummarizedExperiment::assayNames(sobject_ori)[assay_name_index]

              if(length(assay_name_index) != 1){
                stop("More than one assay selected, please consider only one assay.")
              }

              sobject_exprs <- SummarizedExperiment::assays(sobject_ori)[[assay_name]]
              row_data <- SummarizedExperiment::rowData(sobject_ori)
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
              sobject_exprs_new <- sobject_exprs %>% dplyr::as_tibble() %>%
                dplyr::mutate(SYMBOL=row_data$SYMBOL_NEW) %>%
                dplyr::filter(.data$SYMBOL != 'NA')

              ## Expand probe sets for non-specific probes if apllicable
              if(length(grep("///",sobject_exprs_new$SYMBOL)) != 0){
                sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = "///")
              }else if(length(grep(";",sobject_exprs_new$SYMBOL)) != 0){
                sobject_exprs_new <- expandProbesets(sobject_exprs_new, sep = ";")
              }

              # Create expression matrix with gene symbol, take mean values for same genes
              # sobject_exprs_new1 <- dtplyr::lazy_dt(sobject_exprs_new)
              sobject_exprs_new1 <- sobject_exprs_new

              #sobject_exprs_symbol <- sobject_exprs_new1 %>% dplyr::group_by(.data$SYMBOL) %>%
              #  dplyr::summarise_all(mean) %>% as.data.frame()
              # %>% increases time, try aggregate from stats package
              sobject_exprs_symbol <- stats::aggregate(. ~ SYMBOL, data = sobject_exprs_new1,
                                                      FUN = mean)

              row.names(sobject_exprs_symbol) <- sobject_exprs_symbol$SYMBOL

              sobject_exprs_symbol <- sobject_exprs_symbol[,-which(colnames(sobject_exprs_symbol) %in% 'SYMBOL')] %>% as.matrix()

              if(only.matrix){
                return(sobject_exprs_symbol)
              }

              simpleMultiAssay <- c(theObject,assay_reduce=sobject_exprs_symbol)
              index_name <- which(names(simpleMultiAssay) %in% "assay_reduce")
              names(simpleMultiAssay)[index_name] <- createExperimentName

              return(simpleMultiAssay)

          })

#' Expand probe set for non-specific probes.
#' @name expandProbesets
#' @param dat A matrix with one column represents gene symbol with
#' column name "SYMBOL", and rest of columns are gene expression value from each sample.
#' @param sep A character string that separates gene symbols.
#' @return A matrix that split non-uniquely mapped features to one per row.
#' dat_example <- data.frame(SYMBOL=c("OR7E14P///OR7E12P","CEP104///LILRA6",
#' "SNAR-A1///SNAR-A2///SNAR-A12","ALG6"),sample1=rnorm(4),sample2=rnorm(4))
#' expandProbesets(dat_example, sep="///")
#' @export
expandProbesets <- function(dat, sep){

  times <- NULL
  # Get index with duplicated symbol
  index <- grep(sep,dat$SYMBOL)
  sobject_exprs_dup <- dat[index,]

  x_list <- strsplit(as.character(sobject_exprs_dup$SYMBOL), sep)
  symbol_dup <- gsub(" ","",unlist(x_list))
  sobject_exprs_dup$times <- unlist(lapply(x_list, length))

  # Expand rows based on their frequency
  sobject_exprs_expand <- as.data.frame(lapply(sobject_exprs_dup, rep,
                                               sobject_exprs_dup$times))
  sobject_exprs_expand$SYMBOL <- symbol_dup
  sobject_exprs_expand_result <- sobject_exprs_expand %>% dplyr::select(-times)

  # double-check SYMBOL does not contain NA
  if(any(is.na(sobject_exprs_expand_result$SYMBOL))){
    stop("Expand probe sets contain NA's, please check")
  }

  sobject_exprs_result <- rbind(dat[-index,],sobject_exprs_expand_result)
  return(sobject_exprs_result)

}

