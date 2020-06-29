#' Sobject Class
#' Create SummarizedExperiment object for curatedTBData
#'
#' @slot assay A matrix contatins gene expression data
#' @slot row_data A DataFrame class contains gene information
#' @slot column_data A DataFrame class contains sample information
#' @slot meta_data A MIAME class contains experiment information
#' @rdname Sobject-class
#' @importClassesFrom Biobase MIAME
#'
#' @exportClass Sobject
Sobject <- setClass("Sobject",
                    slots = c(assay = "matrix",row_data = "DataFrame",
                              column_data  = "DataFrame", meta_data = "MIAME"),
                    prototype=list(assay = matrix(c(1, 2, 3, 11, 12, 13),
                                                  nrow = 2, ncol = 3, byrow = TRUE,
                                                  dimnames = list(c("111_at", "222_at"),
                                                                  c("S.1", "S.2", "S.3"))),
                                   row_data = S4Vectors::DataFrame(ID_REF=c("111_at", "222_at"),
                                                                   Symbol=c("A1BC","ZAC"),
                                                                   row.names = c("111_at", "222_at")),
                                   column_data = S4Vectors::DataFrame(Gender=c("Male", "Female","Female"),
                                                                      TBStatus=c("PTB","Latent", "Control"),
                                                                      row.names = c("S.1", "S.2", "S.3")),
                                   # Create  new class in biobase
                                   meta_data = methods::new('MIAME', name="XXXXX", lab="XXXXXX", contact="XXXXX",
                                                   title="A title",abstract="An abstract", url="XXXXXXX",
                                                   pubMedIds = '0000000', other=list(Platform = '000000'))),
                    validity = function(object){
                      if(!all(row.names(object@assay)==row.names(object@row_data))) {
                        return("row names in the assay must be the same as row names
                                          in the row data")
                      }
                      else if (!all(colnames(object@assay)==row.names(object@column_data))) {
                        return("column names in the assay must be the same as row names
                                          in the column data")
                      }
                      else {TRUE}
                    })

#' Mobject Class
#' Create MultiAssayExperiment object for curatedTBData
#'
#' @slot assay_reprocess A matrix contatins gene expression data
#' @slot assay_raw A matrix contatins gene expression data
#' @slot row_data A DataFrame class contains gene expression data with different dimensions
#' @slot primary A DataFrame class contains sample information
#' @slot meta_data A MIAME class contains experiment information
#' @rdname Mobject-class
#'
#' @exportClass Mobject
Mobject <- setClass("Mobject",
                    slots = c(assay_reprocess = "matrix", assay_raw = "matrix",
                              row_data = "DataFrame", primary = "DataFrame",
                              meta_data = "MIAME"),
                    prototype = list(assay_reprocess = matrix(seq_len(12),nrow = 3,
                                                              byrow = TRUE,
                                                              dimnames = list(c("AZA","BBD","CCS"),
                                                                              c("S.1","S.2","S.3","S.4"))),
                                     assay_raw = matrix(c(1, 2, 3, 11, 12, 13),
                                                        nrow = 2, ncol = 3, byrow = TRUE,
                                                        dimnames = list(c("111_at", "222_at"),
                                                                        c("S.1", "S.2", "S.3"))),
                                     row_data = S4Vectors::DataFrame(ID_REF=c("111_at", "222_at"),
                                                                     Symbol=c("A1BC","ZAC"),
                                                                     row.names = c("111_at", "222_at")),
                                     primary = S4Vectors::DataFrame(Gender=c("Male", "Female","Female","Female"),
                                                                    TBStatus=c("PTB","Latent", "Control","PTB"),
                                                                    row.names = c("S.1", "S.2", "S.3","S.4")),
                                     meta_data = methods::new('MIAME', name="XXXXX", lab="XXXXXX",
                                                     contact="XXXXX", title="A title",abstract="An abstract",
                                                     url="XXXXXXX", pubMedIds = '0000000',
                                                     other=list(Platform = '000000'))),
                    validity = function(object){
                      if(!all(row.names(object@assay_raw)==(object@row_data$ID_REF))) {
                        return("row names in the assay must be the same as
                               ID_REF in the row data")
                      }
                    })


#' Create SummarizedExperiment/MultiAssayExperiment object
#' @name CreateObject
#' @param theObject A class either Sobject or Mobject
#' @param createExperimentName A character indicates the name of the new experiment
#' after creating the MultiAssayExperiment Object.
#' @param ... Extra named arguments passed to function
#' @rdname CreateObject-methods
#' @exportMethod CreateObject
setGeneric(name="CreateObject", function(theObject,...){
  standardGeneric("CreateObject")
})

#' @rdname CreateObject-methods
setMethod("CreateObject", signature="Sobject",
          function(theObject){
            results <- SummarizedExperiment::SummarizedExperiment(assays = list(theObject@assay),
                                                                  colData = theObject@column_data,
                                                                  rowData = theObject@row_data,
                                                                  metadata = list(theObject@meta_data))
            return(results)
          }
)

#' @rdname CreateObject-methods
setMethod("CreateObject",
          signature="Mobject",
          function(theObject,createExperimentName = "assay_reprocess"){
            objlist1 <- list(createExperimentName = theObject@assay_reprocess,
                             assay_raw = SummarizedExperiment::SummarizedExperiment(
                               theObject@assay_raw,rowData=theObject@row_data))
            names(objlist1) <- c(createExperimentName,"assay_raw")
            assay_reprocess_map <- data.frame(assay = rep(createExperimentName,ncol(theObject@assay_reprocess)),
                                              primary = colnames(theObject@assay_reprocess),
                                              colname = colnames(theObject@assay_reprocess),
                                              stringsAsFactors = FALSE)
            assay_raw_map <- data.frame(assay = rep("assay_raw",ncol(theObject@assay_raw)),
                                        primary = colnames(theObject@assay_raw),
                                        colname = colnames(theObject@assay_raw),
                                        stringsAsFactors = FALSE)
            dfmap1 <- rbind(assay_reprocess_map,assay_raw_map)

            results <- MultiAssayExperiment::MultiAssayExperiment(
              objlist1,theObject@primary, dfmap1,
              metadata = list(theObject@meta_data))
            return(results)
          }
)

#' @title Combine data into SummarizedExperiment/MultiAssayExperiment object.
#' @description \code{SCAN_reprocess_TRUE()} combines individual data into SummarizedExperiment
#' object for microarray study and into MultiAssayExperiment object for RNA-seq study.
#' This function will include reprocessed RNA-seq information.
#' @param geo_access A character/vector that contains geo accession number. If All, get all avaible studies.
#' @param include.SCAN A logical value indicates whether include normalized data processed by SCAN into the final output when available.
#' The default in FALSE.
#' @return A SummarizedExperiment object for microarray data, or a MultiAssayEpxeriment object
#' for RNA-seq data.
#' @examples
#' SCAN_reprocess_TRUE(geo_access = "GSE39939", include.SCAN = TRUE)
#' SCAN_reprocess_TRUE(geo_access = "GSE101705", include.SCAN = TRUE)
#' @export
SCAN_reprocess_TRUE <- function(geo_access, include.SCAN){
  param <- BiocParallel::SerialParam(progressbar=TRUE)
  if(geo_access[1] == "All"){

    # Get all available studies
    file_names_full <- utils::data(package="curatedTBData")[["results"]][,"Item"]
    file_names_full <- file_names_full[grep("GSE",file_names_full)]

    geo_access <- unique(gsub("_.*","",file_names_full))
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    objects_list <- BiocParallel::bplapply(seq_len(length(geo_index_list)), function(x){

    # Load Data into the Environment
    data_load <-  utils::data(list=file_names_full[geo_index_list[[x]]])
    data_list <- lapply(data_load, function(y) get(y))

    names(data_list) <- sub("^[^_]*_", "", data_load)

    # Remove data from environment
    objs <- ls(pos = ".GlobalEnv")
    rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

    # Check whether assemble into Summarized or MultiAssayExperiment Object
    # If no reporcess, then goes to SummarizedExperiment

    check_type <- grep("reprocess",data_load)
    if(length(check_type) == 0){ # combine into SummarizedExperiment

        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                        row_data = data_list$row_data,
                        column_data  = data_list$column_data,
                        meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assay_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assay_name]] <- data_list$SCAN_counts
        }
        return(sobject1_final)

      }

      else { # Combine into MultiAssayExperiment

        mobject1 <- methods::new("Mobject", assay_reprocess = as.matrix(data_list$assay_reprocess),
                        assay_raw = as.matrix(data_list$assay_raw_counts),
                        row_data = data_list$row_data,
                        primary = data_list$column_data,
                        meta_data = data_list$meta_data)

        # No need to check for SCAN in the RNA-seq. RNA-seq data does not have SCAN-normalized data
        mobject1_final <- CreateObject(mobject1)
        sobject1_final <- mobject1_final[["assay_raw"]]

        # Provide assay names
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        mobject1_final[["assay_raw"]] <- sobject1_final

        return(mobject1_final)

      }

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)
    return(objects_list)

  }
  else{

    file_names_full <- utils::data(package="curatedTBData")[["results"]][,"Item"]
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    # Check whether geo accession is available in the pakcage

    index <-  which(unlist(lapply(geo_index_list, length)) == 0)

    if(length(index)!=0){
      message(paste0(names(geo_index_list)[index])," is/are unavailable in the package")
      geo_index_list[index] <- NULL
    }

    if(length(geo_index_list)==0){stop("No available data found in the paackage")}

    objects_list <- BiocParallel::bplapply(seq_len(length(geo_index_list)), function(x){

      # Load Data into the Environment
      data_load <-  utils::data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- sub("^[^_]*_", "", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized/MultiAssayExperiment Object

      check_type <- grep("reprocess",data_load) # If no reporcess, then goes to SummarizedExperiment

      if(length(check_type) == 0){ # combine into SummarizedExperiment

        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                        row_data = data_list$row_data,
                        column_data = data_list$column_data,
                        meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assay_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assay_name]] <- data_list$SCAN_counts
        }

        return(sobject1_final)

      }

      else { # Combine into MultiAssayExperiment

        mobject1 <- methods::new("Mobject", assay_reprocess = as.matrix(data_list$assay_reprocess),
                        assay_raw = as.matrix(data_list$assay_raw_counts), row_data = data_list$row_data,
                        primary = data_list$column_data,meta_data = data_list$meta_data)

        mobject1_final <- CreateObject(mobject1, createExperimentName = "assay_reprocess")

        sobject1_final <- mobject1_final[["assay_raw"]]

        # Provide assay names
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        mobject1_final[["assay_raw"]] <- sobject1_final

        return(mobject1_final)

      }

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)

    return(objects_list)
  }
}

#' @title Combine individual data into SummarizedExperiment object
#' @description \code{SCAN_reprocess_FALSE()} combines individual data into SummarizedExperiment
#' object for all available study.
#' This function will not include reprocessed RNA-seq information.
#' @param geo_access A character/vector that contains geo accession number. If All, get all avaible studies.
#' @param include.SCAN A logical value indicates whether include normalized data processed by SCAN into the final output.
#' The default in FALSE.
#' @examples
#' SCAN_reprocess_TRUE(geo_access = "GSE39939", include.SCAN = TRUE)
#' SCAN_reprocess_TRUE(geo_access = "GSE101705", include.SCAN = TRUE)
#' @export
SCAN_reprocess_FALSE <- function(geo_access, include.SCAN){
  param <- BiocParallel::SerialParam(progressbar=TRUE)

  if(geo_access[1] == "All"){
    # Get all available studies
    file_names_full <- utils::data(package="curatedTBData")[["results"]][,"Item"]
    file_names_full <- file_names_full[grep("GSE",file_names_full)]

    geo_access <- unique(gsub("_.*","",file_names_full))
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    objects_list <- BiocParallel::bplapply(seq_len(length(geo_index_list)), function(x){

      # Load Data into the Environment
      data_load <-  utils::data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- sub("^[^_]*_", "", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized or MultiAssayExperiment Object
      # If no reporcess, goes to SummarizedExperiment directly


      # Check whether all samples from column data are included in the count matrix e.g. GSE79362
      if (suppressWarnings(all(colnames(data_list$assay_raw_counts) ==
                               row.names(data_list$column_data)))){
        # If sample names in count matrix are the same as sample names in the column data
        # combine individual study into SummarizedExperiment object
        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                                 row_data = data_list$row_data,
                                 column_data  = data_list$column_data,
                                 meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)
      }
      else{

        index <- match(colnames(data_list$assay_raw_counts),
                       row.names(data_list$column_data))
        data_list$column_data <- data_list$column_data[index,]
        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                                 row_data = data_list$row_data,
                                 column_data  = data_list$column_data,
                                 meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)
      }

      # Give assay name in SummarizedExperiment Object
      names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
      if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assay_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assay_name]] <- data_list$SCAN_counts
        }

      return(sobject1_final)

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)
    return(objects_list)
  }
  else{

    file_names_full <- utils::data(package="curatedTBData")[["results"]][,"Item"]
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    # Check whether geo accession is available in the pakcage

    index <-  which(unlist(lapply(geo_index_list, length)) == 0)

    if(length(index)!=0){
      message(paste0(names(geo_index_list)[index])," is/are unavailable in the package")
      geo_index_list[index] <- NULL
    }

    if(length(geo_index_list)==0){stop("No available data found in the paackage")}

    objects_list <- BiocParallel::bplapply(seq_len(length(geo_index_list)), function(x){

      # Load Data into the Environment
      data_load <-  utils::data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- sub("^[^_]*_", "", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")


      if (suppressWarnings(all(colnames(data_list$assay_raw_counts) ==
                               row.names(data_list$column_data)))){
        # If sample names in count matrix are the same as sample names in the column data
        # combine individual study into SummarizedExperiment object
        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                                 row_data = data_list$row_data,
                                 column_data  = data_list$column_data,
                                 meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)
      }
      else{

        index <- match(colnames(data_list$assay_raw_counts),
                       row.names(data_list$column_data))
        data_list$column_data <- data_list$column_data[index,]
        sobject1 <- methods::new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                                 row_data = data_list$row_data,
                                 column_data  = data_list$column_data,
                                 meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)
      }

        # Give assay name in SummarizedExperiment Object
      names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")

      if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assay_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assay_name]] <- data_list$SCAN_counts
        }

        return(sobject1_final)

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)

    return(objects_list)
  }
}

#' @title Load curatedTBData into SummarizedExperiment/MultiAssayExperiment object.
#' @description \code{get_curatedTBData} loads available curatedTBData.
#' @name get_curatedTBData
#' @param geo_access A character/vector that contains geo accession number. If All, get all avaible studies.
#' @param include.SCAN A logical value indicates whether include normalized data processed by SCAN into the final output.
#' The default in FALSE. SCAN-normalized results are available for
#' Affymetrix Microarray: GSE31348, GSE36238, GSE41055, GSE54992, GSE73408
#' Agilent two-color Microarray: GSE25534
#' @param include.reprocess A logical value indicated whether include reprocessed RNA-seq data when available.
#' The default is TRUE. Reprocess results are available for
#' Illumina RNA-seq: GSE94438, GSE79362, GSE89403, GSE107991, GSE107992, GSE107993,
#' GSE107994, GSE101705, GSE107104, GSE112104
#' @return A list of SummarizedExperiment/MultiAssayExperiment objects
#' @examples
#' get_curatedTBData(geo_access = "GSE39939")
#' get_curatedTBData(geo_access = c("GSE39939","GSE107993"))
#' get_curatedTBData(geo_access = "All")
#' @export
get_curatedTBData <- function(geo_access, include.SCAN=FALSE, include.reprocess=TRUE){

  # Include reprocessed RNA-seq
  if(include.reprocess){
    final_object <- SCAN_reprocess_TRUE(geo_access, include.SCAN)
    return(final_object)

  }

  # Not include reprocessed RNA-seq
  # All results are in the form of SummarizedExperiment object
  if(!include.reprocess){
    final_object <- SCAN_reprocess_FALSE(geo_access, include.SCAN)
    return(final_object)
  }

}

