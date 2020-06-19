#########################################################
#' Combine individual data to SummarizedExperiment/MultiAssayExperiment object
#' when include.reprocess = TRUE
#' @param geo_access A character/vector that contains geo accession number. If All, get all avaible studies.
#' @param include.SCAN A logical value indicates whether include normalized data processed by SCAN into the final output when available.
#' The default in FALSE.
#' @export
SCAN_reprocess_TRUE <- function(geo_access, include.SCAN){
  param <- BiocParallel::SerialParam(progressbar=TRUE)
  if(geo_access[1] == "All"){

    # Get all available studies
    file_names_full <- data(package="curatedTBData")[["results"]][,"Item"]
    file_names_full <- file_names_full[grep("GSE",file_names_full)]

    geo_access <- unique(gsub("_.*","",file_names_full))
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    objects_list <- BiocParallel::bplapply(1:length(geo_index_list), function(x){

      # Load Data into the Environment
      data_load <-  data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- gsub(paste0(".*",names(geo_index_list)[x],"_",
                                      "([^.]+)[.].*"),"\\1", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized or MultiAssayExperiment Object
      # If no reporcess, then goes to SummarizedExperiment

      check_type <- grep("reprocess",data_load)
      if(length(check_type) == 0){ # combine into SummarizedExperiment

        sobject1 <- new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                        row_data = data_list$row_data,
                        column_data  = data_list$column_data,
                        meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assy_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assy_name]] <- data_list$SCAN_counts
        }
        return(sobject1_final)

      }

      else { # Combine into MultiAssayExperiment

        mobject1 <- new("Mobject", assay_reprocess = as.matrix(data_list$assay_reprocess),
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

    file_names_full <- data(package="curatedTBData")[["results"]][,"Item"]
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    # Check whether geo accession is available in the pakcage

    index <-  which(unlist(lapply(geo_index_list, length)) == 0)

    if(length(index)!=0){
      message(paste0(names(geo_index_list)[index])," is/are unavailable in the package")
      geo_index_list[index] <- NULL
    }

    if(length(geo_index_list)==0){stop("No available data found in the paackage")}

    objects_list <- BiocParallel::bplapply(1:length(geo_index_list), function(x){

      # Load Data into the Environment
      data_load <-  data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- gsub(paste0(".*",names(geo_index_list)[x],"_","([^.]+)[.].*"),"\\1", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized/MultiAssayExperiment Object

      check_type <- grep("reprocess",data_load) # If no reporcess, then goes to SummarizedExperiment

      if(length(check_type) == 0){ # combine into SummarizedExperiment

        sobject1 <- new("Sobject", assay = as.matrix(data_list$assay_raw_counts), row_data = data_list$row_data,
                        column_data  = data_list$column_data, meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assy_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assy_name]] <- data_list$SCAN_counts
        }

        return(sobject1_final)

      }

      else { # Combine into MultiAssayExperiment

        mobject1 <- new("Mobject", assay_reprocess = as.matrix(data_list$assay_reprocess),
                        assay_raw = as.matrix(data_list$assay_raw_counts), row_data = data_list$row_data,
                        primary = data_list$column_data,meta_data = data_list$meta_data)

        mobject1_final <- CreateObject(mobject1)

        return(mobject1_final)

      }

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)

    return(objects_list)
  }
}

#' Combine individual data to SummarizedExperiment/MultiAssayExperiment object
#' when include.reprocess = FALSE
#' @param geo_access A character/vector that contains geo accession number. If All, get all avaible studies.
#' @param include.SCAN A logical value indicates whether include normalized data processed by SCAN into the final output.
#' The default in FALSE.
#' @export
SCAN_reprocess_FALSE <- function(geo_access, include.SCAN){
  param <- BiocParallel::SerialParam(progressbar=TRUE)
  if(geo_access[1] == "All"){
    # Get all available studies
    file_names_full <- data(package="curatedTBData")[["results"]][,"Item"]
    file_names_full <- file_names_full[grep("GSE",file_names_full)]

    geo_access <- unique(gsub("_.*","",file_names_full))
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    objects_list <- BiocParallel::bplapply(1:length(geo_index_list), function(x){

      # Load Data into the Environment
      data_load <-  data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- gsub(paste0(".*",names(geo_index_list)[x],"_","([^.]+)[.].*"),
                               "\\1", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized or MultiAssayExperiment Object
      # If no reporcess, then goes to SummarizedExperiment

      # check_type <- grep("reprocess",data_load)
      # combine studies into SummarizedExperiment object

        sobject1 <- new("Sobject", assay = as.matrix(data_list$assay_raw_counts),
                        row_data = data_list$row_data,
                        column_data  = data_list$column_data,
                        meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")
        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assy_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assy_name]] <- data_list$SCAN_counts
        }

        return(sobject1_final)

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)
    return(objects_list)
  }
  else{

    file_names_full <- data(package="curatedTBData")[["results"]][,"Item"]
    geo_index_list <- lapply(geo_access, function(x) grep(x,file_names_full))
    names(geo_index_list) <- geo_access

    # Check whether geo accession is available in the pakcage

    index <-  which(unlist(lapply(geo_index_list, length)) == 0)

    if(length(index)!=0){
      message(paste0(names(geo_index_list)[index])," is/are unavailable in the package")
      geo_index_list[index] <- NULL
    }

    if(length(geo_index_list)==0){stop("No available data found in the paackage")}

    objects_list <- BiocParallel::bplapply(1:length(geo_index_list), function(x){

      # Load Data into the Environment
      data_load <-  data(list=file_names_full[geo_index_list[[x]]])
      data_list <- lapply(data_load, function(y) get(y))

      names(data_list) <- gsub(paste0(".*",names(geo_index_list)[x],"_","([^.]+)[.].*"),
                               "\\1", data_load)

      # Remove data from environment
      objs <- ls(pos = ".GlobalEnv")
      rm(list = objs[grep(names(geo_index_list)[x], objs)], pos = ".GlobalEnv")

      # Check whether assemble into Summarized/MultiAssayExperiment Object

      # check_type <- grep("reprocess",data_load) # If no reporcess, then goes to SummarizedExperiment

      # combine studues into SummarizedExperiment Object

        sobject1 <- new("Sobject", assay = as.matrix(data_list$assay_raw_counts), row_data = data_list$row_data,
                        column_data  = data_list$column_data, meta_data = data_list$meta_data)

        sobject1_final <- CreateObject(sobject1)

        # Give assay name in SummarizedExperiment Object
        names(SummarizedExperiment::assays(sobject1_final)) <- paste0(names(geo_index_list)[x],
                                                                      "_raw")

        if(include.SCAN && !is.null(data_list$SCAN_counts)){
          assy_name <- paste0(names(geo_index_list)[x],"_SCAN")
          SummarizedExperiment::assays(sobject1_final)[[assy_name]] <- data_list$SCAN_counts
        }

        return(sobject1_final)

    }, BPPARAM = param)

    names(objects_list) <- names(geo_index_list)

    return(objects_list)
  }
}

#' Combine individual data to SummarizedExperiment/MultiAssayExperiment object
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
#' get_curatedTBData("GSE39939")
#' get_curatedTBData(c("GSE39939","GSE107993"))
#' get_curatedTBData("All")
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

