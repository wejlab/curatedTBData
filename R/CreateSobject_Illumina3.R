#' Download data from Illumina HumanHT-12 V3.0 expression beadchip.
#' @name download_data_Illumina3
#'
#' @param geo_access A character that contains GEO accession number.
#' @return A list that contains non-normalized read counts.
#' @examples
#' download_data_Illumina3("GSE19442")
#' @export
download_data_Illumina3<- function(geo_access){
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
 # if (! "GEOquery" %in% installed.packages()) BiocManager::install("GEOquery")
  #if (! "dplyr" %in% installed.packages()) install.packages("dplyr")
  #library(dplyr)
  #library(GEOquery)
  urls <- GEOquery::getGEOSuppFiles(geo_access, fetch_files = FALSE)
  url_temp <- as.character(urls$url[grep(".tar",urls$url)])
  if(!is.null(url_temp)){
    temp <- tempfile()
    tempd <- tempdir()
    download.file(url_temp,temp)
    untar(temp,exdir = tempd)
    files <- list.files(tempd, pattern = "txt.*")
    # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
    Non_normalized_list <- lapply(files, function(x)
      read.delim(paste0(tempd,'/',x), header=TRUE, col.names = c('ID_REF',gsub('_.*','',x),paste0(gsub('_.*','',x),'.Pval')), stringsAsFactors = FALSE))
    return(Non_normalized_list)

  }
  else{stop(paste0("No valid url for ",geo_access))}

}
##################################################
#' Create non-normalized data from Illumina HumanHT-12 V3.0 expression beadchip
#' @name create_rawdata_Illumina3
#' @param Non_normalized_list A list that contains non_normalized reads from each sample.
#' @return A dataframe that contains non-normalized read counts.
#' @examples
#' Non_normalized_list <- download_data_Illumina3("GSE19442")
#' Non_normalized_data <- create_data_Illumina3(Non_normalized_lists)
#' @export
create_data_Illumina3 <- function(Non_normalized_list){
  Non_normalized_list_ID <- sapply(Non_normalized_list, function(x) x[,1]) %>% data.frame()
  check_1 <- sapply(1:ncol(Non_normalized_list_ID),function(x) all(Non_normalized_list_ID[,1]==Non_normalized_list_ID[,x]))
  if (sum(check_1 = F)==0){
    # Combine list into data frame
    Non_normalized <- do.call(cbind,Non_normalized_list)
    # Remove redundant ID_REF
    Non_normalized <- Non_normalized[,-which(colnames(Non_normalized)=='ID_REF')[-1]]
    row.names(Non_normalized) <- Non_normalized$ID_REF

    # Remove .Pval
    Non_pvalue <- Non_normalized[,-c(1,grep('.Pval',colnames(Non_normalized)))]
    return(Non_pvalue)

  }
  else{stop(paste0("Probe ID from each data is not the same, please check"))}
}
##################################################





