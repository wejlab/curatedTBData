#' Download data from Illumina HumanHT-12 V4.0 expression beadchip.
#' @name download_data_Affy
#'
#' @param geo_access A character that contains GEO accession number.
#' @return A dataframe that contains non-normalized read counts.
#' @examples
#' download_data_Affy("GSE36238")
#' @export
download_data_Affy<- function(geo_access){
  urls <- GEOquery::getGEOSuppFiles(geo_access, fetch_files = FALSE)
  url_temp <- as.character(urls$url[grep(".tar",urls$url)])
  if(!is.null(url_temp)){
    temp <- tempfile()
    tempd <- tempdir()
    download.file(url_temp,temp)
    untar(temp,exdir = tempd)
    celFiles <- list.files(path = tempd, pattern = '*.CEL',full.names=TRUE)
    # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
    # Conducr RMA
    data_affy <- affy::ReadAffy(filenames = celFiles)
    normalized_rma <- exprs(affy::rma(data_affy))
    return(normalized_rma)

  }
  else{stop(paste0("No valid url for ",geo_access))}

}
#####################################


