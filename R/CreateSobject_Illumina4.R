#' Download data from Illumina HumanHT-12 V4.0 expression beadchip.
#' @name download_data_Illumina4
#'
#' @param geo_access A character that contains GEO accession number.
#' @return A dataframe that contains non-normalized read counts.
#' @examples
#' download_data_Illumina4("GSE39939")
#' @export
download_data_Illumina4<- function(geo_access){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (! "GEOquery" %in% installed.packages()) BiocManager::install("GEOquery")
  if (! "dplyr" %in% installed.packages()) install.packages("dplyr")
  library(dplyr)
  library(GEOquery)
  urls <- GEOquery::getGEOSuppFiles(geo_access, fetch_files = FALSE)
  url_non_normalized <- as.character(urls$url[grep(".txt",urls$url)])
  if(!is.null(url_non_normalized)){
    temp <- tempfile()
    download.file(url_non_normalized,temp)
    Non_normalized <- read.delim(gzfile(temp),row.names = 1, header = TRUE) %>% dplyr::select_if(~sum(!is.na(.)) > 0)

    # unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
    if (!is.null(grep('.Pval',colnames(Non_normalized)))){
      dat_Non_pvalue <- Non_normalized[,-grep('.Pval',colnames(Non_normalized))]
      return(dat_Non_pvalue)
    }
    else {return(Non_normalized)}

  }
  else{stop(paste0("No valid url for ",geo_access))}

}
#####################################

#' Create row data information from Illumina HumanHT-12 V3.0/4.0 expression beadchip and Affymetrix Human Genome U133 Plus 2.0 Array
#' @name create_RowData
#' @param dat_download A dataframe object read from `download_data_Illumina4` or `create_data_Illumina3`
#' @param plat_access A string that indicates platform access number
#'
#' @return A DataFrame object contains the information for each probe ID
#' @examples
#' Non_normalized_data <- download_data_Illumina4("GSE19442")
#' plat_access <- "GPL10558"
#' row_data <- create_RowData(Non_normalized_data,"GPL10558")
#' @export
create_RowData <- function(dat_download, plat_access){

  dat_download <- data.frame(dat_download)
  # Obatian sample names
  samplename <- colnames(dat_download)

  # Read in data from vendor's information
  HumanHT <- GEOquery::getGEO(plat_access, GSEMatrix = F)
  HumanHT_dat <- HumanHT@dataTable@table

  # Annotation
  PROBES <- row.names(dat_download)
  library(illuminaHumanv3.db)
  library(illuminaHumanv4.db)
  library(hgu133plus2.db)

  if (plat_access == "GPL6947"){
    OUT <- AnnotationDbi::select(illuminaHumanv3.db, PROBES, "SYMBOL")
    OUT[is.na(OUT)] <- NA
  }
  if (plat_access== "GPL10558"){
    OUT <- AnnotationDbi::select(illuminaHumanv4.db, PROBES, "SYMBOL")
    OUT[is.na(OUT)] <- NA
  }

  # Affymetrix hgu133plus2.db
  if (plat_access== "GPL570"){
    OUT <- AnnotationDbi::select(hgu133plus2.db, PROBES, "SYMBOL")
    OUT[is.na(OUT)] <- NA
    OUT_collapse <- OUT %>% dplyr::group_by(PROBEID) %>%
      dplyr::summarise(SYMBOL = paste(SYMBOL, collapse="///"),
                       times = length(unlist(strsplit(SYMBOL, "///"))))
    dat_download$ID_REF <- row.names(dat_download)
    dat_download_final <- dat_download %>% dplyr::left_join(OUT_collapse, by=c('ID_REF' = "PROBEID")) %>% dplyr::left_join(HumanHT_dat, by = c('ID_REF'='ID'))

    # Create Row data annotation
    row_data <- dat_download_final %>% dplyr::select(-samplename)
    Symbol_R <- row_data$SYMBOL
    Symbol_plat <- row_data$`Gene Symbol`

    # Use platform annotation as reference
    Symbol_plat_new <- Symbol_plat
    for (i in 1:length(Symbol_R)){
      if (Symbol_plat_new[i]==""){
        Symbol_plat_new[i] = Symbol_R[i]
      }
    }
    # Create new variable called SYMBOL_NEW, used later in creating multi-assay pbject
    row_data$SYMBOL_NEW <- Symbol_plat_new

    return(DataFrame(row_data))

  }

  # Continue part illumina 3.0/4.0
  # Map ProbeID to Gene Symbol
  OUT_collapse <- OUT %>% dplyr::group_by(PROBEID) %>%
    dplyr::summarise(SYMBOL = paste(SYMBOL, collapse="///"),
                     times = length(unlist(strsplit(SYMBOL, "///"))))
  dat_download$ID_REF <- row.names(dat_download)
  dat_download_final <- dat_download %>% dplyr::left_join(OUT_collapse, by=c('ID_REF' = "PROBEID")) %>% dplyr::left_join(HumanHT_dat, by = c('ID_REF'='ID'))

  # Create Row data annotation
  row_data <- dat_download_final %>% dplyr::select(-samplename)
  Symbol_R <- row_data$SYMBOL
  Symbol_plat <- row_data$Symbol

  # Use platform annotation as reference
  Symbol_plat_new <- Symbol_plat
  for (i in 1:length(Symbol_R)){
    if (Symbol_plat_new[i]==""){
      Symbol_plat_new[i] = Symbol_R[i]
    }
  }
  # Create new variable called SYMBOL_NEW, used later in creating multi-assay pbject
  row_data$SYMBOL_NEW <- Symbol_plat_new

  return(DataFrame(row_data))

}
###################################

#' Create Column data information from Illumina HumanHT-12 V3.0/4.0 expression beadchip and Affymetrix Human Genome U133 Plus 2.0 Array
#' @name create_ColData
#' @inheritParams download_data_Illumina4
#'
#' @return A DataFrame object contains the patient information for each sample
#'
#' @examples
#' gse <- GEOquery::getGEO(geo_access, GSEMatrix =  F)
#' create_ColData_Illumina4("GSE39939",gse)
#'
#' @export
create_ColData <- function(geo_access,gse){

  # Obatain sample characteristics
  data_characteristic <- lapply(1:length(GSMList(gse)), function(x)
    GSMList(gse)[[x]]@header$characteristics_ch1)
  characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
    sapply(data_characteristic, '[[',x))

  # get string after :
  characteristic_data_frame <- sub('(.*?): ','',characteristic_table) %>% as_tibble()
  colnames(characteristic_data_frame) <- sub(':.*','',characteristic_table)[1,]
  row.names(characteristic_data_frame) <- names(GSMList(gse))
  return(DataFrame(characteristic_data_frame))

}

##############################################
#' Make column information consitent across dataset
#' @name create_new_ColData
#' @inheritParams create_ColData
#'
#' @return A DataFrame object contains the patient information for each sample
#'
#' @examples
#' gse <- GEOquery::getGEO(geo_access, GSEMatrix =  F)
#' col_data <- create_ColData("GSE39939",gse)
#' col_data_new <- create_new_ColData(col_data)
#' @export
create_new_ColData <- function(col_data){
  ## Import standard name sequence
  standard_name_seq <- c("Age","Gender","Ethnicity","TBStatus","GeographicalRegion","BcgVaccinated",
                         "BirthRegion","TST","exposure_latent", "index_case_disease_site","smear_of_index_case",
                         "modal_x_ray_grade","SputumSmearStatus","sputum_culture","bal_smear",
                         "bal_culture","isolate_sensitivity")

  col_info <- col_data %>% data.frame()
  dat_NA_new <- matrix(c(rep('n.a.',length(standard_name_seq)*nrow(col_info))),ncol=length(standard_name_seq),
                       byrow=T,dimnames=list(row.names(col_info),
                                             standard_name_seq)) %>% data.frame()

  # fill in with existing data
  overlap_name <- standard_name_seq[which(colnames(dat_NA_new) %in% colnames(col_info))]
  for (i in overlap_name){
    dat_NA_new[i] <- col_info[i]
  }

  # Append the rest of cloumns to the dataframe
  col_info_rest <- col_info %>% dplyr::select(-overlap_name)
  dat_final <- cbind(dat_NA_new,col_info_rest)

  dat_final <- DataFrame(dat_final)

  return(dat_final)
}

#########################################################
#' Create Summarized Experiment object
#' @name create_sobject_Illumina4
#' @param data A datafrmae that contains read count, column names must be sample ID,
#' row name must be probe ID
#' @param column_data A DataFrame that contains column information for each sample
#' @param row_data A DataFrame that contains row information for each probe
#' @param meta_data Formal class 'MIAME' that contains information abou the GEO dataset
#'
#' @return A summarized experiment oject for the GEO dataset
#' @examples
#' data <- download_data_Illumina4("GSE39939")
#' row_data <- create_RowData_Illumina4("GSE39939")
#' column_data <- create_ColData_Illumina4("GSE39939")
#' meta_data <- new('MIAME', name="Victoria Wright", lab="Wright Fleming Institute",
#' contact="v.wright_at_imperial.ac.uk",
#' title="Diagnosis of childhood tuberculosis and host RNA expression in Africa")



