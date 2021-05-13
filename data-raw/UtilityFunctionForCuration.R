readRawData <- function(geo, sequencePlatform, urlIndex = NULL, matrix.only = FALSE) {
  urls <- GEOquery::getGEOSuppFiles(geo, fetch_files = FALSE)
  if (sequencePlatform == "GPL6947") {
    # Illumina Microarry V3
    result <- readRawDataGPL6947(urls, urlIndex)
    return(result)
  } else if (sequencePlatform == "GPL6102") {
    result <- readRawDataGPL6102(urls, urlIndex)
    return(result)
  } else if (sequencePlatform == "GPL10558") {
    # Illumina Microarry V4
    result <- readRawDataGPL10558(urls, urlIndex, matrix.only)
    return(result)
  } else if (sequencePlatform == "GPL570") {
    result <- readRawDataGPL570(urls, urlIndex)
    return(result)
  } else if (sequencePlatform == "GPL6883") {
    result <- readRawDataGPL6883(urls, urlIndex)
    return(result)
  }
}

readRawDataGPL6947 <- function(urls, urlIndex = NULL) {
  temp <- tempfile()
  tempd <- tempdir()
  url_sub <- as.character(urls$url[1])
  utils::download.file(url_sub, temp)
  utils::untar(temp, exdir = tempd)
  files <- list.files(tempd, pattern = "txt.*")
  data_Non_normalized_list <- lapply(files, function(x)
    read.delim(paste0(tempd, "/" ,x), header = TRUE,
               col.names = c("ID_REF", gsub("_.*", "", x),
                             paste0(gsub("_.*", "", x), ".Pval")),
               stringsAsFactors = FALSE))
  # Remove temporary files
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  # data_Non_normalized_list_noPvalue <- lapply(data_Non_normalized_list, function(x)
  #   x[, -grep('pval', colnames(x), ignore.case = TRUE)])
  message("Merge list into one matrix based on probe ID")
  data_Non_normalized <- Reduce(function (x, y)
    merge(x, y, by = "ID_REF", all = FALSE),
    lapply(data_Non_normalized_list, function(x) {x}))
  row.names(data_Non_normalized) <- data_Non_normalized$ID_REF
  data_Non_normalized <- as.matrix(data_Non_normalized[, -1])
  indexPvalue <- grep("pval", colnames(data_Non_normalized), ignore.case=TRUE)
  xr <- new("EListRaw", list(E = data_Non_normalized[, -indexPvalue],
                             other = list(Detection = data_Non_normalized[, indexPvalue])))
  yr <- limma::neqc(xr)
  return(list(data_Non_normalized = xr$E, data_normalized = yr$E))
}

readRawDataGPL6102 <- function(urls, urlIndex = NULL) {
  temp <- tempfile()
  tempd <- tempdir()
  url_sub <- as.character(urls$url[1])
  utils::download.file(url_sub, temp)
  utils::untar(temp, exdir = tempd)
  files <- list.files(tempd, pattern = "txt.*")
  data_Non_normalized_list <- lapply(files, function(x)
    read.delim(paste0(tempd, "/", x), header = TRUE,
               col.names = c("ID_REF", gsub("_.*", "", x),
                             paste0(gsub("_.*", "", x), ".Pval")),
               stringsAsFactors = FALSE))
  # Remove temporary files
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  data_Non_normalized_list_noPvalue <- lapply(data_Non_normalized_list, function(x)
    x[, -grep('pval', colnames(x), ignore.case = TRUE)])
  return(data_Non_normalized_list_noPvalue)
}


readRawDataGPL10558 <- function(urls, urlIndex = NULL, matrix.only = FALSE) {
  temp <- tempfile()
  # tempd <- tempdir()
  if (is.null(urlIndex)) {
    index <- grep("*.txt.gz", unlist(urls$url))
    if (length(index) == 0) {
      stop("Raw data with txt.gz not found from supplementary file. Exit the function.")
    }
    if (length(index) != 1) {
      message("More than one link selected, looking for non-nomarlized file")
      url_temp <- urls$url[index]
      index1 <- grep("non-normalized", url_temp)
      if (length(index1) != 1) {
        stop("Cannot identify the right urls. Plases specifcy it manually using the
             urlIndex parameter.")
      }
      url_sub <- as.character(url_temp[index1])
    } else {
      url_sub <- as.character(urls$url[index])
    }
  } else {
    url_sub <- as.character(urls$url[urlIndex])
  }
  # Illumina Microarray V4
  utils::download.file(url_sub, temp)
  data_Non_normalized <- read.delim(gzfile(temp), row.names = 1, header = TRUE) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0)  # delete columns that contain ONLY NAs
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  if (matrix.only) {
    return(data_Non_normalized)
  }
  indexPvalue <- grep("pval", colnames(data_Non_normalized), ignore.case=TRUE)
  if(length(indexPvalue) == 0) {
    message("Column(s) with p-value not found, return full datasets")
    return(data_Non_normalized)
  } else {
    xr <- new("EListRaw", list(E = data_Non_normalized[, -indexPvalue],
                               other = list(Detection = data_Non_normalized[, indexPvalue])))
    yr <- limma::neqc(xr)
    return(list(data_Non_normalized = xr$E, data_normalized = yr$E))
  }
}
readRawDataGPL6883 <- function(urls, urlIndex = NULL) {
  temp <- tempfile()
  tempd <- tempdir()
  if (is.null(urlIndex)) {
    index <- grep("*.txt.gz", unlist(urls$url))
    if (length(index) == 0) {
      stop("Raw data with txt.gz not found from supplementary file. Exit the function.")
    }
    if (length(index) != 1) {
      message("More than one link selected, looking for non-nomarlized file")
      url_temp <- urls$url[index]
      index1 <- grep("non-normalized", url_temp)
      if (length(index1) != 1) {
        stop("Cannot identify the right urls. Plases specifcy it manually using the
             urlIndex parameter.")
      }
      url_sub <- as.character(url_temp[index1])
    } else {
      url_sub <- as.character(urls$url[index])
    }
  } else {
    url_sub <- as.character(urls$url[urlIndex])
  }
  utils::download.file(url_sub, temp)
  data_Non_normalized <- read.delim(gzfile(temp), row.names = 1, header = TRUE) %>%
    dplyr::select_if(~sum(!is.na(.)) > 0)  # delete columns that contain ONLY NAs
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  indexPvalue <- grep("pval", colnames(data_Non_normalized), ignore.case=TRUE)
  if(length(indexPvalue) == 0) {
    message("Column(s) with p-value not found, return full datasets")
    return(data_Non_normalized)
  } else {
    data_non_pvalue <- data_Non_normalized[, -indexPvalue]
    return(data_non_pvalue)
  }
}

readRawDataGPL570 <- function(urls, urlIndex = NULL) {
  temp <- tempfile()
  tempd <- tempdir()
  url_cel <- as.character(urls$url[1])
  utils::download.file(url_cel, temp)
  utils::untar(temp, exdir = tempd)
  celFiles <- list.files(path = tempd, pattern = "*.CEL", full.names = TRUE)
  dataAffy <- affy::ReadAffy(filenames = celFiles)
  data_normalized_rma <- Biobase::exprs(affy::rma(dataAffy))
  # colnames(data_normalized_rma) <- gsub("_.*", "", colnames(data_normalized_rma))
  return(data_normalized_rma)
}

################################################################################
# Functions to create column data
readRawColData <- function(gse) {
  data_characteristic <- lapply(1:length(GEOquery::GSMList(gse)), function(x)
    GEOquery::GSMList(gse)[[x]]@header$characteristics_ch1)
  characteristic_table <- sapply(1:length(data_characteristic[[1]]), function(x)
    sapply(data_characteristic, "[[",x))
  characteristic_data_frame <- sub("(.*?): ","",characteristic_table) %>%
    S4Vectors::DataFrame()
  row.names(characteristic_data_frame) <- names(GEOquery::GSMList(gse))
  return(characteristic_data_frame)
}
create_standard_coldata <- function(col_data) {
  standard_name_seq <- c("Age", "Gender", "Ethnicity", "TBStatus", "GeographicalRegion",
                         "BcgVaccinated", "BirthRegion", "TST", "Tissue", "HIVStatus",
                         "MeasurementTime", "PatientID", "PneumoniaStatus",
                         "exposure_latent", "index_case_disease_site",
                         "smear_of_index_case", "modal_x_ray_grade", "SputumSmearStatus", "sputum_culture",
                         "bal_smear", "bal_culture", "isolate_sensitivity")
    # colData(readRDS('~/Desktop/RA work/build_TB_data/Signatures/Berry/GSE19443/GSE19443_summarizedExperiemnt.RDS')) %>% data.frame() %>% colnames
  col_info <- col_data %>% data.frame()
  dat_NA_new <- matrix(c(rep(NA, length(standard_name_seq) * nrow(col_info))),
                       ncol = length(standard_name_seq),
                       byrow = T,dimnames = list(row.names(col_info),
                                            standard_name_seq)) %>% data.frame()
  # fill in with existing data
  overlap_name <- standard_name_seq[which(colnames(dat_NA_new) %in% colnames(col_info))]
  for (i in overlap_name){
    dat_NA_new[i] <- col_info[i]
  }

  # Append the rest of cloumns to the dataframe
  col_info_rest <- col_info %>% dplyr::select(-all_of(overlap_name))
  dat_final <- cbind(dat_NA_new, col_info_rest)
  return(dat_final)
}

# Functions to create row data
map_gene_symbol <- function(data_non_pvalue, sequencePlatform) {
  data_non_pvalue <- data.frame(data_non_pvalue)
  sequence_result <- GEOquery::getGEO(sequencePlatform, GSEMatrix = FALSE)
  sequence_result_dat <- sequence_result@dataTable@table
  PROBES <- row.names(data_non_pvalue)
  if (sequencePlatform == "GPL6102") {
    OUT <- AnnotationDbi::select(illuminaHumanv2.db::illuminaHumanv2.db, PROBES, "SYMBOL")
  } else if (sequencePlatform == "GPL6947" || sequencePlatform == "GPL6883") {
    OUT <- AnnotationDbi::select(illuminaHumanv3.db::illuminaHumanv3.db, PROBES, "SYMBOL")
  } else if (sequencePlatform == "GPL10558") {
    OUT <- AnnotationDbi::select(illuminaHumanv4.db::illuminaHumanv4.db, PROBES, "SYMBOL")
  } else if (sequencePlatform == "GPL570") {
    index <- which(colnames(sequence_result_dat) == "Gene Symbol")
    colnames(sequence_result_dat)[index] <- "Symbol"
    OUT <- AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db, PROBES, "SYMBOL")
  }
  OUT[is.na(OUT)] <- NA
  # Map ProbeID to Gene Symbol
  OUT_collapse <- OUT %>%
    dplyr::group_by(PROBEID) %>%
    dplyr::summarise(SYMBOL = paste(SYMBOL, collapse = "///"),
                     times = length(unlist(strsplit(SYMBOL, "///"))))
  data_non_pvalue$ID_REF <- row.names(data_non_pvalue)
  data_final <- data_non_pvalue %>%
    dplyr::left_join(OUT_collapse, by=c("ID_REF" = "PROBEID")) %>%
    dplyr::left_join(sequence_result_dat, by = c("ID_REF" = "ID"))
  # Create row data
  row_data <- data_final %>%
    dplyr::select(-grep("GSM", colnames(data_final))) %>%
    S4Vectors::DataFrame()
  return(row_data)
}

match_gene_symbol <- function(row_data) {
  Symbol_R <- row_data$SYMBOL
  # Use platform annotation as reference
  Symbol_plat_new <- row_data$Symbol
  for (i in 1:length(Symbol_R)) {
    if (Symbol_plat_new[i] == "" || is.na(Symbol_plat_new[i]) ||
        Symbol_plat_new[i] == "NA") {
      Symbol_plat_new[i] = Symbol_R[i]
    }
  }
  row_data$SYMBOL_NEW <- Symbol_plat_new
  return(row_data)
}

# Perform normalization on probe level for the raw data
norm_probeToGenes_Agilent <- function(EListRaw, FUN = median) {
  y <- limma::backgroundCorrect(EListRaw, method = "normexp")
  y <- limma::normalizeBetweenArrays(y, method = "quantile")
  # Filter control probes and probes with no symbol
  Control <- y$genes$ControlType==1L
  NoSymbol <- is.na(y$genes$GeneName)
  yfilt <- y[!Control & !NoSymbol, ]
  dataNorm <- data.frame(yfilt$E)
  dataNorm$SYMBOL <- yfilt$genes$GeneName
  exprs2 <- stats::aggregate(. ~ SYMBOL, data = dataNorm,
                             FUN = FUN, na.action = na.pass)
  row.names(exprs2) <- exprs2$SYMBOL

  final <- as.matrix(exprs2[, -which(colnames(exprs2) == "SYMBOL")])
  return(final)
}
#normalizeIllumina
#normalizeHiseq
normalizeExprs <- function(data_Non_normalized, dataType, platform = NULL,
                          method = NULL) {
  if (dataType == "Microarray") {
    # data_Non_normalized[data_Non_normalized < 10] <- 10
    # datLog <- log(data_Non_normalized, base = 2) # log2 transformed data
    if (platform == "Agilent") {
      datBackground <- limma::backgroundCorrect.matrix(data_Non_normalized,
                                                       method = "normexp")
      datNormed <- limma::normalizeBetweenArrays(datBackground, method = method)
      return(datNormed)
    }
  } else if (dataType == "RNA-seq") {
    return(NULL)
  }
}
probesetsToGenes <- function(row_data, data_normalized, FUN) {
  row_data_sub <- row_data[which(row_data$ID_REF %in% row.names(data_normalized)),]
  if (!all(row.names(data_normalized) == row_data_sub$ID_REF)){
    stop("Row names are not exactly the same with ID_REF from row Data, consider change")
  }
  exprs1 <- data_normalized %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(SYMBOL = row_data_sub$SYMBOL_NEW) %>%
    dplyr::filter(.data$SYMBOL != "NA") %>%
    dplyr::filter(.data$SYMBOL != "")
  if(length(grep("///", exprs1$SYMBOL)) != 0){
    exprs1 <- expandProbesets(exprs1, sep = "///")
  }
  exprs2 <- stats::aggregate(. ~ SYMBOL, data = exprs1,
                             FUN = FUN, na.action = na.pass)
  row.names(exprs2) <- exprs2$SYMBOL

  final <- as.matrix(exprs2[, -which(colnames(exprs2) == "SYMBOL")])
  return(final)
}
# Match gene symbol to normalized data
makeCuratedExprs <- function(row_data, data_Non_normalized, dataType,
                             platform = NULL, method = NULL, FUN = median) {
  data_normalized <- normalizeExprs(data_Non_normalized, dataType, platform, method)
  dat_curated <- probesetsToGenes(row_data, data_normalized, FUN)
  return(dat_curated)
}

# Save files in separate RDS file
save_raw_files <- function(sobject, path, geo) {
  column_data <- SummarizedExperiment::colData(sobject)
  row_data <- SummarizedExperiment::rowData(sobject)
  assay_raw <- SummarizedExperiment::assay(sobject)
  meta_data <- S4Vectors::metadata(sobject)[[1]]
  lst <- list(column_data = column_data, row_data = row_data,
              assay_raw = assay_raw,
              meta_data = meta_data)
  lapply(1:length(lst), function(i) {
    saveRDS(lst[[i]], paste0(path, geo, "_", names(lst)[i], ".RDS"))
  })
}

expandProbesets <- function(dat, sep){
  times <- NULL
  # Get index with duplicated symbol
  index <- grep(sep, dat$SYMBOL)
  sobject_exprs_dup <- dat[index,]

  x_list <- strsplit(as.character(sobject_exprs_dup$SYMBOL), sep)
  symbol_dup <- gsub(" ", "", unlist(x_list))
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




