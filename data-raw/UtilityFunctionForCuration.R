
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

relabel_TB <- function(dat_new){
  dat_new$CultureStatus <- sub(".*\\((.*)\\).*", "\\1", dat_new$TBStatus)

  dat_new$CultureStatus <- ifelse(dat_new$CultureStatus==dat_new$TBStatus,
                                  NA, dat_new$CultureStatus)


  TBStatus <- TBStatus_temp <- dat_new$TBStatus

  for (i in 1:length(TBStatus)){
    if (TBStatus[i] == unique(TBStatus_temp)[1] || TBStatus[i] == unique(TBStatus_temp)[4]){
      TBStatus[i] = 'PTB'
    }
    if (TBStatus[i] == unique(TBStatus_temp)[2] || TBStatus[i] == unique(TBStatus_temp)[5]){
      TBStatus[i] = 'OD'
    }
    if (TBStatus[i] == unique(TBStatus_temp)[3]){
      TBStatus[i] = 'LTBI'
    }
  }
  dat_new$TBStatus <- TBStatus
  return(S4Vecotrs::DataFrame(dat_new))
}

match_gene_symbol <- function(row_data){
  Symbol_R <- row_data$SYMBOL
  Symbol_plat <- row_data$Symbol

  # Use platform annotation as reference
  Symbol_plat_new <- Symbol_plat
  for (i in 1:length(Symbol_R)){
    if (Symbol_plat_new[i] == ""){
      Symbol_plat_new[i] = Symbol_R[i]
    }
  }
  row_data$SYMBOL_NEW <- Symbol_plat_new

  return(row_data)
}

save_raw_files <- function(sobject, path, geo) {
  column_data <- SummarizedExperiment::colData(sobject)
  row_data <- SummarizedExperiment::rowData(sobject)
  assay_raw_counts <- SummarizedExperiment::assay(sobject)
  meta_data <- S4Vectors::metadata(sobject)[[1]]
  lst <- list(column_data = column_data, row_data = row_data,
              assay_raw_counts = assay_raw_counts,
              meta_data = meta_data)
  lapply(1:length(lst), function(i) {
    saveRDS(lst[[i]], paste0(path, geo, "_", names(lst)[i], ".RDS"))
  })
}





