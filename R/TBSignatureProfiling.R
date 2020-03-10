#' Remove empty objects from list contains both SummariexExperiment and MultiAssayExpriment objects
#' @name remove_empty_object
#' @param multi_object A list contains both SummariexExperiment/MultiAssayExpriment objects
#' @return A list contains non-empty SummariexExperiment/MultiAssayExpriment object
#' @export
remove_empty_object <- function(k){
  # should assign empty list with NA NOT NULL, because assigning NULL to list items, removes them...
  x <- k
  for (i in which(sapply(x, function(x) class(x) == "SummarizedExperiment"))){
    if(nrow(colData(x[[i]]))==0){
      x[[i]] <- NA
    }
  }
  for (j in which(sapply(x, function(x) class(x) == "MultiAssayExperiment"))){
    if(length(experiments(x[[j]]))==0){
      x[[j]] <- NA
    }
  }
  x <- x[!is.na(x)]
  return(x)
}
#######################
#' Create sobjects for TB signature profiling based on specific TB status.
#' @name get_sobject_TBSig
#'
#' @param multi_object an multi-assay summarized experiment
#' @param disease1  A character chooses from "Control", "Latent", "PTB", "OD"
#' @param disease2  A character chooses from "Control", "Latent", "PTB", "OD".
#' TBStatus1 and TBStatus2 should be different from each other
#' @return A summarized experiment object for TB signature profiling
#'
#' @examples
#' sobject <- dat("GSE39939_sobject")
#' mobject <- MatchProbe(sobject)
#' sobject_TBSig <- get_sobject_TBSig(mobject,"PTB","Latent")
#' @export
get_sobject_TBSig <- function(multi_object,disease1,disease2,assay_type = "assay_reduce"){

  # check input type
  if (!all(c(disease1,disease2) %in% c("Control", "Latent", "PTB", "OD"))){
    stop(cat("Invalid disease type. Two disease must be diffrent and choose from", "\"Control\"", "\"Latent\"",
             "\"PTB\"", "\"OD\""))
  }
  if (disease1==disease2){
    stop(cat("Invalid disease type. Two diseases must be diffrent and choose from", "\"Control\"",
             "\"Latent\"", "\"PTB\"", "\"OD\""))
  }
  if (class(multi_object) != "MultiAssayExperiment" ){
    stop(cat("Invalid input data type. Only supported for MultiAssayExperiment objects generated from \"Create_MultiAssay_object()\".
        Your input:", class(test_sobject)))
  }
  col_info <- colData(multi_object)
  col_data <- data.frame(Sample=row.names(col_info) %>% as.factor(),
                         Disease = col_info$TBStatus %>% as.factor())
  row.names(col_data) <- row.names(col_info)

  # when not all samples are included in the expression matrix
  if (ncol(multi_object[[assay_type]]) != nrow(col_info)){
    index <- sapply(1:length(colnames(multi_object[[assay_type]])), function (i)
      which(row.names(col_data) %in% colnames(multi_object[[assay_type]])[i]))

    col_data <- col_data[index,]
  }

  sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(multi_object[[assay_type]])), colData = col_data)

  # subsetting disease1 and disease2
  sobject_TBSig_filter <- sobject_TBSig[,sobject_TBSig$Disease %in% c(disease1,disease2)]
  TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)$Disease
  # check if both status are in the column data
  if(length(unique(TB_status)) == 2){
    return(sobject_TBSig_filter)
  }

}

#####################################
#' A function to show boxplot across signatures from a list of SummarizedExperiment Object
#' @name Boxplot_TBSig
#' @param result_list A list of SummarizedExperiment Object from `runTBsigProfiler`
#' @param sig_name Name of the TB Signatures
#' @return Boxplot
#'
#' @export
Boxplot_TBSig <- function(result_list,sig_name){

  p_boxplot <- lapply(1:length(result_list), function(x)
    TBSignatureProfiler::signatureBoxplot(inputData = result_list[[x]],
                                          name = names(result_list)[x],
                                          signatureColNames = sig_name,
                                          annotationColName = "TBStatus", rotateLabels = FALSE,fill_colors = c("#4E84C4", "#FC4E07")))

  library(gridExtra)
  library(ggplot2)
  p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
  return(p_combine)

}

#######################################
#' Get signatures with only scores and disease status


#' Boxplot functions for list with inconsistent Signatures matched for each object
#' Also applies to consistent sigantures column names
Boxplot_TBSig1 <- function(result_list, gset, sig_name, annotationName = "TBStatus"){
  sig_list <- lapply(1:length(result_list), function(y,gset){

    x <- result_list[[y]]
    GSE <- rep(names(result_list[y]), nrow(colData(x)))
    TBStatus <- colData(x)[,"TBStatus"]
    index <- na.omit(match(names(gset),names(colData(x))))
    cbind(TBStatus,colData(x)[,index],GSE)

  }, gset)

  # Merge them into one large dataset
  sig_data <- plyr::rbind.fill(lapply(sig_list,function(x){as.data.frame(x)}))

  p_boxplot <- lapply(unique(sig_data$GSE), function(x,sig_name){
    sig_data_gse <- sig_data %>% filter(GSE == x)
    if(sig_data_gse %>% dplyr::select(sig_name) %>% is.na() %>% all()){return(NULL)}
    sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)


    if(length(unique(colData(sig_data1)$TBStatus)) == 2){
      p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                  name = x,
                                                  signatureColNames = sig_name,
                                                  annotationColName = annotationName,
                                                  rotateLabels = FALSE,fill_colors = c("#4E84C4", "#FC4E07"))
      return(p)
    }
    if(length(unique(colData(sig_data1)$TBStatus)) ==3){
      p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                  name = x,
                                                  signatureColNames = sig_name,
                                                  annotationColName = annotationName,
                                                  rotateLabels = FALSE,fill_colors = c("#4E84C4", "#FC4E07","#999999"))

      return(p)
    }

  },sig_name)

  p_boxplot <- plyr::compact(p_boxplot)

  library(gridExtra)
  library(ggplot2)
  p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
  return(p_combine)


}

###############################################

#' Obtain two-sample t-test pvalues and emprirical AUC for signature scores.
#' @name get_pvalue_auc
#' @param SE_scored A summarized experiment from TB signature profiling
#' @param annotationColName A character, column name for TB status (default is "Disease")
#' @param signatureColNames A character or vector that contains gene signature name
#' @return A data frame contains p-value from two-sample t-test and AUC value for each signature
#'
#' @export
get_pvalue_auc <- function(SE_scored, annotationColName = "TBStaus", signatureColNames){

  # check signatureColNames

  pvals <- aucs <- NULL
  annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][,1] %>% as.character() %>% as.factor()
  for (i in signatureColNames) {
    score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
    pvals <- c(pvals, stats::t.test(score ~ annotationData)$p.value)
    pred <- ROCit::rocit(score, annotationData)
    auc <- pred$AUC
    aucs <- c(aucs, max(auc, 1 - auc))
  }
  return(data.frame(Signature=signatureColNames,P.value=round(pvals,4),AUC=round(aucs,4)))

}

################################
#' Combine results from list. Calculate p-value and AUC values
#' @name combine_auc
#' @param result_list A list of SummarizedExperiment Object from `runTBsigProfiler`
#' @return A data frame with signatures, p-value, and AUC
#' @export
combine_auc <- function(result_list){
  aucs_result <- lapply(result_list, function(x){
    index <- na.omit(match(names(TBsignatures),names(colData(x))))
    get_pvalue_auc(x,
                   annotationColName = "TBStatus",
                   signatureColNames = names(colData(x))[index])
  }
    )
  aucs_result_dat <- do.call(rbind,aucs_result)

  # re-order data based on their median AUC
  aucs_result_dat_median <- aucs_result_dat %>% group_by(Signature) %>% summarise_all(median) %>% arrange(desc(AUC))

  # New addition: order signatures based on median AUC values

  Signature_order <- as.character(aucs_result_dat_median$Signature)
  # Re-order gene siganture, re-level
  aucs_result_dat$Signature <- factor(aucs_result_dat$Signature, levels = Signature_order)

  return(aucs_result_dat)
}

###########################

#' Obtain ridge plots for emprirical AUC distribution for signature scores.
#' @name get_auc_distribution
#' @param aucs_result_dat A dataframe contains AUC and p-value for certain TB signatures across different datasets
#' @return Ridge plot with median line
#'
#' @examples
#' aucs_result_dat <- data.frame(Signature=c("Anderson_42", "Anderson_OD_51", "Berry_393"), AUC=rnorm(3))
#' p_ridge <- get_auc_distribution(aucs_result_dat)
#' @export
get_auc_distribution <- function(aucs_result_dat){
  library(ggplot2)
  library(gridExtra)
  library(ggridges)

  # add 50% AUC line
  aucs_result_dat_lines <- data.frame(Signature = aucs_result_dat$Signature,x0=0.5)

  p_ridge <- ggplot(aucs_result_dat,aes(x=AUC,y=Signature)) + geom_density_ridges(jittered_points=TRUE,alpha=0.7,quantile_lines = TRUE, quantiles = 2) +
    geom_segment(data = aucs_result_dat_lines, aes(x = x0, xend = x0, y = as.numeric(Signature),
                                                   yend = as.numeric(Signature) + .9), color = "red")
  return(p_ridge)
}
