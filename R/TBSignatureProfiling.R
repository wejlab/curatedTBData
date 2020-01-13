#' Create sobjects for TB signature profiling based on specific TB status.
#' @name get_sobject_TBSig
#'
#' @param multi_object an multi-assay summarized experiment
#' @param disease1  A character chooses from "Control", "Latent", "PTB", "OD"
#' @param disease2  A character chooses from "Control", Latent, "PTB", "OD".
#' TBStatus1 and TBStatus2 should be different from each other
#' @return A summarized experiment object for TB signature profiling
#'
#' @examples
#' sobject <- dat("GSE39939_sobject")
#' mobject <- get_expr_symbol(sobject)
#' sobject_TBSig <- get_sobject_TBSig(mobject,"PTB","Latent")
#' @export
get_sobject_TBSig <- function(multi_object,disease1,disease2){

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
  col_info <- multi_object[[1]] %>% colData() %>% data.frame()
  col_data <- data.frame(Sample=row.names(col_info) %>% as.factor(),Disease = col_info$TBStatus %>% as.factor())
  row.names(col_data) <- row.names(col_info)

  sobject_TBSig <- SummarizedExperiment::SummarizedExperiment(assays = list(counts= as.matrix(multi_object[[2]])), colData = col_data)

  # subseeting disease1 and disease2
  sobject_TBSig_filter <- sobject_TBSig[,sobject_TBSig$Disease %in% c(disease1,disease2)]
  TB_status <- SummarizedExperiment::colData(sobject_TBSig_filter)$Disease
  # check if both status are in the column data
  if(length(unique(TB_status)) == 2){
    return(sobject_TBSig_filter)
  }

}
################################


########################
#' Get boxplot plot comparison across all signatures
#' @name get_boxplot
#' @param result_list A list that contains TB signature scores for each dataset
#' @param sig_name Character(s) stands for signature name(s)
#' @return Combined boxplot that shows specific signature across all dataset
#' @examples
#' p_boxplot <- get_boxplot(TB_sobject,"Sloot_HIV_2")
#' @export
get_boxplot <- function(result_list,sig_name){
  if (class(result_list) == "SummarizedExperiment"){
    p_plot <- TBSignaturePro::signatureBoxplot(inputData = result_list,
                                               name = names(result_list),
                                               signatureColNames = sig_name,
                                               annotationColName = "Disease", rotateLabels = FALSE,fill_colors = c("#00AFBB", "#E7B800"))
    return(p_plot)
  }
  if (class(result_list) == "list" & class(result_list[[1]]) == "SummarizedExperiment"){
    p_boxplot <- lapply(1:length(result_list), function(x)
      TBSignaturePro::signatureBoxplot(inputData = result_list[[x]],
                                       name = names(result_list)[x],
                                       signatureColNames = sig_name,
                                       annotationColName = "Disease", rotateLabels = FALSE,fill_colors = c("#00AFBB", "#E7B800")))

    library(gridExtra)
    library(ggplot2)
    p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
    return(p_combine)
  }
  else{ stop(cat("Invalid input data type. Only supported for SummarizedExperiment objects or a list of SummarizedExperiment objects.
        generated from \"Create_MultiAssay_object()\". Your input:", class(result_list)))}

}
########################
#' Obtain two-sample t-test pvalues and emprirical AUC for signature scores.
#' @name get_pvalue_auc
#' @param SE_scored A summarized experiment from TB signature profiling
#' @param annotationColName A character, column name for TB status (default is "Disease")
#' @param signatureColNames A character or vector that contains gene signature name
#' @return A data frame contains p-value from two-sample t-test and AUC value for each signature
#'
#' @examples
#' sobject_TBSig <- get_sobject_TBSig(mobject,"PTB","Latent")
#' SE_scored <- runTBsigProfiler(input = sobject_TBSig, useAssay = "counts", signatures = TBsignatures, algorithm = "ssGSEA", combineSigAndAlgorithm = TRUE, parallel.sz = 1)
#' annotationColName <- "Disease"
#' signatureColNames <- c("Anderson_42","Anderson_OD_51" )
#' results <- get_pvalue_auc(SE_scored,"Disease",signatureColNames)
#' @export
get_pvalue_auc <- function(SE_scored, annotationColName = "Disease", signatureColNames){

  # check annotationColName
  if (length(unique(colData(SE_scored)[annotationColName][,1]))!=2){
    stop(cat("Invalid \"annotationColName\". Expect two-level factor"))
  }

  # check signatureColNames
  if (!all(signatureColNames %in% colnames(colData(SE_scored)))){
    stop(cat("Cannot find score information for",
             signatureColNames[!signatureColNames %in% colnames(colData(SE_scored))]))
  }
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

  aucs_result_dat_lines <- data.frame(Signature = aucs_result_dat$Signature,x0=0.5)

  p_ridge <- ggplot(aucs_result_dat,aes(x=AUC,y=Signature)) + geom_density_ridges(jittered_points=TRUE,alpha=0.7,quantile_lines = TRUE, quantiles = 2) +
    geom_segment(data = aucs_result_dat_lines, aes(x = x0, xend = x0, y = as.numeric(Signature),
                                                   yend = as.numeric(Signature) + .9), color = "red")
  return(p_ridge)
}
