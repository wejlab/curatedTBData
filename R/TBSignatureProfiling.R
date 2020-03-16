#######################################
#' Get signatures with only siganture scores and disease status

setGeneric(name="SignatureFilter", function(sig_list, gset,...){
  standardGeneric("SignatureFilter")
})

setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "list"),
          function(sig_list, gset){
            sig_list <- lapply(1:length(sig_list), function(y,gset){

              x <- sig_list[[y]]
              GSE <- rep(names(sig_list[y]), nrow(colData(x)))
              TBStatus <- colData(x)[,"TBStatus"]
              index <- na.omit(match(names(gset),names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            return(sig_list)

          }
)

setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset){
            sig_list <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(colData(x)))
              TBStatus <- colData(x)[,"TBStatus"]
              index <- na.omit(match(gset,names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            return(sig_list)

          }
)
# tt = SignatureFilter(ssgsea_PTB_Latent,TBsignatures)

setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset){
            sig_list <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(colData(x)))
              TBStatus <- colData(x)[,"TBStatus"]
              index <- na.omit(match(gset,names(colData(x))))
              result <- DataFrame(TBStatus,colData(x)[,index],GSE)
              colnames(result)[2] <- gset
              result

            }, gset)
            return(sig_list)

          }
)

# k = SignatureFilter(ssgsea_PTB_Latent,"Anderson_42")

setMethod("SignatureFilter",
          signature (sig_list = "SummarizedExperiment", gset = "list"),
          function(sig_list, gset, GSE){
            index <- na.omit(match(names(gset),names(colData(sig_list))))
            TBStatus <- colData(sig_list)[,"TBStatus"]
            result <- DataFrame(cbind(TBStatus,colData(sig_list)[,index],GSE))

          }
)

# kk = SignatureFilter(ssgsea_PTB_Latent[[1]],TBsignatures, GSE="test")

setMethod("SignatureFilter",
          signature(sig_list = "SummarizedExperiment", gset = "character"),
          function(sig_list, gset, GSE){
            index <- na.omit(match(gset,names(colData(sig_list))))
            TBStatus <- colData(sig_list)[,"TBStatus"]
            result <- DataFrame(TBStatus,colData(sig_list)[,index], GSE)
            colnames(result)[2] <- gset
            return(result)

          }
)

# kkk = SignatureFilter(ssgsea_PTB_Latent[[1]],"Anderson_42", GSE="test")

#####################################

#' Boxplot functions for list with Signatures matched for each object (particularly, for inconsistent mattching)
#' Also applies to consistent sigantures column names

# sig_list <- MDP_result_NULL;x <- "GSE107993";annotationName = "TBStatus"

setGeneric("BoxplotTBSig", function(sig_list, sig_name, ...) standardGeneric("BoxplotTBSig"))

setMethod("BoxplotTBSig", signature (sig_list = "list", sig_name = "NULL"),
          function(sig_list, sig_name = NULL, annotationName = "TBStatus"){
            sig_data <- plyr::rbind.fill(lapply(sig_list,function(x){as.data.frame(x)}))
            sig_name <- "sMDP"

            p_boxplot <- lapply(unique(sig_data$GSE), function(x, sig_name){
              sig_data_gse <- sig_data %>% filter(GSE == x)

              #sig_data_gse_order <- c("Conrol", "Latent","PTB", "OD","NA")
              # Re-order gene siganture, re-level
              #aucs_result_dat$Signature <- factor(aucs_result_dat$Signature, levels = Signature_order)

              # if(sig_data_gse %>% dplyr::select(sig_name) %>% is.na() %>% all()){return(NULL)}
              sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)

              p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                            name = x,
                                                            signatureColNames = sig_name,
                                                            annotationColName = annotationName,
                                                            rotateLabels = FALSE,
                                                            fill_colors = c("#999999", "#E69F00", "#56B4E9", "#FC4E07"))

              return(p)


            }, sig_name)

            library(gridExtra)
            library(ggplot2)
            p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
            return(p_combine)
          })

#sig_list <- ssgsea_PTB_Latent1;sig_name = "Zak_RISK_16";annotationName = "TBStatus";x = "GSE107993"

setMethod("BoxplotTBSig", signature (sig_list = "list", sig_name = "character"),
          function(sig_list, sig_name = sig_name, annotationName = "TBStatus"){

            sig_data <- plyr::rbind.fill(lapply(sig_list,function(x){as.data.frame(x)}))

            p_boxplot <- lapply(unique(sig_data$GSE), function(x, sig_name){
              sig_data_gse <- sig_data %>% filter(GSE == x)


              if(sig_data_gse %>% dplyr::select(sig_name) %>% is.na() %>% all()){return(NULL)}
              sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)

              p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                            name = x,
                                                            signatureColNames = sig_name,
                                                            annotationColName = annotationName,
                                                            rotateLabels = FALSE,
                                                            fill_colors = c("#999999", "#E69F00", "#56B4E9", "#FC4E07"))

              return(p)


            }, sig_name)

            # Remove empty element from list

            p_boxplot <- plyr::compact(p_boxplot)

            library(gridExtra)
            library(ggplot2)
            p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
            return(p_combine)
          })


################################################

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
combine_auc <- function(result_list, gset, annotationName = "TBStatus"){
  aucs_result <- lapply(result_list, function(x){
    index <- na.omit(match(names(gset),names(colData(x))))
    get_pvalue_auc(x,
                   annotationColName = annotationName,
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
