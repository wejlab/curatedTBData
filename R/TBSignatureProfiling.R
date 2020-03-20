#######################################
#' Subset signatures scores and disease status from Coldata of SummarizedExperiment Objects.
#' @param sig_list SummarizedExperiment object(s) produced from `TBSignatureProfiler::runTBsigProfiler`.
#' @param gset A vector contians name(s) of the signatures. See `TBSignatureProfiler::TBSignatureProfiler` for example.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param ... Extra named arguments passed to function.
#' @rdname SignatureFilter
#' @export
setGeneric(name="SignatureFilter", function(sig_list, gset,...){
  standardGeneric("SignatureFilter")
})

#' @rdname SignatureFilter
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "list"),
          function(sig_list, gset, annotationColName="TBStatus"){
            sig_list <- lapply(1:length(sig_list), function(y,gset){

              x <- sig_list[[y]]
              GSE <- rep(names(sig_list[y]), nrow(colData(x)))
              TBStatus <- colData(x)[,annotationColName]
              index <- na.omit(match(names(gset),names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            return(sig_list)

          }
)

#' @rdname SignatureFilter
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset, annotationColName="TBStatus"){
            sig_list <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(colData(x)))
              TBStatus <- colData(x)[,annotationColName]
              index <- na.omit(match(gset,names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            return(sig_list)

          }
)
# tt = SignatureFilter(ssgsea_PTB_Latent,TBsignatures)

#' @rdname SignatureFilter
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset, annotationColName="TBStatus"){
            sig_list <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(colData(x)))
              TBStatus <- colData(x)[,annotationColName]
              index <- na.omit(match(gset,names(colData(x))))
              result <- DataFrame(TBStatus,colData(x)[,index],GSE)
              colnames(result)[2] <- gset
              result

            }, gset)
            return(sig_list)

          }
)

# k = SignatureFilter(ssgsea_PTB_Latent,"Anderson_42")

#' @rdname SignatureFilter
setMethod("SignatureFilter",
          signature (sig_list = "SummarizedExperiment", gset = "list"),
          function(sig_list, gset, GSE, annotationColName="TBStatus"){
            index <- na.omit(match(names(gset),names(colData(sig_list))))
            TBStatus <- colData(sig_list)[,annotationColName]
            result <- DataFrame(cbind(TBStatus,colData(sig_list)[,index],GSE))

          }
)

# kk = SignatureFilter(ssgsea_PTB_Latent[[1]],TBsignatures, GSE="test")

#' @rdname SignatureFilter
setMethod("SignatureFilter",
          signature(sig_list = "SummarizedExperiment", gset = "character"),
          function(sig_list, gset, GSE, annotationColName="TBStatus"){
            index <- na.omit(match(gset,names(colData(sig_list))))
            TBStatus <- colData(sig_list)[,annotationColName]
            result <- DataFrame(TBStatus,colData(sig_list)[,index], GSE)
            colnames(result)[2] <- gset
            return(result)

          }
)

# kkk = SignatureFilter(ssgsea_PTB_Latent[[1]],"Anderson_42", GSE="test")

#########################################

#' Boxplot functions for list of signature scores across studies (particularly, for inconsistent signatures within study)
#' Also applies to consistent siganture column names.
#' @param sig_list List of dataframes contain signature scores across studies. Prduced from `SignatureFilter`
#' @param gset A vector contians name(s) of the signatures. See `TBSignatureProfiler::TBSignatureProfiler` for examples.
#' @param annotationColName A character indicates feature of interest in column names.
#' @param ... Extra named arguments passed to function
#' @rdname BoxplotTBSig
#' @export

setGeneric("BoxplotTBSig", function(sig_list, gset, ...) standardGeneric("BoxplotTBSig"))

# sig_list <- MDP_result_NULL;x <- "GSE107993";annotationName = "TBStatus"
# sig_list <- MDP_result;gset = "Anderson_42";x = "GSE56153"

#' @rdname BoxplotTBSig
setMethod("BoxplotTBSig", signature (sig_list = "list", gset = "character"),
          function(sig_list, gset = gset, annotationName = "TBStatus"){

            sig_data <- plyr::rbind.fill(lapply(sig_list,function(x){as.data.frame(x)}))

            p_boxplot <- lapply(unique(sig_data$GSE), function(x, gset){
              sig_data_gse <- sig_data %>% filter(GSE == x)
              sig_data_gse$annotationNameLevels <- factor(sig_data_gse[,annotationName], levels = c("Control", "Latent", "PTB", "OD"))

              if(sig_data_gse %>% dplyr::select(gset) %>% is.na() %>% all()){return(NULL)}

              sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)

              # Create a custom color scale to deal with different factors
              myColors <- RColorBrewer::brewer.pal(4,"Set1")
              names(myColors) <- levels(sig_data_gse$annotationNameLevels)

              p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                            name = x,
                                                            signatureColNames = gset,
                                                            annotationColName = "annotationNameLevels",
                                                            rotateLabels = FALSE,
                                                            fill_colors = myColors)
                                                            # c("#999999", "#E69F00", "#56B4E9", "#FC4E07"))

              return(p)


            }, gset)

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
#' @param SE_scored A SummarizedExperiment Object from TB signature profiling
#' @param annotationColName A character indicates feature of interest in the object's column data
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
#' @param result_list A list of SummarizedExperiment Object from `TBSignatureProfiler::TBSignatureProfiler`.
#' @return A data frame with features including signatures, p-value, and AUC.
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
#' @param aucs_result_dat A dataframe contains signatures, p-value, and AUC, can be obtained from `combine_auc`.
#' @return Ridge plot with median line
#'
#' @examples
#' aucs_result <- data.frame(Signature=c("Anderson_42", "Anderson_OD_51", "Berry_393"), AUC=rnorm(3))
#' p_ridge <- get_auc_distribution(aucs_result)
#' @export
get_auc_distribution <- function(aucs_result){
  library(ggplot2)
  library(gridExtra)
  library(ggridges)

  # add 50% AUC line
  aucs_result_dat_lines <- data.frame(Signature = aucs_result$Signature,x0=0.5)

  p_ridge <- ggplot(aucs_result,aes(x=AUC,y=Signature)) + geom_density_ridges(jittered_points=TRUE,alpha=0.7,quantile_lines = TRUE, quantiles = 2) + geom_segment(data = aucs_result_dat_lines, aes(x = x0, xend = x0, y = as.numeric(Signature),
                                                   yend = as.numeric(Signature) + .9), color = "red")
  return(p_ridge)
}
