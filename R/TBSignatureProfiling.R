#######################################
#' Subset signatures scores and disease status from Coldata of SummarizedExperiment Objects.
#' @name SignatureFilter
#' @param sig_list SummarizedExperiment object(s) produced from `TBSignatureProfiler::runTBsigProfiler`.
#' @param gset A vector contians name(s) of the signatures. See `TBSignatureProfiler::TBSignatureProfiler` for example.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param ... Extra named arguments passed to function.
#' @rdname SignatureFilter-methods
#' @exportMethod SignatureFilter
setGeneric(name="SignatureFilter", function(sig_list, gset,...){
  standardGeneric("SignatureFilter")
})

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "list"),
          function(sig_list, gset, annotationColName="TBStatus"){
            if(!any(is(sig_list[[1]])=="SummarizedExperiment")){
              stop(cat("Should be a list of SummarizedExperiment","You list element is ",class(sig_list[[1]])))
            }
            sig_list1 <- lapply(1:length(sig_list), function(y,gset){

              x <- sig_list[[y]]
              GSE <- rep(names(sig_list[y]), nrow(colData(x)))
              TBStatus <- colData(x)[,annotationColName]
              index <- na.omit(match(names(gset),names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            names(sig_list1) <- names(sig_list)
            return(sig_list1)

          }
)

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset, annotationColName="TBStatus"){
            if(!any(is(sig_list[[1]])=="SummarizedExperiment")){
              stop(cat("Should be a list of SummarizedExperiment","You list element is ",class(sig_list[[1]])))
            }
            sig_list1 <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(colData(x)))
              TBStatus <- colData(x)[,annotationColName]
              index <- na.omit(match(gset,names(colData(x))))
              cbind(TBStatus,colData(x)[,index],GSE)

            }, gset)
            names(sig_list1) <- names(sig_list)
            return(sig_list1)

          }
)
# tt = SignatureFilter(ssgsea_PTB_Latent,TBsignatures)

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature (sig_list = "SummarizedExperiment", gset = "list"),
          function(sig_list, gset, GSE, annotationColName="TBStatus"){
            index <- na.omit(match(names(gset),names(colData(sig_list))))
            TBStatus <- colData(sig_list)[,annotationColName]
            result <- DataFrame(cbind(TBStatus,colData(sig_list)[,index],GSE))

          }
)

# kk = SignatureFilter(ssgsea_PTB_Latent[[1]],TBsignatures, GSE="test")

#' @rdname SignatureFilter-methods
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
#' @name BoxplotTBSig
#' @param sig_list List of dataframes contain signature scores across studies. Prduced from `SignatureFilter`
#' @param gset A vector contians name(s) of the signatures. See `TBSignatureProfiler::TBSignatureProfiler` for examples.
#' @param annotationColName A character indicates feature of interest in column names.
#' @param ... Extra named arguments passed to function
#' @rdname BoxplotTBSig-methods
#' @exportMethod BoxplotTBSig

setGeneric("BoxplotTBSig", function(sig_list, gset, ...) standardGeneric("BoxplotTBSig"))

# sig_list <- MDP_result_NULL;x <- "GSE107993";annotationColName = "TBStatus"
# sig_list <- MDP_result;gset = "Anderson_42";x = "GSE56153"

#' @rdname BoxplotTBSig-methods
setMethod("BoxplotTBSig", signature (sig_list = "list", gset = "character"),
          function(sig_list, gset = gset, annotationColName = "TBStatus"){

            # Clean column data, extracting information about annotation, signature scores and name of the each dataset
            if(any(is(sig_list[[1]]) == "SummarizedExperiment")){
              sig_list <- SignatureFilter(sig_list, gset = gset, annotationColName = annotationColName)
            }

            sig_data <- plyr::rbind.fill(lapply(sig_list,function(x){as.data.frame(x)}))

            p_boxplot <- lapply(unique(sig_data$GSE), function(x, gset){
              sig_data_gse <- sig_data %>% filter(GSE == x)
              sig_data_gse$annotationNameLevels <- factor(sig_data_gse[,annotationColName], levels = c("Control", "Latent", "PTB", "OD"))

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
              p1 <- p + theme(plot.title = ggplot2::element_text(size=26, face="bold"),
                              legend.title = element_blank(),
                              legend.position = "none",
                              legend.text = ggplot2::element_text(size=20),
                              axis.text.x = ggplot2::element_text(colour="Black", size=26, hjust = 0.5, face="bold"),
                              axis.text.y = ggplot2::element_text(size=20, angle = 0, hjust = 0.5))

              return(p1)


            }, gset)

            # Remove empty element from list

            p_boxplot <- plyr::compact(p_boxplot)

            library(gridExtra)
            library(ggplot2)
            p_combine <- do.call("grid.arrange", c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
            return(p_combine)
          })


#' @rdname BoxplotTBSig-methods
setMethod("BoxplotTBSig", signature (sig_list = "data.frame", gset = "character"),
          function(sig_list, gset = gset, annotationColName = "TBStatus"){

            sig_data <- sig_list
            signatureColNames <- colnames(sig_data)[na.omit(match(gset, colnames(sig_data)))]

            sig_data$annotationNameLevels <- factor(sig_data[,annotationColName], levels = c("Control", "Latent", "PTB", "OD"))
            myColors <- RColorBrewer::brewer.pal(4,"Set1")
            names(myColors) <- levels(sig_data$annotationNameLevels)

            sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data)


            p <- TBSignatureProfiler::signatureBoxplot(
              inputData = sig_data1,
              name = NULL,
              signatureColNames = signatureColNames,
              annotationColName = "annotationNameLevels",
              rotateLabels = FALSE,
              fill_colors = myColors)
            p1 <- p + theme(plot.title = ggplot2::element_text(size=26, face="bold"),
                            strip.text = element_text(size=26, face="bold"),
                            legend.title = element_blank(),
                            legend.position = "none",
                            legend.text = ggplot2::element_text(size=20),
                            axis.text.x = ggplot2::element_text(colour="Black", size=24, hjust = 0.5, face="bold"),
                            axis.text.y = ggplot2::element_text(size=22, angle = 0, hjust = 0.5))

            return(p1)


          })


################################################

#' Obtain pvalue, emprirical AUC, and 95% CI for each signature using two-sample t-test, ROCit::rocit, and bootstraping
#' @name get_stats
#' @param SE_scored A SummarizedExperiment Object from TB signature profiling.
#' @param annotationColName A character indicates feature of interest in the object's column data
#' @param signatureColNames A character/vector contains name of gene signature.
#' @param num.boot Number of bootstrapping.
#' @param output A character specifies types of output, either data.frame or datatable from `DT`
#' @return A data frame/datatable contains p-value from two-sample t-test and AUC value for each signature.
#' @export
get_stats <- function(SE_scored, annotationColName = "TBStatus", signatureColNames,
                      num.boot=NULL, output="data.frame"){
  # check signatureColNames
  index <- na.omit(match(signatureColNames,colnames(SummarizedExperiment::colData(SE_scored))))
  signatureColNames <-  colnames(SummarizedExperiment::colData(SE_scored))[index]

  annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][,1] %>% as.character() %>% as.factor()

  if (is.null(num.boot)){

    sig_result <- lapply(signatureColNames, function(i, SE_scored, annotationData){
      score <- SummarizedExperiment::colData(SE_scored)[i][,1]
      pvals <- stats::t.test(score ~ annotationData)$p.value
      pred <- ROCit::rocit(score, annotationData)
      aucs <- max(pred$AUC, 1 - pred$AUC)
      data.frame(Signature=i,P.value=round(pvals,4),AUC=round(aucs,4))


    }, SE_scored, annotationData)

    result <- data.frame(do.call(rbind, sig_result))
    row.names(result) <- NULL
    if(output == "DataTabble"){
      return(DT::datatable(result))
    }
    else{
      return(result)
    }

  }
  else{

    sig_result <- lapply(signatureColNames, function(i, SE_scored, annotationData){
      score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
      pvals <- stats::t.test(score ~ annotationData)$p.value
      pred <- ROCit::rocit(score, annotationData)
      aucs <- max(pred$AUC, 1 - pred$AUC)
      bootCI <- sapply(1:num.boot, function(j, score, annotationData){

        index <- sample(1:length(score), replace = TRUE)
        tmp_score <- score[index]
        tmp_annotationData <- annotationData[index]
        # Consider when resampling only has 1 cases, remove it
        if(length(unique(tmp_annotationData)) == 2){
          tmp_pred <- ROCit::rocit(tmp_score, tmp_annotationData)
          tmp_auc <- max(tmp_pred$AUC, 1 - tmp_pred$AUC)
          tmp_auc
        }else{NA}

      }, score, annotationData)

      bootCI <- na.omit(bootCI)

      LowerAUC <- stats::quantile(bootCI,prob=0.05)
      UpperAUC <- stats::quantile(bootCI,prob=0.95)
      data.frame(Signature=i,P.value=round(pvals,4),`neg10xLog(P.value)` = round(-10 * log(pvals), 4),
                 AUC=round(aucs,4), LowerAUC=round(LowerAUC,4), UpperAUC=round(UpperAUC,4))

    }, SE_scored, annotationData)
    result <- data.frame(do.call(rbind, sig_result))
    row.names(result) <- NULL
    if(output == "DataTabble"){
      return(DT::datatable(result))
    }
    else{
      return(result)
    }


  }

}

################################
#' Combine results from list. Calculate p-value and AUC values
#' @name combine_auc
#' @param SE_scored_list A list of SummarizedExperiment Object from `TBSignatureProfiler::TBSignatureProfiler`.
#' @param annotationColName A character indicates feature of interest in the object's column data
#' @param signatureColNames A character/vector contains name of gene signature.
#' @param num.boot Number of bootstrapping.
#' @param output A character specifies types of output, either data.frame or datatable from `DT`
#' @return A data frame/datatable with features including signatures, p-value, and AUC for each signature across datasets.
#' @export
combine_auc <- function(SE_scored_list, annotationColName = "TBStatus", signatureColNames,
                        num.boot=NULL, output="data.frame"){
  aucs_result <- lapply(SE_scored_list, function(x){
    get_stats(x,
                   annotationColName = annotationColName,
                   signatureColNames = signatureColNames,
                   num.boot = num.boot)
  }
    )
  aucs_result_dat <- do.call(rbind,aucs_result)

  # re-order data based on their median AUC
  aucs_result_dat_median <- aucs_result_dat %>% dplyr::group_by(Signature) %>% dplyr::summarise_all(median) %>% dplyr::arrange(desc(AUC))

  # New addition: order signatures based on median AUC values

  Signature_order <- as.character(aucs_result_dat_median$Signature)
  # Re-order gene siganture, re-level
  aucs_result_dat$Signature <- factor(aucs_result_dat$Signature, levels = Signature_order)

  return(aucs_result_dat)
}

##############################################################################

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
