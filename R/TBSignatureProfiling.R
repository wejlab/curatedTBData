#######################################
#' Subset signatures scores and disease status from Coldata of SummarizedExperiment Objects.
#' @name SignatureFilter
#' @param sig_list SummarizedExperiment object(s) produced from `TBSignatureProfiler::runTBsigProfiler`.
#' @param gset A character/vector contians name(s) of the signatures. See `TBSignatureProfiler::TBSignatureProfiler` for example.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param ... Extra named arguments passed to function.
#' @rdname SignatureFilter-methods
#' @exportMethod SignatureFilter
setGeneric(name="SignatureFilter", function(sig_list, gset,...){
  standardGeneric("SignatureFilter")
})

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset, annotationColName="TBStatus"){
            if(!any(is(sig_list[[1]])=="SummarizedExperiment")){
              stop(cat("Should be a list of SummarizedExperiment","You list element is ",class(sig_list[[1]])))
            }
            sig_list1 <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(SummarizedExperiment::colData(x)))
              index <- na.omit(match(gset,names(SummarizedExperiment::colData(x))))
              if(length(index) == 0) {
                result <- data.frame(SummarizedExperiment::colData(x)[,annotationColName],NA,GSE=GSE)
                colnames(result) <- c(annotationColName, gset,"GSE")
                return(result)
              }
              result <- data.frame(SummarizedExperiment::colData(x)[,annotationColName],SummarizedExperiment::colData(x)[,index],GSE=GSE)
              colnames(result) <- c(annotationColName, gset,"GSE")
              result

            }, gset)
            names(sig_list1) <- names(sig_list)
            return(sig_list1)

          }
)

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature(sig_list = "SummarizedExperiment", gset = "character"),
          function(sig_list, gset, GSE, annotationColName = "TBStatus"){
            index <- na.omit(match(gset,names(colData(sig_list))))
            if(length(index)==0){
              result <- data.frame(SummarizedExperiment::colData(sig_list)[,annotationColName],NA, GSE)
              colnames(result) <- c(annotationColName, gset,"GSE")
              return(result)
            }

            result <- data.frame(SummarizedExperiment::colData(sig_list)[,annotationColName],SummarizedExperiment::colData(sig_list)[,index], GSE)
            colnames(result) <- c(annotationColName, gset,"GSE")
            return(result)

          }
)


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


#sig_list <- ssgsea_progress;gset = "Anderson_42";x = "GSE107994";annotationColName = "Progression"
#BoxplotTBSig(sig_list, gset = gset, annotationColName = "TBStatus")
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
              sig_data_gse$annotationNameLevels <- factor(sig_data_gse[,annotationColName],
                                                          levels = c("Control", "Latent", "PTB", "OD", "Positive", "Negative"))

              if(sig_data_gse %>% dplyr::select(gset) %>% is.na() %>% all()){return(NULL)}

              sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)

              # Create a custom color scale to deal with different factors
              myColors <- RColorBrewer::brewer.pal(length(levels(sig_data_gse$annotationNameLevels)),"Set1")
              names(myColors) <- levels(sig_data_gse$annotationNameLevels)

              p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                            name = x,
                                                            signatureColNames = gset,
                                                            annotationColName = "annotationNameLevels",
                                                            rotateLabels = FALSE,
                                                            fill_colors = myColors)
                                                            # c("#999999", "#E69F00", "#56B4E9", "#FC4E07"))
              p1 <- p + ggplot2::theme(plot.title = ggplot2::element_text(size=26, face="bold"),
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

            sig_data$annotationNameLevels <- factor(sig_data[,annotationColName],
                                                    levels = c("Control", "Latent", "PTB", "OD", "Positive", "Negative"))

            # Create a custom color scale to deal with different factors
            myColors <- RColorBrewer::brewer.pal(length(levels(sig_data$annotationNameLevels)),"Set1")

            names(myColors) <- levels(sig_data$annotationNameLevels)

            sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data)


            p <- TBSignatureProfiler::signatureBoxplot(
              inputData = sig_data1,
              name = NULL,
              signatureColNames = signatureColNames,
              annotationColName = "annotationNameLevels",
              rotateLabels = FALSE,
              fill_colors = myColors)
            p1 <- p + ggplot2::theme(plot.title = element_text(size=26, face="bold"),
                            strip.text = element_text(size=26, face="bold"),
                            legend.title = element_blank(),
                            legend.position = "none",
                            legend.text = element_text(size=20),
                            axis.text.x = element_text(colour="Black", size=24, hjust = 0.5, face="bold"),
                            axis.text.y = element_text(size=22, angle = 0, hjust = 0.5))

            return(p1)


          })


################################################

#' Obtain pvalue, emprirical AUC, and 95% CI for each signature using two-sample t-test, ROCit::rocit, and bootstraping
#' @name get_auc_stats
#' @param SE_scored A SummarizedExperiment Object from TB signature profiling.
#' @param annotationColName A character indicates feature of interest in the object's column data
#' @param signatureColNames A character/vector contains name of gene signature.
#' @param num.boot Number of bootstrapping.
#' @return A data frame/datatable contains p-value from two-sample t-test and AUC value for each signature.
#' @export
get_auc_stats <- function(SE_scored, annotationColName = "TBStatus", signatureColNames,
                      num.boot=NULL, percent=0.95){

  # check signatureColNames
  index <- na.omit(match(signatureColNames,colnames(SummarizedExperiment::colData(SE_scored))))
  signatureColNames <-  colnames(SummarizedExperiment::colData(SE_scored))[index]

  annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][,1] %>% as.character() %>% as.factor()

  # Get lower and upper quantile
  lower <- (1-percent)/2
  upper <- 1-lower

  # get AUC value for each signature along with corresponding datasets
  if (is.null(num.boot)){

    sig_result <- lapply(signatureColNames, function(i, SE_scored, annotationData){
      score <- SummarizedExperiment::colData(SE_scored)[i][,1]

      # Deal with scores that have constant value (mostly from Sloot_HIV_2)
      if (length(unique(score))==1){
        dat <- data.frame(Signature=i,P.value=NA,AUC=NA)

        return(dat)
      }

      pvals <- stats::t.test(score ~ annotationData)$p.value
      pred <- ROCit::rocit(score, annotationData)
      aucs <- max(pred$AUC, 1 - pred$AUC)
      data.frame(Signature=i,P.value=round(pvals,4),AUC=round(aucs,4))


    }, SE_scored, annotationData)

    result <- data.frame(do.call(rbind, sig_result))
    row.names(result) <- NULL

    return(result)

  }

  else{
    # Initialize parallel
    param <- BiocParallel::SerialParam(progressbar=TRUE)

    sig_result <- BiocParallel::bplapply(signatureColNames, function(i, SE_scored, annotationData, lower, upper){
      score <- SummarizedExperiment::colData(SE_scored)[i][, 1]

      # Deal with PLAGE that have constant score (mostly from Sloot_HIV_2)
      if (length(unique(score))==1){
        dat <- data.frame(i,NA,NA,NA,NA)
        colnames(dat) <- c("Signature","P.value","AUC",
                           paste0("CI lower.",lower*100,"%"),paste0("CI upper.",upper*100,"%"))
        return(dat)
      }
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

      LowerAUC <- stats::quantile(bootCI, prob=lower, na.rm=TRUE)
      UpperAUC <- stats::quantile(bootCI, prob=upper, na.rm=TRUE)
      dat <- data.frame(i,round(pvals,4), round(aucs,4),
                        round(LowerAUC,4), round(UpperAUC,4))
      colnames(dat) <- c("Signature","P.value","AUC",
                         paste0("CI lower.",lower*100,"%"),paste0("CI upper.",upper*100,"%"))
      dat
    }, SE_scored, annotationData, lower, upper,BPPARAM = param)

    result <- do.call(rbind, sig_result)
    row.names(result) <- NULL

    return(result)

  }

}
################################
#' Combine results from list. Calculate p-value and AUC values
#' @name combine_auc
#' @param SE_scored_list A list of SummarizedExperiment Object from `TBSignatureProfiler::TBSignatureProfiler`.
#' @param annotationColName A character indicates feature of interest in the object's column data
#' @param signatureColNames A character/vector contains name of gene signature.
#' @param num.boot Number of bootstrapping.
#' @return A data frame with features including Signatures, P.value, neg10xLog(P.value) and AUC for each signature across datasets.
#' @export
combine_auc <- function(SE_scored_list, annotationColName = "TBStatus", signatureColNames,
                        num.boot=NULL, percent=0.95){
  param <- BiocParallel::SerialParam(progressbar=TRUE)
  aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x){
    get_auc_stats(x,annotationColName,
                signatureColNames,num.boot, percent)
  },BPPARAM = param)
  aucs_result_dat <- do.call(rbind,aucs_result)

  # re-order data based on their median AUC
  # Remove NA value
  aucs_result_dat1 <- na.omit(aucs_result_dat)
  aucs_result_dat_median <- aucs_result_dat1 %>% dplyr::group_by(Signature) %>% dplyr::summarise_all(median) %>% dplyr::arrange(desc(AUC))

  # New addition: order signatures based on median AUC values

  Signature_order <- as.character(aucs_result_dat_median$Signature)

  # Re-order gene siganture, re-level
  # this step is to let ridge plot ordered based on median value
  aucs_result_dat$Signature <- factor(aucs_result_dat$Signature, levels = Signature_order)

  # label name of each dataset as "GSE"
  aucs_result_dat$GSE <- gsub("\\..*","",row.names(aucs_result_dat))

  return(aucs_result_dat)
}


##############################################################################
#' Combine results from list. Calculate p-value and AUC values
#' @name bootstrap_mean_CI
#' @param data A data frame contains the interested numeric vector .
#' @param percent A number indicates the percentage of confidence interval.
#' @param method A character indicates the method used for computing bootstrap confidence interval. By default, this is set to `empirical`.
#'
#' @param num.boot Number of bootstrap times.
#' @return A data frame with lower and upper bootstrap confidence interval.
#' @export
#'
bootstrap_mean_CI <- function(data,colName, percent=0.95, method=c("percentile","empirical"), num.boot){

  if(missing(method)){method="empirical"}
  method <- match.arg(method)
  # cat("The method used for bootstrap confidence interval is ",method)
  lower <- (1-percent)/2
  upper <- 1-lower

  x <- unlist(data[,colName])
  x <- na.omit(x) # Remove NA's in PLAGE method


  names(x) <- NULL
  n <- length(x)
  if (n==1){
    xbar <- x
    ci <- data.frame(xbar,NA, NA)
    colnames(ci) <- c("Mean AUC",paste0("CI lower.",lower*100,"%"),paste0("CI upper.",upper*100,"%"))
    row.names(ci) <- NULL

    return(ci)
  }
  # sample mean
  xbar  <-  mean(x)

  # random resamples from x
  bootstrapsample <- sapply(1:num.boot, function(i) sample(x,n, replace=TRUE))

  # Compute the means x∗
  bsmeans <-  colMeans(bootstrapsample)

  if (method == "empirical"){
    # Compute deltastar for each bootstrap sample
    deltastar <-  bsmeans - xbar

    # Find the 0.0.25 and 0.975 quantile for deltastar
    d <-  quantile(deltastar, c(lower, upper),na.rm=TRUE)

    # Calculate the confidence interval for the mean.
    ci  <-  xbar - c(d[2], d[1])

    ci <- data.frame(xbar,ci[1], ci[2])
    colnames(ci) <- c("Mean AUC",paste0("CI lower.",lower*100,"%"),paste0("CI upper.",upper*100,"%"))
    row.names(ci) <- NULL

    return(ci)
  }

  if (method == "percentile"){
    ci_percent <- quantile(bsmeans, c(lower, upper), na.rm=TRUE)
    ci_percent <- data.frame(xbar,ci_percent[1], ci_percent[2])
    colnames(ci_percent) <- c("Mean",paste0("CI lower.",lower*100,"%"),paste0("CI upper.",upper*100,"%"))
    row.names(ci_percent) <- NULL

    return(ci_percent)
  }

}
##############################################################################

#' Obtain ridge plots for emprirical AUC distribution for signature scores.
#' @name get_auc_distribution
#' @param aucs_result_dat A dataframe contains signatures, p-value, and AUC, can be obtained from `combine_auc`.
#' @return Ridge plot with median line
#'
#' @examples
#' aucs_result <- data.frame(Signature=c("Anderson_42", "Anderson_OD_51", "Berry_393"), AUC=ruinf(3,0,0.5))
#' p_ridge <- get_auc_distribution(aucs_result)
#' @export
get_auc_distribution <- function(aucs_result){

  library(gridExtra)

  # add 50% AUC line
  aucs_result_dat_lines <- data.frame(Signature = aucs_result$Signature,x0=0.5)

  p_ridge <- ggplot2::ggplot(aucs_result,aes(x=AUC,y=Signature)) + ggridges::geom_density_ridges(jittered_points=TRUE,alpha=0.7,quantile_lines = TRUE, quantiles = 2) + geom_segment(data = aucs_result_dat_lines, aes(x = x0, xend = x0, y = as.numeric(Signature),
                                                   yend = as.numeric(Signature) + .9), color = "red")
  return(p_ridge)
}

##############################################################################

#' Obtain ridge plots for emprirical AUC distribution for signature scores.
#' @name heatmap_auc
#' @param combine_dat A dataframe contains signatures, datsets name, and AUC, can be obtained from `combine_auc`.
#' @param GSE_sig A dataframe contains information about eacch signature and its traning dataset name.
#' @param signatureColNames A character vector. Expect in the format "Name_SignatureType_Number". e.g. "Anderson_OD_51"
#' @param facet Logic. True if want to group signatures into groups. Default is FLASE.
#' @return Heatmap with AUC values. x axis represents expression data, y axis represents signatures.
#' @examples
#' combine_dat_exp <- data.frame(Signature=rep(c("Anderson_42", "Anderson_OD_51", "Berry_393","Berry_OD_86","Blankley_5"),2),
#'                AUC=runif(10,0.5,1), GSE=rep(c("GSE39939","GSE19442"), each=5))
#' GSE_sig_exp <- data.frame(Signature=c("Anderson","Anderson","Berry","Berry"),GSE=c("GSE39939","GSE39940","GSE19442","GSE19443"))
#' TBsignatures_exp <- c("Anderson_42", "Anderson_OD_51", "Berry_393","Berry_OD_86","Blankley_5")
#' heatmap_auc(combine_dat_exp,GSE_sig_exp, TBsignatures_exp, facet=FALSE)
#' heatmap_auc(combine_dat_exp,GSE_sig_exp, TBsignatures_exp, facet=TRUE)
#'@export
heatmap_auc <- function(combine_dat,GSE_sig, signatureColNames, facet=FALSE){
  dat <- cbind(combine_dat[,c("Signature","GSE","AUC")])
  data_wide <- tidyr::spread(dat, Signature, AUC)
  row.names(data_wide) <- data_wide$GSE
  dat_input <- data_wide[,-1] %>% as.matrix
  dat_input[is.na(dat_input)] <- NA

  # Clustering AUC values
  dd <- dist(dat_input)
  hc <- hclust(dd)
  dat_input1 <-dat_input[hc$order,]

  # Get mean AUC for each dataset
  dat_input1<- cbind(dat_input1,Avg=rowMeans(dat_input1, na.rm = TRUE))

  # Trasform into long format
  datta <- reshape2::melt(dat_input1)

  # Get traning data position index
  datta$trian <- FALSE

  index <- NULL
  for (i in 1:nrow(GSE_sig)){
    kk <- datta[grep(GSE_sig$Signature[i],datta$Var2),]
    kk$indx <- row.names(kk)
    indx <- kk[which(as.character(kk$Var1) %in% GSE_sig$GSE[i]),"indx"]
    index <- c(index,indx)
  }

  # Label signature type
  sig_type_temp <- sig_type_temp2 <- sapply(strsplit(signatureColNames,"_"), function(x) x[2])
  sig_type <- suppressWarnings(sig_type_temp[which(is.na(as.numeric(sig_type_temp2)))]) %>% unique()
  datta$sig_typek <- "Disease"

  for (i in sig_type){
    datta$sig_typek[grep(i,datta$Var2)] <- i
  }
  datta$sig_typek[grep("Avg",datta$Var2)] <- "Avg"

  datta$sig_typek <- factor(datta$sig_typek, levels = c("Avg",sig_type,"Disease"))

  datta[as.numeric(index),"trian"] <- TRUE
  frames2 <-  frames <- datta[datta$trian, c("Var1", "Var2","sig_typek")]

  frames2$Var1 <- as.integer(frames$Var1)
  frames2$Var2 <- as.integer(frames$Var2)

  # Functions to create correct traning index in facet grid
  facet_rect_position <- function(datta, frames){

    # Split dataframe into list based on different signature type
    frames_list <- frames %>% dplyr::group_split(sig_typek)
    names(frames_list) <- sapply(frames_list, function(x) x$sig_typek[1])

    datta_list <- datta %>% dplyr::group_split(sig_typek)
    names(datta_list) <- sapply(datta_list, function(x) x$sig_typek[1])

    # Get the correct index in for traning dataset
    # chnage sig_type levels from sub list based on characters in full list
    frame_facet1 <- lapply(names(frames_list), function(i){
      frames_list[[i]]$Var1 <- as.integer(frames_list[[i]]$Var1)
      frames_list[[i]]$Var2 <- as.integer(factor(frames_list[[i]]$Var2, levels = unique(datta_list[[i]]$Var2)))

      frames_list[[i]]
    })
    frame_facet <- do.call(rbind,frame_facet1)
    return(frame_facet)

  }

  frame_facet <- facet_rect_position(datta,frames)

  if (facet==FALSE){
    p <- ggplot(data = datta, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)), cex=3.5) +
      scale_fill_distiller(palette = "RdPu", trans = "reverse") +
      geom_rect(data=frames2, size=1, fill=NA, colour="black",
                aes(xmin=Var1 - 0.5, xmax=Var1 + 0.5, ymin=Var2 - 0.5, ymax=Var2 + 0.5)) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12))
    return(p)

  }
  if(facet==TRUE){
    p <- ggplot(data = datta, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      scale_fill_distiller(palette = "RdPu", trans = "reverse") +
      facet_grid(sig_typek~., switch = "y", scales="free", space="free") +
      geom_text(aes(label = round(value, 2)), cex=3.5) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1),
            axis.text.y = element_text(size = 12)) +
      geom_rect(data = data.frame(frame_facet),
                aes(xmin = Var1-0.5, xmax = Var1+0.5, ymin = Var2-0.5, ymax = Var2+0.5),
                size=1, fill=NA, colour="black", inherit.aes = FALSE)

    return(p)
  }

}
