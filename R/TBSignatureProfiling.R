#' Subset signatures scores and disease status from column data of SummarizedExperiment Objects.
#' @name SignatureFilter
#' @param sig_list SummarizedExperiment object(s) obtained from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' Names of each list should be the GEO accession of the study.
#' @param gset A character/vector contians name(s) of the signatures.
#' @param annotationColName A character indicates feature of interest in the object's column data.
#' @param GSE A character/vector for the GEO accession of the SummarizedExperiment Object.
#' Need GSE value when the input is a single object instead of the list.
#' @param ... Extra named arguments passed to function.
#' @seealso TBSignatureProfiler::TBSignatureProfiler.
#' @rdname SignatureFilter-methods
#' @exportMethod SignatureFilter
setGeneric(name="SignatureFilter", function(sig_list, gset,...){
  standardGeneric("SignatureFilter")
})

#' @rdname SignatureFilter-methods
setMethod("SignatureFilter",
          signature (sig_list = "list", gset = "character"),
          function(sig_list, gset, annotationColName = "TBStatus"){
            if(!any(is(sig_list[[1]])=="SummarizedExperiment")){
              stop(paste("Should be a list of SummarizedExperiment",
                         "Your list element is ",class(sig_list[[1]])))
            }
            sig_list1 <- lapply(1:length(sig_list), function(i,gset){

              x <- sig_list[[i]]
              GSE <- rep(names(sig_list[i]), nrow(SummarizedExperiment::colData(x)))
              index <- stats::na.omit(match(gset,names(SummarizedExperiment::colData(x))))
              if(length(index) == 0) {
                result <- data.frame(SummarizedExperiment::colData(x)[,annotationColName],
                                     NA,GSE=GSE)
                colnames(result) <- c(annotationColName, gset,"GSE")
                return(result)
              }
              result <- data.frame(SummarizedExperiment::colData(x)[,annotationColName],
                                   SummarizedExperiment::colData(x)[,index],GSE=GSE)
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
            index <- stats::na.omit(match(gset,
                                          names(SummarizedExperiment::colData(sig_list))))
            if(length(index) == 0){
              result <- data.frame(SummarizedExperiment::colData(sig_list)
                                   [,annotationColName],
                                   NA, GSE)
              colnames(result) <- c(annotationColName, gset,"GSE")
              return(result)
            }

            result <- data.frame(SummarizedExperiment::colData(sig_list)[,annotationColName],
                                 SummarizedExperiment::colData(sig_list)[,index],
                                 GSE)
            colnames(result) <- c(annotationColName, gset,"GSE")
            return(result)

          }
)


#' Boxplot functions for list of signature scores across studies
#' Also applies to consistent siganture column names.
#' @name BoxplotTBSig
#' @param sig_list List of dataframes contain signature scores across studies.
#' Obtained from \code{\link{SignatureFilter}}.
#' @param gset A vector contians name(s) of the signatures.
#' See \code{\link[TBSignatureProfiler]{TBsignatures}} for details.
#' @param annotationColName A character indicates feature of interest in column names.
#' @param ... Extra named arguments passed to function
#' @return A list of boxplot in the form of ggplot object for each TB gene sigantures across TB studies.
#' @rdname BoxplotTBSig-methods
#' @exportMethod BoxplotTBSig

setGeneric("BoxplotTBSig", function(sig_list, gset, ...)
  standardGeneric("BoxplotTBSig"))

#' @rdname BoxplotTBSig-methods
setMethod("BoxplotTBSig", signature (sig_list = "list", gset = "character"),
          function(sig_list, gset = gset, annotationColName = "TBStatus"){

            # Clean column data, extracting information about annotation, signature scores and name of the each dataset
            if(any(methods::is(sig_list[[1]]) == "SummarizedExperiment")){
              sig_list <- SignatureFilter(sig_list, gset = gset,
                                          annotationColName = annotationColName)
            }

            sig_data <- plyr::rbind.fill(lapply(sig_list,function(x)
                                                            {as.data.frame(x)}))

            p_boxplot <- lapply(unique(sig_data$GSE), function(x, gset){
              sig_data_gse <- sig_data %>% dplyr::filter(.data$GSE == x)
              sig_data_gse$annotationNameLevels <- factor(sig_data_gse[,annotationColName],
                                   levels = c("Control", "LTBI", "PTB", "OD",
                                              "Positive", "Negative", "Others"))

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

              p1 <- p + ggplot2::theme(plot.title = ggplot2::element_text(size=26, face="bold"),
                              legend.title = ggplot2::element_blank(),
                              legend.position = "none",
                              legend.text = ggplot2::element_text(size=20),
                              axis.text.x = ggplot2::element_text(colour="Black",
                                                                  size=26, hjust = 0.5,
                                                                  face="bold"),
                              axis.text.y = ggplot2::element_text(size=20, angle = 0,
                                                                  hjust = 0.5))

              return(p1)


            }, gset)

            # Remove empty element from list

            p_boxplot <- plyr::compact(p_boxplot)


            p_combine <- do.call("grid.arrange",
                                  c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
            return(p_combine)



          })


#' @rdname BoxplotTBSig-methods
setMethod("BoxplotTBSig", signature (sig_list = "data.frame", gset = "character"),
          function(sig_list, gset = gset, annotationColName = "TBStatus"){

            sig_data <- sig_list
            signatureColNames <- colnames(sig_data)[stats::na.omit(match(gset,
                                                                         colnames(sig_data)))]

            sig_data$annotationNameLevels <- factor(sig_data[,annotationColName],
                                                    levels = c("Control", "LTBI",
                                                               "PTB", "OD",
                                                               "Positive", "Negative",
                                                               "Others"))

            # Create a custom color scale to deal with different factors
            myColors <- RColorBrewer::brewer.pal(length(levels(sig_data$annotationNameLevels)),"Set1")

            names(myColors) <- levels(sig_data$annotationNameLevels)

            sig_data1 <- SummarizedExperiment::SummarizedExperiment(colData = sig_data)

            p <- TBSignatureProfiler::signatureBoxplot(
                                          inputData = sig_data1,
                                          name = NULL,
                                          signatureColNames = signatureColNames,
                                          annotationColName = "annotationNameLevels",
                                          rotateLabels = FALSE,
                                          fill_colors = myColors)
            p1 <- p + ggplot2::theme(plot.title = ggplot2::element_text(size=26, face="bold"),
                            strip.text = ggplot2::element_text(size=26, face="bold"),
                            legend.title = ggplot2::element_blank(),
                            legend.position = "none",
                            legend.text = ggplot2::element_text(size=20),
                            axis.text.x = ggplot2::element_text(colour="Black",
                                                                size=24,
                                                                hjust = 0.5,
                                                                face="bold"),
                            axis.text.y = ggplot2::element_text(size=22, angle = 0,
                                                                hjust = 0.5))

            return(p1)


          })



#' Obtain ridge plots for emprirical AUC distribution for signature scores.
#' @name get_auc_distribution
#' @param aucs_result A dataframe contains signatures, p-value, and AUC.
#' This data frame can be obtained from \code{\link{combine_auc}}.
#' @return Ridge plot with median line
#' @examples
#' aucs_result <- data.frame(Signature = c("Anderson_42", "Anderson_OD_51", "Berry_393"),
#'                           AUC = stats::runif(3,0,0.5))
#' p_ridge <- get_auc_distribution(aucs_result)
#' @export
get_auc_distribution <- function(aucs_result){

  # add 50% AUC line
  aucs_result_dat_lines <- data.frame(Signature = aucs_result$Signature,x0=0.5)

  p_ridge <- ggplot2::ggplot(aucs_result, ggplot2::aes(x = .data$AUC,
                                                       y = .data$Signature)) +
    ggridges::geom_density_ridges(jittered_points=TRUE,alpha=0.7,
                                  quantile_lines = TRUE, quantiles = 2) +
    ggplot2::geom_segment(data = aucs_result_dat_lines,
                 ggplot2::aes(x = .data$x0,
                              xend = .data$x0,
                              y = as.numeric(.data$Signature),
                              yend = as.numeric(.data$Signature) + .9),
                              color = "red")
  return(p_ridge)
}


#' Obtain ridge plots for emprirical AUC distribution for signature scores.
#' @name heatmap_auc
#' @param combine_dat A dataframe contains signatures, datsets name, and AUC.
#' Such dataset can be obtained from \code{\link[curatedTBData]{combine_auc}}.
#' @param GSE_sig A dataframe contains information about each signature and its traning dataset name.
#' Defult is NULL.
#' @param signatureColNames A character vector. Expect in the format "Name_SignatureType_Number". e.g. "Anderson_OD_51"
#' @param facet Boolean. TRUE if the users want to group signatures into groups.
#' Default is TRUE.
#' @param clustering Boolena. TRUE if the users want to perform clustering of the heatmap using hierarchical clustering.
#' Default is TRUE.
#' @return Heatmap with AUC values. x axis is the expression data, y axis represents signatures.
#' @examples
#' combine_dat_exp <- data.frame(Signature=rep(c("Anderson_42", "Anderson_OD_51",
#'                               "Berry_393","Berry_OD_86","Blankley_5"),2),
#'                AUC = stats::runif(10,0.5,1), GSE=rep(c("GSE39939","GSE19442"), each=5))
#' GSE_sig_exp <- data.frame(Signature=c("Anderson","Anderson","Berry","Berry"),
#'                    GSE=c("GSE39939","GSE39940","GSE19442","GSE19443"))
#' TBsignatures_exp <- c("Anderson_42", "Anderson_OD_51", "Berry_393","Berry_OD_86",
#'                       "Blankley_5")
#' heatmap_auc(combine_dat_exp,GSE_sig_exp, TBsignatures_exp, facet = FALSE)
#' heatmap_auc(combine_dat_exp,GSE_sig_exp, TBsignatures_exp, facet = TRUE)
#'@export
heatmap_auc <- function(combine_dat, GSE_sig = NULL, signatureColNames, facet = TRUE,
                        clustering = TRUE, order_increase_avg = FALSE,
                        x_axis_name = NULL) {

  dat <- cbind(combine_dat[,c("Signature","GSE","AUC")])
  data_wide <- tidyr::spread(dat, .data$Signature, .data$AUC)
  row.names(data_wide) <- data_wide$GSE
  dat_input <- as.matrix(data_wide[,-1])
  dat_input[is.na(dat_input)] <- NA
  if (order_increase_avg){
    datasets_order_name <- apply(dat_input, 1,function(x) mean(x,na.rm=T)) %>%
      sort(decreasing = TRUE) %>% names()
    dat_input <- dat_input[match(datasets_order_name,row.names(dat_input)),]
  }

  if(!is.null(x_axis_name)){
    dat_input <- dat_input[x_axis_name,]
  }

  if(unique(length(data_wide$GSE)) > 1 && clustering){
    # Clustering AUC values if unique(GSE)>1
    dd <- stats::dist(dat_input)
    hc <- stats::hclust(dd)
    dat_input1 <- dat_input[hc$order,]
  } else {
    dat_input1 <- dat_input
  }


  # Get mean AUC for each dataset
  dat_input1<- cbind(dat_input1,Avg=rowMeans(dat_input1, na.rm = TRUE))

  # Trasform into long format
  datta <- reshape2::melt(dat_input1)

  # Get traning data position index
  datta$trian <- FALSE
  index <- NULL
  if (is.null(GSE_sig)){
    index <- NULL
  }
  else{
    for (i in seq_len(nrow(GSE_sig))){
      kk <- datta[grep(GSE_sig$Signature[i],datta$Var2),]
      kk$indx <- row.names(kk)
      indx <- kk[which(as.character(kk$Var1) %in% GSE_sig$GSE[i]),"indx"]
      index <- c(index,indx)
    }
  }

  # Label signature type
  sig_type_temp <- sig_type_temp2 <- unlist(lapply(strsplit(signatureColNames,"_"),
                                                   function(x) x[2]))
  sig_type <- unique(suppressWarnings(sig_type_temp[which(is.na(as.numeric(sig_type_temp2)))]))

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

  if(facet){

    # Functions to create correct traning index in facet grid
    facet_rect_position <- function(datta, frames){

      # Split dataframe into list based on different signature type
      frames_list <- frames %>% dplyr::group_split(.data$sig_typek)
      names(frames_list) <- unlist(lapply(frames_list, function(x) x$sig_typek[1]))

      datta_list <- datta %>% dplyr::group_split(.data$sig_typek)
      names(datta_list) <- unlist(lapply(datta_list, function(x) x$sig_typek[1]))

      # Get the correct index in for traning dataset
      # change sig_type levels from sub list based on characters in full list
      frame_facet1 <- lapply(names(frames_list), function(i){
        num_Var1 <- as.integer(frames_list[[i]]$Var1)
        num_Var2 <- as.integer(as.integer(factor(frames_list[[i]]$Var2,
                                                 levels = unique(datta_list[[i]]$Var2))))
        frames_list[[i]] %>% dplyr::mutate(Var1 = num_Var1, Var2 = num_Var2)
      })

      frame_facet <- do.call(rbind,frame_facet1)
      return(frame_facet)
    }
    frame_facet <- data.frame(facet_rect_position(datta,frames))
    if (nrow(frame_facet) == 0){

      p <- ggplot2::ggplot(data = datta, ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                                      fill = .data$value)) +
        ggplot2::geom_tile() + ggplot2::scale_fill_distiller(palette = "RdPu",
                                                             trans = "reverse") +
        ggplot2::facet_grid(.data$sig_typek ~ ., switch = "y", scales="free", space="free") +
        ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2)), cex=3.5) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                           size = 12, hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 12))
      return(p)
    }

    else{
      p <- ggplot2::ggplot(data = datta, ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                                      fill = .data$value)) +
        ggplot2::geom_tile() + ggplot2::scale_fill_distiller(palette = "RdPu",
                                                             trans = "reverse") +
        ggplot2::facet_grid(.data$sig_typek ~ ., switch = "y", scales="free", space="free") +
        ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2)), cex=3.5) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                           size = 12, hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 12)) +
        ggplot2::geom_rect(data = frame_facet,
                           ggplot2::aes(xmin = .data$Var1-0.5, xmax = .data$Var1+0.5,
                                        ymin = .data$Var2-0.5, ymax = .data$Var2+0.5),
                           size=1, fill=NA, colour="black", inherit.aes = FALSE)

      return(p)
    }

  }
  else{

    p <- ggplot2::ggplot(data = datta, ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                                    fill = .data$value)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2)), cex=3.5) +
      ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
      ggplot2::geom_rect(data=frames2, size=1, fill=NA, colour="black",
                         ggplot2::aes(xmin = .data$Var1 - 0.5,
                                      xmax = .data$Var1 + 0.5,
                                      ymin = .data$Var2 - 0.5,
                                      ymax = .data$Var2 + 0.5)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                         size = 12, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 12))
    return(p)

  }
}

