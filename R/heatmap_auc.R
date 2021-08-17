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

