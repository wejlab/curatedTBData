#' Obtain signature gene rank
#' @param data_sobject A SummarizedExperiment Object contains assays and column data with disease status.
#' @param TBSig A list of signature genes with list names as signature name.
#' @param annotationColNames A character indicates feature of interest in the object's column data.
#' @return A list of dataframes (single dataframe if only one signature) with Disease status and signature genes rank.
#' @examples
#' data_sobject <- TBSignatureProfiler::TB_hiv
#' TBSig <- TBSignatureProfiler::TBsignatures[c("Sweeney_OD_3","Roe_3")]
#' annotationColNames <- "Disease"
#' get_gsea_rank(data_sobject,TBSig,annotationColNames)
#'
#' @export
get_gsea_rank <- function(data_sobject,TBSig,annotationColNames="TBStatus"){

  data_sobject_assay <- SummarizedExperiment::assay(data_sobject)
  # gsea_result <- GSVA::gsva(data_sobject_assay,TBsignatures,method="ssgsea",ssgsea.norm=F,parallel.sz=size)
  get_gene_rank <- function(data,index){

    geneRanking <- order(data, decreasing=TRUE)
    indicatorFunInsideGeneSet <- match(geneRanking, index)
    indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
    indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
    position <- which(indicatorFunInsideGeneSet==1)
    return(data[geneRanking][position])

  }

  result_sig <- BiocParallel::bplapply(names(TBSig), function(signatureColNames){
    index <- stats::na.omit(match(TBSig[[signatureColNames]],row.names(data_sobject_assay)))
    if(length(index)==0){
      gene_position_summary <- data.frame(Disease = SummarizedExperiment::colData(data_sobject)
                                          [,annotationColNames])
      return(gene_position_summary)

    }
    len <- length(index)

    # Rank the assay first
    data_rank <- apply(data_sobject_assay, 2, rank)
    position <- lapply(seq_len(ncol(data_rank)), function(i)
                         get_gene_rank(data_rank[,i],index=index)) %>% unlist()
    if(len==1){
      gene_position_summary <- data.frame(Disease = SummarizedExperiment::colData(data_sobject)
                                          [,annotationColNames], position)
      return(gene_position_summary)
    }
    gene_position_summary <- data.frame(Disease = SummarizedExperiment::colData(data_sobject)
                                        [,annotationColNames], t(position))
    gene_position_summary
  },BPPARAM = BiocParallel::SerialParam(progressbar=TRUE))

  names(result_sig) <- names(TBSig)
  return(result_sig)

}

#' Obtain boxplot of gene rank boxplot for each signature
#' @param sig_list A (nested) list of data frame
#' @param gset A character indicates the name of the signature
#' @return A group of boxplots show signature gene ranks under difference disease status.
#' @examples
#' data_sobject <- TBSignatureProfiler::TB_hiv
#' TBSig <- TBSignatureProfiler::TBsignatures[c("Sweeney_OD_3","Roe_3")]
#' annotationColNames <- "Disease"
#' sig_list <- get_gsea_rank(data_sobject,TBSig,annotationColNames)
#' get_rank_boxplot(sig_list,"Sweeney_OD_3")
#' @export
get_rank_boxplot <- function(sig_list,gset){

  if(is(sig_list[[1]])=="data.frame"){
    result <- sig_list[[gset]]
    data_long <- suppressMessages({reshape2::melt(result)})

    p <- ggplot2::ggplot(data_long, ggplot2::aes(x = .data$variable,
                                                 y = .data$value,
                                                 fill = .data$Disease)) +
      ggplot2::geom_boxplot() + ggplot2::theme_bw() +
      ggplot2::theme(
            plot.title = ggplot2::element_text(size=12, face="bold"),
            legend.title = ggplot2::element_blank(),
            legend.position = "none",
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(colour="Black", size=12, hjust = 0.5,
                                                angle = 30,face = "bold"),
            axis.text.y = ggplot2::element_text(size=12, angle = 0, hjust = 0.5),
            axis.title.y = ggplot2::element_blank())
    return(p)

  }
  data_list <- lapply(sig_list, function(x) {

    result <- x[[gset]]
    if(ncol(result)==1){result=NULL}
    return(result)
  })
  data_list <- plyr::compact(data_list)

  p_boxplot <- lapply(seq_along(data_list), function(x){
    data_long <- suppressMessages({
      reshape2::melt(data_list[[x]])
    })
    p <- ggplot2::ggplot(data_long, ggplot2::aes(x = .data$variable,
                                                 y = .data$value,
                                                 fill = .data$Disease)) +
      ggplot2::geom_boxplot() + ggplot2::ggtitle(names(data_list)[x]) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(size=33, face="bold"),
                     legend.title = ggplot2::element_blank(),
                     legend.position = "none",
                     axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(colour="Black", size=27, hjust = 0.5, angle = 30,face = "bold"),
            axis.text.y = ggplot2::element_text(size=25, angle = 0, hjust = 0.5),
            axis.title.y = ggplot2::element_blank())
    p
  })

    p_combine <- do.call("grid.arrange",
                         c(p_boxplot, ncol=floor(sqrt(length(p_boxplot)))))
    return(p_combine)


}

