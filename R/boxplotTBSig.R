#' Subset signatures scores and disease status from column data of SummarizedExperiment objects.
#' @name .signature_filter
#' @param object_list A `list` of [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class] objects.
#' Usually output from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param gset A character indicates the name of the signatures.
#' @param annotationColName A character indicates the name of interest in the object's column data.
#' @return A `list` of `data.frame` contains the `annotationColName`, prediction score and study name of single object.
.signature_filter <- function(object_list, gset, annotationColName) {
  obj_name <- base::names(object_list)
  sig_list1 <- lapply(base::seq_len(base::length(object_list)), function(i, gset) {
    x <- object_list[[i]]
    GSE <- base::rep(base::names(object_list[i]),
                     base::nrow(SummarizedExperiment::colData(x)))
    index <- stats::na.omit(base::match(gset,
                                        base::colnames(SummarizedExperiment::colData(x))))
    if (base::length(index) == 0) {
      message(sprintf("Gene signature: %s not found in study: %s, NA is returned.",
                      gset, obj_name[i]))
      result <- base::data.frame(SummarizedExperiment::colData(x)[, annotationColName],
                                  NA, GSE = GSE)
    } else {
      result <- base::data.frame(SummarizedExperiment::colData(x)[, annotationColName],
                                 SummarizedExperiment::colData(x)[, index], GSE = GSE)
    }
    base::colnames(result) <- c(annotationColName, gset, "Study")
    result
  }, gset)
  base::names(sig_list1) <- obj_name
  return(sig_list1)
}

#' Boxplot functions for list of signature scores across studies
#' Also applies to consistent signature column names.
#' @name BoxplotTBSig
#' @param object_list A `list` of [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class] objects.
#' Usually output from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}
#' @param gset A character of vector of characters represent name(s) of the signatures.
#' See \code{\link[TBSignatureProfiler]{TBsignatures}} for details.
#' @param annotationColName A character indicates the name of interest in the object's column data.
#' @return A `gtable` contains single signature's performance across multiple studies.
#' @export
boxplotTBSig <- function(object_list, gset, annotationColName = "TBStatus") {
  # check names of the input object list
  obj_name <- base::names(object_list)
  if (base::is.null(obj_name)) {
    base::stop("Names of the input list should not be NULL. Add unique names for each object within the list.")
  } else if (!base::is.na(base::match("", obj_name))) {
    base::stop(base::sprintf("Names of the input contains \"\". Replace \"\" with a non-empty character."))
  }
  sig_list <- .signature_filter(object_list, gset, annotationColName)
  # Combine list of data frame
  rbindx <- function(dfs) {
    ns <- base::unique(base::unlist(base::sapply(dfs, base::colnames)))
    base::do.call(base::rbind, base::lapply(dfs, function(x) {
      for (i in ns[!ns %in% base::colnames(x)]) {
        x[[n]] <- NA
      }
      x
    }
    ))
  }
  sig_data <- rbindx(sig_list)
  p_boxplot <- base::lapply(base::unique(sig_data$Study), function(x, gset) {
    sig_data_gse <- sig_data %>%
      dplyr::filter(.data$Study == x)
    sig_data_gse$annotationNameLevels <- base::factor(sig_data_gse[, annotationColName])
    is_gset_na <- sig_data_gse %>%
      dplyr::select(gset) %>%
      base::is.na() %>%
      base::all()

    if (is_gset_na) {
      message(sprintf("Gene signature: %s not found from the input list. NULL is returned",
                      gset))
      return(NULL)
    }
    sig_data1 <-  SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)

    # Create a custom color scale to deal with different factors
    n <- base::length(base::levels(sig_data_gse$annotationNameLevels))
    myColors <- base::suppressWarnings(RColorBrewer::brewer.pal(n, "Set1"))
    base::names(myColors) <- base::levels(sig_data_gse$annotationNameLevels)

    p <-  TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                name = x,
                                                signatureColNames = gset,
                                                annotationColName = "annotationNameLevels",
                                                rotateLabels = FALSE,
                                                fill_colors = myColors)

    p1 <- p + ggplot2::theme(plot.title = ggplot2::element_text(size = 12,
                                                                face = "bold"),
                             legend.position = "none",
                             axis.title.x = ggplot2::element_blank(),
                             axis.title.y = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_text(colour = "Black",
                                                                 size = 12, hjust = 0.5,
                                                                 face = "bold"),
                             axis.text.y = ggplot2::element_text(size = 12, angle = 0,
                                                                 hjust = 0.5))
    p1
  }, gset)
  # Remove empty element from list
  p_boxplot <- p_boxplot[!base::sapply(p_boxplot, base::is.null)]
  n_sqrt <- base::sqrt(base::length(p_boxplot))
  p_combine <- base::do.call("grid.arrange",
                             c(p_boxplot, ncol = base::floor(n_sqrt)))
  return(p_combine)
}
