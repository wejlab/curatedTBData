#' Generate Boxplot for single gene signature scores across multiple studies
#'
#' @name BoxplotTBSig
#' @param object_list A \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects. Usually output from
#'   \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param annotationColName A character indicates the name of interest in the
#'   object's column data.
#' @param signatureColNames A character of vector of characters
#'   represent name(s) of the signatures.
#'   See \code{\link[TBSignatureProfiler]{TBsignatures}} for details.
#' @return A \code{gtable} object that contains multiple Boxplots in the form of
#'   \code{ggplot} objects. Results show single signature's performance
#'   across multiple studies.
#' @export
#' @examples
#' returned_resources <- curatedTBData(c("GSE107104", "GSE19435"),
#'                                     dryrun = FALSE, curated.only = TRUE)
#' mysignatures <- list(Sweeney_OD_3 = c("DUSP3", "GBP5", "KLF2"))
#' re1 <- lapply(returned_resources, function(x)
#'                    subset_curatedTBData(x, "TBStatus", c("Control","PTB")))
#' re2 <- lapply(re1, function(x)
#'              TBSignatureProfiler::runTBsigProfiler(input = x,
#'                                                    useAssay = 1,
#'                                                    signatures = mysignatures,
#'                                                    algorithm = "ssGSEA",
#'                                                    update_genes = FALSE))
#' boxplotTBSig(object_list = re2,
#'              annotationColName = "TBStatus",
#'              signatureColNames = "Sweeney_OD_3")
boxplotTBSig <- function(object_list, annotationColName, signatureColNames) {
    .check_input(object_list)
    sig_list <- .signature_filter(object_list, annotationColName, signatureColNames)
    ## Combine list of data frame
    rbindx <- function(dfs) {
        ns <- lapply(dfs, colnames) |>
            unlist() |>
            unique()
        do.call(rbind, lapply(dfs, function(x) {
            for (n in ns[!ns %in% colnames(x)]) {
                x[[n]] <- NA
            }
            x
        }))
    }
    sig_data <- rbindx(sig_list)
    ## Check whether annotationColName or signatureColNames are all NA's
    if (all(is.na(sig_data[, signatureColNames]))) {
        stop(sprintf("Gene signature: %s is not found from the entire input.",
                     signatureColNames))
    } else if (all(is.na(sig_data[, annotationColName]))) {
        stop(sprintf("Annotation name: %s is not found from the entire input.",
                     annotationColName))
    }
    study_name <- unique(sig_data$Study)
    p_boxplot <- lapply(study_name, function(x, signatureColNames) {
        sig_data_gse <- sig_data |>
            dplyr::filter(.data$Study == x)
        sig_data_gse$anno_names <- sig_data_gse[, annotationColName] |>
            factor()
        is_signatureColNames_na <- sig_data_gse |>
            dplyr::select(signatureColNames) |>
            is.na() |>
            all()
        if (is_signatureColNames_na) {
            sprintf("Gene signature: %s not found from the input list.",
                    signatureColNames) |>
                paste("NULL is returnted") |>
                message()
            return(NULL)
        }
        sig_data1 <-
            SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)
        ## Create a custom color scale to deal with different factors
        anno_levels <- levels(sig_data1$anno_names)
        n <- length(anno_levels)
        if (n > 9L) {
            sprintf("The number of levels under %s",
                          annotationColName) |>
                paste("is greater than 9.")  |>
                paste("Only first 9 levels is included.") |>
                message()
            n <- 9
            ## Remove observations for the last level
            sig_data1 <- sig_data1[, sig_data1$anno_names %in% anno_levels[1:9]]
            sig_data1$anno_names <- factor(sig_data1$anno_names)
        } else if (n <= 3L) {
            n <- 3
        }
        myColors <- RColorBrewer::brewer.pal(n, "Set1")
        names(myColors) <- levels(sig_data1$anno_names)
        p <-
            TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                  name = x,
                                                  signatureColNames = signatureColNames,
                                                  annotationColName = "anno_names",
                                                  rotateLabels = FALSE,
                                                  fill_colors = myColors)
        p1 <- p +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 12,
                                                              face = "bold"),
                           legend.position = "none",
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(colour = "Black",
                                                               size = 12,
                                                               hjust = 0.5,
                                                               face = "bold"),
                           axis.text.y = ggplot2::element_text(size = 12,
                                                               angle = 0,
                                                               hjust = 0.5))
        p1
    }, signatureColNames)
    ## Remove empty element from list
    p_boxplot <- p_boxplot[!vapply(p_boxplot, is.null, TRUE)]
    n_sqrt <- sqrt(length(p_boxplot))
    p_combine <- do.call(gridExtra::"grid.arrange",
                         c(p_boxplot, ncol = floor(n_sqrt)))
    return(p_combine)
}

#' Subset signatures scores and disease status from \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} objects
#' @name .signature_filter
#' @param object_list A \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects. Usually output from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param annotationColName A character indicates the name of interest in the
#'   object's column data.
#' @param signatureColNames A character indicates the name of the signatures.
#' @return A \code{list} of \code{data.frame} contains the
#'   \code{annotationColName}, prediction score,
#'   and study name of single object.
.signature_filter <- function(object_list, annotationColName,
                              signatureColNames) {
    obj_name <- names(object_list)
    object_list_seq <- seq_len(length(object_list))
    sig_list1 <- lapply(object_list_seq, function(i, signatureColNames) {
        x <- object_list[[i]]
        col_info <- SummarizedExperiment::colData(x)
        GSE <- rep(names(object_list[i]), nrow(col_info))
        index <- match(signatureColNames, colnames(col_info)) |>
            stats::na.omit()
        anno <- col_info[, annotationColName]
        if (length(index) == 0L) {
            sprintf("Gene signature: %s not found in study: %s,",
                    signatureColNames, obj_name[i]) |>
                paste("NA is returned.") |>
                message()
            result <- data.frame(anno, NA, GSE)
        } else {
            result <- data.frame(anno, col_info[, index], GSE)
        }
        colnames(result) <- c(annotationColName, signatureColNames, "Study")
        result
    }, signatureColNames)
    names(sig_list1) <- obj_name
    return(sig_list1)
}
