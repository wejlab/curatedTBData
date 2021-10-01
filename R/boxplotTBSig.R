#' Generate Boxplot for single gene signature scores across multiple studies
#'
#' @name BoxplotTBSig
#' @param object_list A \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects. Usually output from
#'   \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param gset A character of vector of characters
#'   represent name(s) of the signatures.
#'   See \code{\link[TBSignatureProfiler]{TBsignatures}} for details.
#' @param annotationColName A character indicates the name of interest in the
#'   object's column data.
#' @return A \code{gtable} object that contains multiple Boxplots in the form of
#'   \code{ggplot} objects. Results show single signature's performance
#'   across multiple studies.
#' @export
#' @examples
#' returned_resources <- curatedTBData(c("GSE107104", "GSE19435", "GSE19443"),
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
#'              gset = "Sweeney_OD_3",
#'              annotationColName = "TBStatus")
boxplotTBSig <- function(object_list, gset, annotationColName) {
    ## check names of the input object list
    obj_name <- names(object_list)
    if (is.null(obj_name)) {
        paste("Names of the input list should not be NULL.",
                    "Add unique names for each object within the list.") |>
            stop(call. = FALSE)
    } else if (!is.na(match("", obj_name))) {
        paste("Names of the input contains \"\".",
                    "Replace \"\" with a non-empty character.") |>
            stop(call. = FALSE)
    }
    sig_list <- .signature_filter(object_list, gset, annotationColName)
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
    study_name <- unique(sig_data$Study)
    p_boxplot <- lapply(study_name, function(x, gset) {
        sig_data_gse <- sig_data |>
            filter(.data$Study == x)
        sig_data_gse$anno_names <- sig_data_gse[, annotationColName] |>
            factor()
        is_gset_na <- sig_data_gse |>
            select(gset) |>
            is.na() |>
            all()
        if (is_gset_na) {
            sprintf("Gene signature: %s not found from the input list.",
                    gset) |>
                paste("NULL is returnted") |>
                message()
            return(NULL)
        }
        sig_data1 <-
            SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)
        ## Create a custom color scale to deal with different factors
        n <- sig_data_gse$anno_names |>
            levels() |>
            length()
        if (n > 9L) {
            sprintf("The number of levels under %s",
                          annotationColName) |>
                paste("is greater than 9.")  |>
                paste("Only first 9 levels is included.") |>
                message()
            n <- 9
        } else if (n <= 3L) {
            n <- 3
        }
        myColors <- RColorBrewer::brewer.pal(n, "Set1")
        names(myColors) <- sig_data_gse$anno_names |>
            levels()
        p <-
            TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                  name = x,
                                                  signatureColNames = gset,
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
    }, gset)
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
#' @param gset A character indicates the name of the signatures.
#' @param annotationColName A character indicates the name of interest in the
#'   object's column data.
#' @return A \code{list} of \code{data.frame} contains the
#'   \code{annotationColName}, prediction score,
#'   and study name of single object.
.signature_filter <- function(object_list, gset, annotationColName) {
    obj_name <- names(object_list)
    object_list_seq <- seq_len(length(object_list))
    sig_list1 <- lapply(object_list_seq, function(i, gset) {
        x <- object_list[[i]]
        col_info <- SummarizedExperiment::colData(x)
        GSE <- rep(names(object_list[i]), nrow(col_info))
        index <- match(gset, colnames(col_info)) |>
            stats::na.omit()
        anno <- col_info[, annotationColName]
        if (length(index) == 0L) {
            sprintf("Gene signature: %s not found in study: %s,",
                    gset, obj_name[i]) |>
                paste("NA is returned.") |>
                message()
            result <- data.frame(anno, NA, GSE)
        } else {
            result <- data.frame(anno, col_info[, index], GSE)
        }
        colnames(result) <- c(annotationColName, gset, "Study")
        result
    }, gset)
    names(sig_list1) <- obj_name
    return(sig_list1)
}
