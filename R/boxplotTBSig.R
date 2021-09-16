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
    obj_name <- base::names(object_list)
    if (base::is.null(obj_name)) {
        base::paste("Names of the input list should not be NULL.",
                    "Add unique names for each object within the list.") %>%
            base::stop(call. = FALSE)
    } else if (!base::is.na(base::match("", obj_name))) {
        base::paste("Names of the input contains \"\".",
                    "Replace \"\" with a non-empty character.") %>%
            base::stop(call. = FALSE)
    }
    sig_list <- .signature_filter(object_list, gset, annotationColName)
    ## Combine list of data frame
    rbindx <- function(dfs) {
        ns <- base::lapply(dfs, base::colnames) %>%
            base::unlist() %>%
            base::unique()
        base::do.call(base::rbind, base::lapply(dfs, function(x) {
            for (n in ns[!ns %in% base::colnames(x)]) {
                x[[n]] <- NA
            }
            x
        }))
    }
    sig_data <- rbindx(sig_list)
    study_name <- base::unique(sig_data$Study)
    p_boxplot <- base::lapply(study_name, function(x, gset) {
        sig_data_gse <- sig_data %>%
            dplyr::filter(.data$Study == x)
        sig_data_gse$anno_names <- sig_data_gse[, annotationColName] %>%
            base::factor()
        is_gset_na <- sig_data_gse %>%
            dplyr::select(gset) %>%
            base::is.na() %>%
            base::all()
        if (is_gset_na) {
            base::sprintf("Gene signature: %s not found from the input list.",
                          gset) %>%
                base::paste("NULL is returnted") %>%
                base::message()
            return(NULL)
        }
        sig_data1 <-
            SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)
        ## Create a custom color scale to deal with different factors
        n <- sig_data_gse$anno_names %>%
            base::levels() %>%
            base::length()
        if (n > 9L) {
            base::sprintf("The number of levels under %s",
                          annotationColName) %>%
                base::paste("is greater than 9.") %>%
                base::paste("Only first 9 levels is included.") %>%
                base::message()
            n <- 9
        } else if (n <= 3L) {
            n <- 3
        }
        myColors <- RColorBrewer::brewer.pal(n, "Set1")
        base::names(myColors) <- sig_data_gse$anno_names %>%
            base::levels()
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
    p_boxplot <- p_boxplot[!base::vapply(p_boxplot, base::is.null, TRUE)]
    n_sqrt <- base::sqrt(base::length(p_boxplot))
    p_combine <- base::do.call(gridExtra::"grid.arrange",
                               c(p_boxplot, ncol = base::floor(n_sqrt)))
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
    obj_name <- base::names(object_list)
    object_list_seq <- base::seq_len(base::length(object_list))
    sig_list1 <- lapply(object_list_seq, function(i, gset) {
        x <- object_list[[i]]
        col_info <- SummarizedExperiment::colData(x)
        GSE <- base::rep(base::names(object_list[i]),
                         base::nrow(col_info))
        index <- base::match(gset, base::colnames(col_info)) %>%
            stats::na.omit()
        anno <- col_info[, annotationColName]
        if (base::length(index) == 0L) {
            base::sprintf("Gene signature: %s not found in study: %s,",
                          gset, obj_name[i]) %>%
                base::paste("NA is returned.") %>%
                base::message()
            result <- base::data.frame(anno, NA, GSE)
        } else {
            result <- base::data.frame(anno, col_info[, index], GSE)
        }
        base::colnames(result) <- c(annotationColName, gset, "Study")
        result
    }, gset)
    base::names(sig_list1) <- obj_name
    return(sig_list1)
}
