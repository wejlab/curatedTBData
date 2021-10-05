#' Combine results from a \code{list} of \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} objects
#'
#' Calculate p-value, AUC values, and bootstrapped confidence interval
#'   for TB gene signatures across multiple studies.
#' @name combine_auc
#' @param SE_scored_list A \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects obtained from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param annotationColName A string specifying the feature of interest in the
#'   object's column data.
#' @param signatureColNames A character/vector string
#'   contains name of gene signature.
#' @param num.boot Integer. Number of bootstrapping.
#'   If \code{NULL} (default), confidence interval will not be computed.
#' @param percent Numeric. A number between 0 and 1,
#'   indicating the percentage of confidence interval. Default is \code{0.95}.
#' @param AUC.abs Boolean. If \code{AUC.abs = TRUE},
#'   compute the AUC values based on \code{\link[ROCit]{rocit}}.
#'   If \code{AUC.abs = FALSE},
#'   compute the AUC values for \code{max(AUC, 1-AUC)}.
#' @param BPPARAM An instance inherited from \code{\link[BiocParallel]{bplapply}}.
#' @return A \code{data.frame} with features including
#'   Signatures, P.value, neg10xLog(P.value) and AUC
#'   for each signature across studies.
#' @examples
#' library(TBSignatureProfiler)
#' data(TB_indian, package = "TBSignatureProfiler")
#' TBsignaturesSub <- TBsignatures[1:5]
#' res <- runTBsigProfiler(input = TB_indian, useAssay = "logcounts",
#'                         signatures = TBsignaturesSub, algorithm = "ssGSEA",
#'                         update_genes = FALSE)
#' re <- combine_auc(list(TB_indian = res), annotationColName = "label",
#'                   signatureColNames = names(TBsignaturesSub), num.boot = 100,
#'                   percent = 0.95)
#' @export
combine_auc <- function(SE_scored_list, annotationColName, signatureColNames,
                        num.boot = NULL, percent = 0.95, AUC.abs = FALSE,
                        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    .check_input(SE_scored_list)
    bpparam <- BPPARAM
    if (is.null(num.boot)) {
        paste("\"num.boot\" is NULL",
              "Bootstrap Confidence Interval is not computed.") |>
            message()
    }
    aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x) {
        .get_auc_stats(x, annotationColName, signatureColNames, num.boot,
                       percent, AUC.abs)
    }, BPPARAM = bpparam)
    aucs_result_dat <- do.call(rbind, aucs_result)
    if (nrow(aucs_result_dat) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found in the list",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. Check \"signatureColNames\".") |>
            stop(call. = FALSE)
    }
    ## Re-order data based on their median AUC (from largest to smallest)
    ## Remove NA value
    aucs_result_dat_median <- aucs_result_dat |>
        dplyr::filter(!is.na(.data$AUC)) |>
        dplyr::group_by(.data$Signature) |>
        dplyr::summarise_all(stats::median) |>
        dplyr::arrange(dplyr::desc(.data$AUC))
    ## Order signatures based on median AUC values
    Signature_order <- as.character(aucs_result_dat_median$Signature)
    Sig_NA <- aucs_result_dat |>
        dplyr::filter(is.na(.data$AUC)) |>
        dplyr::select(.data$Signature) |>
        unlist(use.names = FALSE) |>
        unique()
    ## Re-order gene signature
    ## re-level this step is to let ridge plot ordered based on median value
    sig_levels <- unique(c(Signature_order, Sig_NA))
    aucs_result_dat$Signature <- factor(aucs_result_dat$Signature,
                                        levels = sig_levels)
    ## Label name of each data under column 'Study'
    aucs_result_dat$Study <- gsub("\\..*", "", row.names(aucs_result_dat))
    row.names(aucs_result_dat) <- NULL
    return(aucs_result_dat)
}

#' Obtain pvalue, emprirical AUC, and Bootstrap Confidence Interval for each signature
#' @name .get_auc_stats
#' @param SE_scored A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object from TB signature profiling.
#' @param annotationColName A string specifying the feature of interest in
#'   the object's column data.
#' @param signatureColNames A character/vector string contains name of gene signature.
#' @param num.boot Integer. Number of bootstrapping.
#' @param percent Numeric. A number between 0 and 1,
#'   indicating the percentage of confidence interval.
#' @param AUC.abs Boolean. If \code{AUC.abs = TRUE},
#'   compute the AUC values from function \code{\link[ROCit]{rocit}}.
#'   If \code{AUC.abs = FALSE}, compute the AUC values for \code{max(AUC, 1-AUC)}.
#' @return A \code{data.frame} contains
#'   p-value from two-sample t-test and AUC value for each signature.
.get_auc_stats <- function(SE_scored, annotationColName, signatureColNames,
                           num.boot, percent, AUC.abs) {
    ## Check signatureColNames
    col_info <- SummarizedExperiment::colData(SE_scored)
    index <- match(signatureColNames, colnames(col_info)) |>
        stats::na.omit()
    if (length(index) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. NULL is returned.\n") |>
            message()
    }
    signatureColNames <- colnames(col_info)[index]
    ## Check annotationColName
    index_anno <- match(annotationColName, colnames(col_info))
    if (is.na(index_anno)) {
        sprintf("\"annotationColName\": %s is not found from the study.\n",
                annotationColName) |>
            stop(call. = FALSE)
    }
    annotationData <- col_info[annotationColName][, 1] |>
        as.character() |>
        as.factor()
    ## Check levels of annotationData
    anno_level_len <- length(unique(annotationData))
    if (anno_level_len != 2L) {
        paste("Annotation data should have exactly two levels.",
              "The number of input levels is:",
              anno_level_len, ".\n") |>
            stop(call. = FALSE)
    }
    ## Get AUC value for each signature along with corresponding datasets
    if (is.null(num.boot)) {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData) {
            score <- col_info[i][, 1] |>
                as.vector()
            ## Deal with scores that have constant value (e.g. Sloot_HIV_2)
            if (length(unique(score)) == 1L) {
                sprintf(paste("Constant score found for siganture: %s,",
                              "results will be NA.\n"), i) |>
                    message()
                dat <- data.frame(Signature = i, P.value = NA,
                                  neg10xP.value = NA, AUC = NA)
                return(dat)
            }
            pvals <- stats::t.test(score ~ annotationData)$p.value
            neg10log <- -1 * log(pvals + 1e-4)
            pred <- ROCit::rocit(score, annotationData)
            if (AUC.abs) {
                aucs <- pred$AUC
            } else {
                aucs <- max(pred$AUC, 1 - pred$AUC)
            }
            ## Create data frame for With signature, P.value, AUC
            data.frame(Signature = i, P.value = round(pvals, 4),
                       neg10xP.value = round(neg10log, 4),
                       AUC = round(aucs, 4))
        }, SE_scored, annotationData)
        result <- do.call(rbind, sig_result) |>
            as.data.frame()
        row.names(result) <- NULL
        return(result)
    } else {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData, percent) {
            score <- col_info[i][, 1]
            ## Get lower and upper quantile
            lower <- (1 - percent) / 2
            upper <- 1 - lower
            ## Deal with PLAGE that have constant score (e.g. Sloot_HIV_2)
            if (length(unique(score)) == 1L) {
                sprintf(paste("Constant score found for siganture: %s,",
                              "results will be NA.\n"), i) |>
                    message()
                dat <- data.frame(i, NA, NA, NA, NA, NA)
                colnames(dat) <- c("Signature", "P.value", "neg10xP.value",
                                   "AUC",
                                   paste0("CI lower.", lower * 100, "%"),
                                   paste0("CI upper.", upper * 100, "%"))
                return(dat)
            }
            pvals <- stats::t.test(score ~ annotationData)$p.value
            neg10log <- -1 * log(pvals + 1e-4)
            pred <- ROCit::rocit(score, annotationData)
            if (AUC.abs) {
                aucs <- pred$AUC
            } else {
                aucs <- max(pred$AUC, 1 - pred$AUC)
            }
            ## Calculate bootstrapped AUC confidence interval.
            ## Repeated sampling scores and annotationData
            ## compute the AUC for the sampled pairs
            bootCI <- lapply(seq_len(num.boot),
                             function(j, score, annotationData) {
                index <- sample(seq_len(length(score)), replace = TRUE)
                tmp_score <- score[index]
                tmp_annotationData <- annotationData[index]
                ## Consider when re-sampling only has 1 cases, remove it
                if (length(unique(tmp_annotationData)) == 2L) {
                    tmp_pred <- ROCit::rocit(tmp_score, tmp_annotationData)
                    if (AUC.abs) {
                        tmp_auc <- tmp_pred$AUC
                    } else {
                        tmp_auc <- max(tmp_pred$AUC, 1 - tmp_pred$AUC)
                    }
                    tmp_auc
                } else {
                    NA
                }
            }, score, annotationData)
            bootCI <- unlist(bootCI) |>
                stats::na.omit()
            LowerAUC <- stats::quantile(bootCI, prob = lower, na.rm = TRUE)
            UpperAUC <- stats::quantile(bootCI, prob = upper, na.rm = TRUE)
            dat <- data.frame(i, round(pvals, 4),
                              round(neg10log, 4), round(aucs, 4),
                              round(LowerAUC, 4),
                              round(UpperAUC, 4))
            colnames(dat) <- c("Signature", "P.value",
                               "neg10xP.value", "AUC",
                               paste0("CI lower.", lower * 100, "%"),
                               paste0("CI upper.", upper * 100, "%"))
            dat
        }, SE_scored, annotationData, percent)
        result <- do.call(rbind, sig_result)
        row.names(result) <- NULL
        return(result)
    }
}

#' Check the input list format
#'
#' Function to check whether input list has desired format.
#' @param object_list A \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' @return Nothing is returned for this function.
.check_input <- function(object_list) {
    ## Check if the input is a list
    if (!is.list(object_list)) {
        paste("Function only supports a list of",
              "SummarizrdExperiment objetcs.") |>
            stop(call. = FALSE)
    }
    check_element_class <- vapply(object_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE)
    if (!all(check_element_class)) {
        msg <- sprintf("Elements(s) %s ",
                       paste0(which(!check_element_class), collapse = ", "))
        paste("Function only supports class: SummarizedExperiment.", msg,
              "in the list is/are not SummarizedExperiment object.") |>
            stop(call. = FALSE)
    }
    ## Check valid list names
    list_name <- names(object_list)
    if (is.null(list_name)) {
        paste("names of the input list should not be NULL.",
              "Add unique name for each element from the list.") |>
            stop(call. = FALSE)
    } else if (!is.na(match("", list_name))) {
        paste("Names of the input contains \"\".",
              "Replace  \"\" with unique character.") |>
            stop(call. = FALSE)
    }
}
