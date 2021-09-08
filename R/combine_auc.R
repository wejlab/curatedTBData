#' Combine results from a \code{list} of \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} objects
#'
#' Calculate p-value, AUC values, and bootstrapped confidence interval for TB gene signatures across multiple studies
#' @name combine_auc
#' @param SE_scored_list A \code{list} of \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' objects obtained from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}
#' @param annotationColName A string specifying the feature of interest in the object's column data
#' @param signatureColNames A character/vector string contains name of gene signature
#' @param num.boot Integer. Number of bootstrapping. If \code{NULL} (default), no confidence interval is computed
#' @param percent Numeric. A number between 0 and 1, indicating the percentage of
#' confidence interval. Default is \code{0.95}
#' @param AUC.abs Boolean. If \code{AUC.abs = TRUE}, compute the AUC values based on
#' \code{\link[ROCit]{rocit}}. If \code{AUC.abs = FALSE}, compute the AUC values for \code{max(AUC, 1-AUC)}
#' @param BPPARAM An instance inherited from \code{\link[BiocParallel]{bplapply}}
#' @return A \code{data.frame} with features including Signatures, P.value, neg10xLog(P.value) and AUC for each signature across studies
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
    param <- BPPARAM
    ## Check if the input is a list
    if (!methods::is(SE_scored_list, "list")) {
        base::stop(base::sprintf("Function only supports a list of SummarizedExperiment. The input class: %s.",
                                 base::class(SE_scored_list)[1]))
    }
    check_element_class <- base::vapply(SE_scored_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE)
    if (!base::all(check_element_class)) {
        base::stop(base::sprintf("Function only supports class: SummarizedExperiment within the list. Element(s) %s in the list is/are not SummarizedExperiment object.",
                                 paste0(which(!check_element_class), collapse = ", ")))
    }
    list_name <- base::names(SE_scored_list)
    if (base::is.null(list_name)) {
        ## Only 1 study with NULL names within the list, make the names as Study1
        base::names(SE_scored_list) <- "Study1"
    } else if (!base::is.na(base::match("", list_name))) {
        base::stop(base::sprintf("Names of the input contains \"\". Replace \"\" with a non-empty string."))
    }
    aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x) {
        .get_auc_stats(x, annotationColName, signatureColNames, num.boot,
                       percent, AUC.abs)
    }, BPPARAM = param)
    aucs_result_dat <- base::do.call(base::rbind, aucs_result)
    ## Re-order data based on their median AUC Remove NA value
    aucs_result_dat1 <- stats::na.omit(aucs_result_dat)
    aucs_result_dat_median <- aucs_result_dat1 %>%
        dplyr::group_by(.data$Signature) %>%
        dplyr::summarise_all(stats::median) %>%
        dplyr::arrange(dplyr::desc(.data$AUC))

    ## Order signatures based on median AUC values
    Signature_order <- base::as.character(aucs_result_dat_median$Signature)
    ## Re-order gene signature, re-level this step is to let ridge plot ordered based on median value
    aucs_result_dat$Signature <- base::factor(aucs_result_dat$Signature,
                                              levels = Signature_order)
    ## Label name of each dataset under column 'Study'
    aucs_result_dat$Study <- base::gsub("\\..*", "",
                                        base::row.names(aucs_result_dat))
    base::row.names(aucs_result_dat) <- NULL
    return(aucs_result_dat)
}

#' Obtain pvalue, emprirical AUC, and Bootstrap Confidence Interval for each signature
#' @name .get_auc_stats
#' @param SE_scored A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object from TB signature profiling
#' @param annotationColName A string specifying the feature of interest in the object's column data
#' @param signatureColNames A character/vector string contains name of gene signature
#' @param num.boot Integer. Number of bootstrapping
#' @param percent Numeric. A number between 0 and 1, indicating the percentage of
#' confidence interval
#' @param AUC.abs Boolean. If \code{AUC.abs = TRUE}, compute the AUC values from function \code{\link[ROCit]{rocit}}
#' If \code{AUC.abs = FALSE}, compute the AUC values for \code{max(AUC, 1-AUC)}
#' @return A \code{data.frame} contains p-value from two-sample t-test and AUC value for each signature
.get_auc_stats <- function(SE_scored, annotationColName, signatureColNames,
                           num.boot, percent, AUC.abs) {
    ## Check signatureColNames
    index <- stats::na.omit(base::match(signatureColNames,
                                        base::colnames(SummarizedExperiment::colData(SE_scored))))
    signatureColNames <- base::colnames(SummarizedExperiment::colData(SE_scored))[index]
    annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][, 1] %>%
        base::as.character() %>%
        base::as.factor()
    ## Check levels of annotationData
    anno_level_len <- base::length(base::unique(annotationData))
    if (anno_level_len != 2L) {
        base:: stop("Annotation data should have exactly two levels. Your level: %i",
                    anno_level_len)
    }
    ## Get AUC value for each signature along with corresponding datasets
    if (base::is.null(num.boot)) {
        sig_result <- base::lapply(signatureColNames, function(i, SE_scored, annotationData) {
            score <- SummarizedExperiment::colData(SE_scored)[i][, 1] %>%
                base::as.vector()
            ## Deal with scores that have constant value (e.g. Sloot_HIV_2)
            if (base::length(base::unique(score)) == 1L) {
                dat <- base::data.frame(Signature = i, P.value = NA, AUC = NA)
                return(dat)
            }
            pvals <- stats::t.test(score ~ annotationData)$p.value
            pred <- ROCit::rocit(score, annotationData)
            if (AUC.abs) {
                aucs <- pred$AUC
            } else {
                aucs <- base::max(pred$AUC, 1 - pred$AUC)
            }
            ## Create data frame for With signature, P.value, AUC
            base::data.frame(Signature = i, P.value = round(pvals, 4),
                             AUC = round(aucs, 4))
        }, SE_scored, annotationData)
        result <- base::do.call(base::rbind, sig_result) %>%
            base::as.data.frame()
        base::row.names(result) <- NULL
        return(result)
    } else {
        sig_result <- base::lapply(signatureColNames, function(i, SE_scored, annotationData, percent) {
            score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
            ## Deal with PLAGE that have constant score (e.g. Sloot_HIV_2)
            if (base::length(base::unique(score)) == 1) {
                dat <- base::data.frame(i, NA, NA, NA, NA)
                base::colnames(dat) <- c("Signature", "P.value", "AUC",
                                         paste0("CI lower.", lower * 100, "%"),
                                         paste0("CI upper.", upper * 100, "%"))
                return(dat)
            }
            pvals <- stats::t.test(score ~ annotationData)$p.value
            neg10log <- -1 * log(pvals)
            pred <- ROCit::rocit(score, annotationData)
            if (AUC.abs) {
                aucs <- pred$AUC
            } else {
                aucs <- base::max(pred$AUC, 1 - pred$AUC)
            }
            ## Get lower and upper quantile
            lower <- (1 - percent) / 2
            upper <- 1 - lower
            ## Calculate bootstrapped AUC confidence interval Repeated sampling scores and annotationData, compute the AUC for
            ## the sampled pairs
            bootCI <- base::lapply(base::seq_len(num.boot), function(j, score, annotationData) {
                index <- base::sample(base::seq_len(base::length(score)), replace = TRUE)
                tmp_score <- score[index]
                tmp_annotationData <- annotationData[index]
                ## Consider when re-sampling only has 1 cases, remove it
                if (base::length(base::unique(tmp_annotationData)) == 2) {
                    tmp_pred <- ROCit::rocit(tmp_score, tmp_annotationData)
                    if (AUC.abs) {
                        tmp_auc <- tmp_pred$AUC
                    } else {
                        tmp_auc <- base::max(tmp_pred$AUC, 1 - tmp_pred$AUC)
                    }
                    tmp_auc
                } else {
                    NA
                }
            }, score, annotationData)
            bootCI <- base::unlist(bootCI) %>%
                stats::na.omit()
            LowerAUC <- stats::quantile(bootCI, prob = lower, na.rm = TRUE)
            UpperAUC <- stats::quantile(bootCI, prob = upper, na.rm = TRUE)
            dat <- base::data.frame(i,
                                    base::round(pvals, 4), base::round(neg10log, 4),
                                    base::round(aucs, 4), base::round(LowerAUC, 4),
                                    base::round(UpperAUC, 4))
            base::colnames(dat) <- c("Signature", "P.value", "neg10xP.value", "AUC",
                                     base::paste0("CI lower.", lower * 100, "%"),
                                     base::paste0("CI upper.", upper * 100, "%"))
            dat
        }, SE_scored, annotationData, percent)
        result <- base::do.call(base::rbind, sig_result)
        base::row.names(result) <- NULL
        return(result)
    }
}
