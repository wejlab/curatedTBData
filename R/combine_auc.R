#' Obtain pvalue, emprirical AUC, and Bootstrap Confidence Interval for each signature.
#' @name .get_auc_stats
#' @param SE_scored A \code{SummarizedExperiment} object from TB signature profiling.
#' @param annotationColName A string specifying the feature of interest in the object's column data
#' @param signatureColNames A character/vector string contains name of gene signature.
#' @param num.boot Integer. Number of bootstrapping. If `NULL` (default), no confidence interval is computed.
#' @param percent Numeric. A number between 0 and 1, indicating the percentage of
#' confidence interval. Default is `0.95`.
#' @param AUC.abs Boolean. If AUC.abs = `TRUE`, return the AUC values from function \code{\link[ROCit]{rocit}}.
#' If AUC.abs = `FALSE`, return the AUC values for `max(AUC, 1-AUC)`.
#' @return A data frame/datatable contains p-value from two-sample t-test and AUC value for each signature.
.get_auc_stats <- function(SE_scored, annotationColName = "TBStatus", signatureColNames,
                           num.boot = NULL, percent=0.95, AUC.abs = FALSE) {
  # check signatureColNames
  index <- stats::na.omit(base::match(signatureColNames,
                                      base::colnames(SummarizedExperiment::colData(SE_scored))))
  signatureColNames <-  base::colnames(SummarizedExperiment::colData(SE_scored))[index]

  annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][, 1] %>%
    base::as.character() %>%
    base::as.factor()

  # get AUC value for each signature along with corresponding datasets
  if (base::is.null(num.boot)) {
    sig_result <- base::lapply(signatureColNames, function(i, SE_scored, annotationData) {
      score <- SummarizedExperiment::colData(SE_scored)[i][, 1] %>%
        base::as.vector()

      # Deal with scores that have constant value (e.g. Sloot_HIV_2)
      if (base::length(base::unique(score)) == 1) {
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
      # Create data frame for With signature, P.value, AUC
      base::data.frame(Signature = i, P.value = round(pvals, 4), AUC = round(aucs, 4))
    }, SE_scored, annotationData)

    result <- base::do.call(base::rbind, sig_result) %>%
      base::as.data.frame()
    base::row.names(result) <- NULL
    return(result)
  } else {
    sig_result <- base::lapply(signatureColNames,
                               function(i, SE_scored, annotationData, percent) {
                                 score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
                                 # Deal with PLAGE that have constant score (e.g. Sloot_HIV_2)
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
                                 # Get lower and upper quantile
                                 lower <- (1 - percent) / 2
                                 upper <- 1 - lower
                                 # Calculate bootstrapped AUC confidence interval
                                 # Repeated sampling scores and annotationData, compute the AUC for the sampled pairs
                                 bootCI <- base::lapply(base::seq_len(num.boot), function(j, score, annotationData) {
                                   index <- base::sample(base::seq_len(base::length(score)), replace = TRUE)
                                   tmp_score <- score[index]
                                   tmp_annotationData <- annotationData[index]
                                   # Consider when re-sampling only has 1 cases, remove it
                                   if (base::length(base::unique(tmp_annotationData)) == 2) {
                                     tmp_pred <- ROCit::rocit(tmp_score, tmp_annotationData)
                                     if (AUC.abs) {
                                       tmp_auc <- tmp_pred$AUC
                                     } else {
                                       tmp_auc <- base::max(tmp_pred$AUC, 1 - tmp_pred$AUC)
                                     }
                                     tmp_auc
                                   } else{
                                     NA
                                   }
                                 }, score, annotationData)

                                 bootCI <- base::unlist(bootCI) %>%
                                   stats::na.omit()
                                 LowerAUC <- stats::quantile(bootCI, prob = lower,
                                                             na.rm = TRUE)
                                 UpperAUC <- stats::quantile(bootCI, prob = upper,
                                                             na.rm = TRUE)
                                 dat <- base::data.frame(i, base::round(pvals, 4),
                                                         base::round(neg10log, 4),
                                                         base::round(aucs, 4),
                                                         base::round(LowerAUC, 4),
                                                         base::round(UpperAUC, 4))
                                 base::colnames(dat) <- c("Signature", "P.value",
                                                          "neg10xP.value", "AUC",
                                                          base::paste0("CI lower.", lower * 100, "%"),
                                                          base::paste0("CI upper.", upper * 100, "%"))
                                 dat
                               }, SE_scored, annotationData, percent)
    result <- base::do.call(base::rbind, sig_result)
    base::row.names(result) <- NULL
    return(result)
  }
}

#' Combine results from the `list` of [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class] object.
#'
#' Calculate p-value, AUC values, and bootstrapped confidence interval for the AUC.
#' @name combine_auc
#' @param SE_scored_list A `list` of [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#' objects from \code{\link[TBSignatureProfiler]{runTBsigProfiler}}.
#' @param annotationColName A string specifying the feature of interest in the object's column data
#' @param signatureColNames A character/vector string contains name of gene signature.
#' @param num.boot Integer. Number of bootstrapping. If `NULL` (default), no confidence interval is computed.
#' @param percent Numeric. A number between 0 and 1, indicating the percentage of
#' confidence interval. Default is `0.95`.
#' @param AUC.abs Boolean. If AUC.abs = `TRUE`, return the AUC values from function \code{\link[ROCit]{rocit}}.
#' If AUC.abs = `FALSE`, return the AUC values for `max(AUC, 1-AUC)`.
#' @param BPPARAM An instance inherited from \code{bplappy}.
#' See \code{\link[BiocParallel]{bplapply}} for details.
#' @return A data frame with features including Signatures, P.value, neg10xLog(P.value)
#' and AUC for each signature across studies.
#' @examples
#' data(TB_indian, package = "TBSignatureProfiler")
#' TBsignaturesSub <- TBSignatureProfiler::TBsignatures[1:5]
#' res <- TBSignatureProfiler::runTBsigProfiler(input = TB_indian, useAssay = "logcounts",
#'                                              signatures = TBsignaturesSub, algorithm = "ssGSEA",
#'                                              combineSigAndAlgorithm = TRUE)
#' re <- combine_auc(list(TB_indian = res), annotationColName = "label",
#'                   signatureColNames = names(TBsignaturesSub), num.boot = 100, percent = 0.95)
#' re
#' @export
combine_auc <- function(SE_scored_list, annotationColName, signatureColNames,
                        num.boot = NULL, percent = 0.95, AUC.abs = FALSE,
                        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
  param <- BPPARAM
  # Check if the input is a list
  if (base::class(SE_scored_list)[1] != "list") {
    base::stop(base::sprintf("Function only supports a list of SummarizedExperiment. The input class: %s.",
                             base::class(SE_scored_list)[1]))
  }
  SE_scored_list_class <- base::class(SE_scored_list[[1]])[1]
  if (SE_scored_list_class != "SummarizedExperiment") {
    base::stop(base::sprintf("Function only supports SummarizedExperiment within the list. The input class: %s. ",
                             SE_scored_list_class))
  }
  list_name <- base::names(SE_scored_list)
  if (base::is.null(list_name)) {
    # Only 1 study with NULL names within the list, make the names as Study1
    base::names(SE_scored_list) <- "Study1"
  } else if (!base::is.na(base::match("", list_name))) {
    base::stop(base::sprintf("Names of the input contains \"\". Replace \"\" with a non-empty string."))
  }
  aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x) {
    .get_auc_stats(x, annotationColName, signatureColNames, num.boot, percent, AUC.abs)
  }, BPPARAM = param)
  aucs_result_dat <- base::do.call(base::rbind, aucs_result)

  # re-order data based on their median AUC
  # Remove NA value
  aucs_result_dat1 <- stats::na.omit(aucs_result_dat)
  aucs_result_dat_median <- aucs_result_dat1 %>%
    dplyr::group_by(.data$Signature) %>%
    dplyr::summarise_all(stats::median) %>%
    dplyr::arrange(dplyr::desc(.data$AUC))

  # Order signatures based on median AUC values
  Signature_order <- base::as.character(aucs_result_dat_median$Signature)
  # Re-order gene siganture, re-level
  # this step is to let ridge plot ordered based on median value
  aucs_result_dat$Signature <- base::factor(aucs_result_dat$Signature,
                                            levels = Signature_order)
  # label name of each dataset under column "Study"
  aucs_result_dat$Study <- base::gsub("\\..*", "", base::row.names(aucs_result_dat))
  base::row.names(aucs_result_dat) <- NULL
  return(aucs_result_dat)
}

#' Compute bootstrapped confidence interval for signature's mean AUC across studies
#' @name bootstrap_mean_CI
#' @param data A data frame/matrix contains the interested numeric vector.
#' @param colName A character string specifying the column name of the data frame
#' for bootstrapping.
#' @param percent A number indicates the percentage of confidence interval.
#' @param method A character string specifying the method used for computing bootstrap confidence interval.
#' The choices are c("percentile","empirical"). The default is "empirical".
#' @param num.boot Numeric. Number of bootstrap times.
#' @return A data frame with lower and upper bootstrap confidence interval.
#' @export
#'
bootstrap_mean_CI <- function(data, colName, percent = 0.95,
                              method = c("percentile", "empirical"), num.boot) {

  if (base::missing(method)) {
    method <- "empirical"
    base::message(base::sprintf("Missing method argument. The default method used for bootstrap confidence interval is %s",
                          method))
  }
  method <- base::match.arg(method)
  lower <- (1 - percent) / 2
  upper <- 1 - lower

  x <- base::unlist(data[, colName], use.names = FALSE)
  x <- stats::na.omit(x) # Remove NA's in PLAGE method

  n <- length(x)
  if (n == 1) {
    xbar <- x
    ci <- base::data.frame(base::round(xbar, 4), NA, NA)
    base::colnames(ci) <- c("MeanAUC", base::paste0("CI lower.", lower * 100, "%"),
                            base::paste0("CI upper.", upper * 100, "%"))
    base::row.names(ci) <- NULL
    return(ci)
  }
  # sample mean
  xbar <- base::mean(x)
  # random resamples from x
  bootstrapsample <- base::lapply(base::seq_len(num.boot), function(i)
    base::sample(x, n, replace = TRUE))
  bootstrapsample <- base::do.call(base::cbind, bootstrapsample)

  # Compute the means x∗
  bsmeans <-  base::colMeans(bootstrapsample)

  if (method == "empirical") {
    # Compute deltastar for each bootstrap sample
    deltastar <-  bsmeans - xbar
    # Find the 0.0.25 and 0.975 quantile for deltastar
    d <-  stats::quantile(deltastar, c(lower, upper), na.rm = TRUE)
    # Calculate the confidence interval for the mean.
    ci  <-  xbar - c(d[2], d[1])
  } else if (method == "percentile") {
    ci <- stats::quantile(bsmeans, c(lower, upper), na.rm = TRUE)
  }
  # Set upper and lower bound for the confidence interval
  lower_ci <- base::round(ci[1], 4)
  lower_ci<- base::ifelse(lower_ci <= 0.5, 0.5, lower_ci)
  upper_ci <- base::round(ci[2], 4)
  upper_ci <- base::ifelse(upper_ci >= 1, 1, upper_ci)
  ci <- base::data.frame(base::round(xbar, 4), lower_ci, upper_ci)
  base::colnames(ci) <- c("MeanAUC", base::paste0("CI lower.", lower * 100, "%"),
                          base::paste0("CI upper.", upper * 100, "%"))
  base::row.names(ci) <- NULL
  return(ci)
}
