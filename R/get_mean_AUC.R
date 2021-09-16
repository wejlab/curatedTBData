#' Compute mean AUC value and bootstrapped confidence interval
#'   for multiple signature's mean AUC across studies
#' @name get_mean_AUC
#' @param df A data frame/matrix contains the interested numeric vector obtained
#'   from \code{\link[curatedTBData]{combine_auc}}.
#' @param colName A character string specifying the column name of the
#'   \code{data frame} for bootstrapping.
#' @param method A character string specifying the method used for
#'   computing bootstrap confidence interval.
#'   The choices are \code{c("percentile","empirical")}.
#'   The default is \code{"empirical"}.
#' @param num.boot Numeric. Number of bootstrap times.
#'   Default value is \code{100}.
#' @param percent A number indicates the percentage of confidence interval
#'   Default is value \code{0.95}.
#' @return A \code{data.frame} contains mean AUC,
#'   lower and upper bootstrap confidence interval
#'   for each gene signature across multiple studies.
#' @export
#' @examples
#' returned_resources <- curatedTBData(c("GSE107104", "GSE19435", "GSE19443"),
#'                                     dryrun = FALSE, curated.only = TRUE)
#' TBsignaturesSub <- TBSignatureProfiler::TBsignatures[1:5]
#' re1 <- lapply(returned_resources, function(x)
#'                    subset_curatedTBData(x, "TBStatus", c("Control","PTB")))
#' re2 <- lapply(re1, function(x)
#'          TBSignatureProfiler::runTBsigProfiler(input = x,
#'                                                useAssay = 1,
#'                                                signatures = TBsignaturesSub,
#'                                                algorithm = "ssGSEA",
#'                                                update_genes = FALSE))
#' df <- combine_auc(re2, annotationColName = "TBStatus",
#'                   signatureColNames = names(TBsignaturesSub),
#'                   num.boot = 100, percent = 0.95)
#' get_mean_AUC(df, colName = "AUC", method = "percentile",
#'              num.boot = 1000, percent = 0.95)
get_mean_AUC <- function(df, colName, method = c("percentile", "empirical"),
                         num.boot = 100, percent = 0.95) {
    ## Select signatures and associated AUC
    ## split them into list based on signature
    df_list <- df %>%
        dplyr::select("Signature", "AUC") %>%
        dplyr::group_split(.data$Signature)
    sigInfo <- df_list %>%
        base::vapply(function(x) as.character(x$Signature[1]), character(1))
    ## Get summarized table and
    ## bootstrap 95% Confidence Interval for the mean AUC
    sigInfo <- base::data.frame(Signature = sigInfo)
    meanAUC_list <- base::lapply(df_list, function(x)
        .bootstrap_mean_CI(x, colName, method, num.boot, percent))
    meanAUC <- base::do.call(base::rbind, meanAUC_list)
    return(base::cbind(sigInfo, meanAUC))
}
#' Compute bootstrapped confidence interval for single signature's mean AUC
#'   across multiple studies
#' @name .bootstrap_mean_CI
#' @inheritParams get_mean_AUC
#' @return A \code{data.frame} contains mean AUC,
#'   lower and upper bootstrap confidence interval
#'   for single gene signature across multiple studies.
.bootstrap_mean_CI <- function(df, colName,
                               method = c("percentile", "empirical"),
                               num.boot, percent = 0.95) {
    if (base::missing(method)) {
        method <- "empirical"
        base::paste("Missing method argument.",
                    "Using the default method: empirical") %>%
            base::message()
    }
    method <- base::match.arg(method)
    lower <- (1 - percent) / 2
    upper <- 1 - lower
    x <- base::unlist(df[, colName], use.names = FALSE)
    ## Remove NA's in PLAGE method
    x <- stats::na.omit(x)
    n <- length(x)
    if (n == 1L) {
        xbar <- x
        ci <- base::data.frame(base::round(xbar, 4), NA, NA)
        base::colnames(ci) <- c("MeanAUC",
                                base::paste0("CI lower.", lower * 100, "%"),
                                base::paste0("CI upper.", upper * 100, "%"))
        base::row.names(ci) <- NULL
        return(ci)
    }
    ## Sample mean
    xbar <- base::mean(x)
    ## Random re-samples from x
    bootstrapsample <- base::lapply(base::seq_len(num.boot), function(i)
        base::sample(x, n, replace = TRUE))
    bootstrapsample <- base::do.call(base::cbind, bootstrapsample)
    ## Compute the means x∗
    bsmeans <- base::colMeans(bootstrapsample)
    if (method == "empirical") {
        ## Compute deltastar for each bootstrap sample
        deltastar <- bsmeans - xbar
        # Find the 0.0.25 and 0.975 quantile for deltastar
        d <- stats::quantile(deltastar, c(lower, upper), na.rm = TRUE)
        ## Calculate the confidence interval for the mean.
        ci <- xbar - c(d[2], d[1])
    } else if (method == "percentile") {
        ci <- stats::quantile(bsmeans, c(lower, upper), na.rm = TRUE)
    }
    ## Set upper and lower bound for the confidence interval
    lower_ci <- base::round(ci[1], 4)
    lower_ci <- base::ifelse(lower_ci <= 0.5, 0.5, lower_ci)
    upper_ci <- base::round(ci[2], 4)
    upper_ci <- base::ifelse(upper_ci >= 1, 1, upper_ci)
    ci <- base::data.frame(base::round(xbar, 4), lower_ci, upper_ci)
    base::colnames(ci) <- c("MeanAUC",
                            base::paste0("CI lower.", lower * 100, "%"),
                            base::paste0("CI upper.", upper * 100, "%"))
    base::row.names(ci) <- NULL
    return(ci)
}
