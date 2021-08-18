#' Obtain ridge plots for empirical AUC distribution for signature scores.
#' @name ridgeplot_auc
#' @param aucs_result A `data.frame` contains signatures, p-value, and AUC.
#' Usually the output from \code{\link{combine_auc}}.
#' @return Ridge plot with median line
#' @examples
#' aucs_result <- data.frame(Signature = c("Anderson_42", "Anderson_OD_51", "Berry_393"),
#'                           AUC = stats::runif(3,0,0.5))
#' p_ridge <- ridgeplot_auc(aucs_result)
#' @export
ridgeplot_auc <- function(aucs_result) {
  # add 50% AUC line
  aucs_result_dat_lines <- data.frame(Signature = aucs_result$Signature,x0=0.5)

  p_ridge <- ggplot2::ggplot(aucs_result, ggplot2::aes(x = .data$AUC,
                                                       y = .data$Signature)) +
    ggridges::geom_density_ridges(jittered_points = TRUE,alpha = 0.7,
                                  quantile_lines = TRUE, quantiles = 2) +
    ggplot2::geom_segment(data = aucs_result_dat_lines,
                          ggplot2::aes(x = .data$x0,
                                       xend = .data$x0,
                                       y = as.numeric(.data$Signature),
                                       yend = as.numeric(.data$Signature) + .9),
                          color = "red")
  return(p_ridge)
}
