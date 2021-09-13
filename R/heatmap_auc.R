#' Plot heatmap for multiple signatures across datasets
#'
#' Obtain heatmap plots for empirical AUC distribution of gene signatures across multiple studies
#' @name heatmap_auc
#' @param combine_dat A \code{data.frame} contains at least three columns with name \code{Signatures},
#' \code{Study}, and \code{AUC} respectively. Usually the output from \code{\link[curatedTBData]{combine_auc}}
#' The format of the signature name is expected as: Name_SignatureType_Number (i.e. Anderson_OD_42).
#' See \code{names(TBSignatureProfiler::TBsignatures)} for examples
#' @param GSE_sig A \code{data.frame} contains information about each signature and
#' its training/discovery dataset(s) name. Default is \code{NULL}
#' @param facet Boolean. If \code{TRUE}, grouping signatures into clusters based on their type.
#' If \code{FALSE}, output without grouping. Default is \code{TRUE}
#' @param clustering Boolean. If \code{TRUE}, perform clustering using hierarchical clustering method.
#' If \code{FALSE}, output without clustering. Default is \code{TRUE}
#' @return A heatmap shows the performance of signature's across multiple studies
#' filled with AUC values
#' @export
#' @examples
#' combine_dat_exp <- data.frame(Signature = rep(c("Anderson_42", "Anderson_OD_51",
#'                               "Berry_393", "Berry_OD_86", "Blankley_5"), 2),
#'                               AUC = runif(10, 0.5, 1),
#'                               Study = rep(c("GSE39939", "GSE19442"), each = 5))
#' GSE_sig_exp <- data.frame(TBSignature = c("Anderson", "Anderson", "Berry", "Berry"),
#'                           Study = c("GSE39939", "GSE39940", "GSE19442", "GSE19443"))
#' heatmap_auc(combine_dat_exp, GSE_sig_exp, facet = FALSE)
#' heatmap_auc(combine_dat_exp, GSE_sig_exp, facet = TRUE)
heatmap_auc <- function(combine_dat, GSE_sig = NULL, facet = TRUE, clustering = TRUE) {
    ## Subset input data.frame with the desired column names check whether column contains 'Signature', 'Study', 'AUC'
    expect_name <- c("Signature", "Study", "AUC")
    index_name <- base::match(expect_name, base::colnames(combine_dat))
    if (base::any(base::is.na(index_name))) {
        base::stop(base::sprintf("Column with name(s): %s is/are missing.",
                                 paste0(expect_name[base::is.na(index_name)],
                                        collapse = ", ")))
    }
    dat <- combine_dat[, index_name]
    if (!base::is.factor(dat$Study)) {
        dat$Study <- base::as.factor(dat$Study)
    }
    if (base::length(base::unique(dat$Study)) > 1L && clustering == TRUE) {
        ## Clustering AUC values if the number of studies is greater than 1 Transform form long to wide data: first column is
        ## the study names and column names is signatures This step is necessary for clustering
        data_wide <- reshape2::dcast(dat, stats::formula("Study ~ Signature"))
        base::row.names(data_wide) <- data_wide$Study
        # remove study name column
        dat_input <- base::as.matrix(data_wide[, -1])
        dat_input[base::is.na(dat_input)] <- NA
        dd <- stats::dist(dat_input)
        hc <- stats::hclust(dd)
        dat_input <- dat_input[hc$order, ]
        ## Get mean AUC for each study across multiple gene signatures
        datta <- base::cbind(dat_input,
                             Avg = base::rowMeans(dat_input, na.rm = TRUE)) %>% ## Transform back into long format
            reshape2::melt()
    } else {
        datta <- base::data.frame(Var1 = dat$Study, Var2 = dat$Signature,
                                  value = dat$AUC)
        mean_df <- dat %>%
            dplyr::group_by(.data$Study) %>%
            dplyr::summarise(value = mean(.data$AUC))
        mean_df <- base::data.frame(Var1 = mean_df$Study, Var2 = "Avg",
                                    value = mean_df$value)
        datta <- base::rbind(datta, mean_df)
    }
    ## Get training data position index
    datta$trian <- FALSE
    index <- NULL
    if (!base::is.null(GSE_sig)) {
        GSE_sig <- .expand_study(GSE_sig)
        for (i in base::seq_len(base::nrow(GSE_sig))) {
            kk <- datta[base::grep(GSE_sig$TBSignature[i], datta$Var2), ]
            kk$indx <- base::row.names(kk)
            indx <- kk[base::which(base::as.character(kk$Var1) %in%
                                       GSE_sig$Study[i]), "indx"]
            index <- c(index, indx)
        }
    }
    ## Label signature type based on the input signatures
    signatureColNames <- base::as.character(base::unique(dat$Signature))
    sig_type_temp <- base::vapply(base::strsplit(signatureColNames, "_"),
                                  function(x) x[2], base::character(1))
    sig_type_index <- sig_type_temp %>%
        base::as.numeric() %>%
        base::is.na() %>%
        base::suppressWarnings()
    ## Get signature type
    sig_type <- sig_type_temp[sig_type_index] %>%
        base::unique()
    ## Assign category: Disease for those do not have signature type: e.g. Anderson_42
    datta$sig_typek <- "Disease"
    for (i in sig_type) {
        datta$sig_typek[grep(i, datta$Var2)] <- i
    }
    datta$sig_typek[base::grep("Avg", datta$Var2)] <- "Avg"
    datta$sig_typek <- base::factor(datta$sig_typek,
                                    levels = c("Avg", sig_type, "Disease"))
    datta[base::as.numeric(index), "trian"] <- TRUE
    ## Subset datta with training study and its associated signature(s)
    frames <- datta[datta$trian, c("Var1", "Var2", "sig_typek")]
    p <- ggplot2::ggplot(data = datta,
                         ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                      fill = .data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_distiller(palette = "RdPu", trans = "reverse") +
        ggplot2::geom_text(ggplot2::aes(label = base::round(.data$value, 2)), cex = 3.5)
    if (facet) {
        p <- p + ggplot2::facet_grid(.data$sig_typek ~ ., switch = "y", scales = "free", space = "free")
        frame_facet <- .facet_rect_position(datta, frames)
        if (!base::nrow(frame_facet) == 0L) {
            p <- p + ggplot2::geom_rect(data = frame_facet,
                                        ggplot2::aes(xmin = .data$Var1 - 0.5,
                                                     xmax = .data$Var1 + 0.5,
                                                     ymin = .data$Var2 -  0.5,
                                                     ymax = .data$Var2 + 0.5),
                                        size = 1, fill = NA, colour = "black",
                                        inherit.aes = FALSE)
        }
    } else {
        frames$Var1 <- base::as.integer(frames$Var1)
        frames$Var2 <- base::as.integer(frames$Var2)
        p <- p + ggplot2::geom_rect(data = frames,
                                    ggplot2::aes(xmin = .data$Var1 - 0.5,
                                                 xmax = .data$Var1 + 0.5,
                                                 ymin = .data$Var2 - 0.5,
                                                 ymax = .data$Var2 + 0.5),
                                    size = 1, fill = NA, colour = "black")
    }
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_text(angle = 45,
                                                                vjust = 1,
                                                                size = 12,
                                                                hjust = 1),
                            axis.text.y = ggplot2::element_text(size = 12))
    return(p)
}

#' Obtain training index in facet grid
#' @param datta A \code{data.frame} with column names indicating signatures, study,
#' value, and signature type
#' @param frames A \code{data.frame} with column names indicating signatures, study,
#' value, and signature type for training/discovery studies
#' @return A \code{data.frame} with position index on the x-axis and y-axis, and signature type.
#' The result is mainly used as input for \code{\link[ggplot2]{geom_rect}}
.facet_rect_position <- function(datta, frames) {
    # Split data frame into list based on different signature type
    frames_list <- frames %>%
        dplyr::group_split(.data$sig_typek)
    base::names(frames_list) <- base::lapply(frames_list, function(x) x$sig_typek[1]) %>%
        base::unlist()
    datta_list <- datta %>%
        dplyr::group_split(.data$sig_typek)
    base::names(datta_list) <- base::lapply(datta_list, function(x) x$sig_typek[1]) %>%
        base::unlist()
    ## Get the correct index in for training dataset change sig_type levels from sub list based on characters in the full list
    frame_facet1 <- base::lapply(base::names(frames_list), function(i) {
        num_Var1 <- frames_list[[i]]$Var1 %>%
            base::as.integer()
        frame_sig <- frames_list[[i]]$Var2
        num_Var2 <- base::factor(frame_sig, levels = base::unique(frame_sig)) %>%
            base::as.integer()
        frames_list[[i]] %>%
            dplyr::mutate(Var1 = num_Var1, Var2 = num_Var2)
    })
    re <- base::do.call(base::rbind, frame_facet1) %>%
        base::as.data.frame()
    return(re)
}

#' Expand study section for \code{SignatureInfoTraining}
#' @param GSE_sig A \code{data.frame} contains information about each signature and
#' its training/discovery dataset(s) name. Default is \code{NULL}
#' @return A expand \code{data.frame} for gene signatures and dataset
.expand_study <- function(GSE_sig) {
    n <- base::nrow(GSE_sig)
    col_name <- base::colnames(GSE_sig)
    data_list <- base::lapply(base::seq_len(n), function(i) {
        study_vector <- base::strsplit(GSE_sig[, col_name[2]][i], split = "&")
        df <- base::data.frame(GSE_sig[, col_name[1]][i], study_vector)
        colnames(df) <- col_name
        df
    })
    re <- base::do.call(rbind, data_list) %>%
        base::as.data.frame()
    return(re)
}
