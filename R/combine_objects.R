#' Merge samples with common gene names from selected studies
#' @name combine_objects
#' @param object_list A \code{list} of
#'   \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}
#'   or \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects. The object's assay contain expression data
#'   with probes mapped to gene symbol.
#'   \code{names(object_list)} should not be \code{NULL}.
#' @param experiment_name A character/vector of character to choose
#'   the name of the assay from the input \code{list} of object.
#' @param update_genes Boolean. Indicate whether update the gene symbols using
#'   \code{\link[HGNChelper]{checkGeneSymbols}}. Default is \code{TRUE}.
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object that contains combined data from the input.
#' @examples
#' geo <-  c("GSE19435", "GSE19439")
#' data_list <-  curatedTBData(geo, dryrun = FALSE, curated.only = TRUE)
#' combine_objects(data_list, experiment_name = "assay_curated")
#' @export
combine_objects <- function(object_list, experiment_name, update_genes = TRUE) {
    ## check the experiment_name argument
    if (missing(experiment_name)) {
        stop("Argument \"experiment_name\" is missing, with no default.")
    }
    ## check length of the list, should be greater than 1
    n <- length(object_list)
    if (n <= 1L) {
        sprintf("The length of the input list is %i,", n) |>
            paste("expecting more than 1 elements from the list") |>
            stop(call. = FALSE)
    }
    ## check names of the input object list
    obj_name <- names(object_list)
    if (is.null(obj_name)) {
        paste("Names of the input list should not be NULL.",
                    "Add unique name for each element in the list.") |>
            stop(call. = FALSE)
    } else if (!is.na(match("", obj_name))) {
        paste("Names of the input contains \"\".",
                    "Replace  \"\" with unique character.") |>
            stop(call. = FALSE)
    }
    ## Check whether it is a list of SummarizedExperiment or MultiAssayExperiment objects
    isSummarizedExperiment <- all(vapply(object_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE))
    isMultiAssayExperiement <- all(vapply(object_list, function(x)
        methods::is(x, "MultiAssayExperiment"), TRUE))
    if (isSummarizedExperiment) {
        dat_exprs_match <- .select_assay(object_list, experiment_name,
                                         Sobject = TRUE)
    } else if (isMultiAssayExperiement) {
        dat_exprs_match <- .select_assay(object_list, experiment_name,
                                         Sobject = FALSE)
    } else {
        paste("Input is not a list of MultiAssayExperiment",
                    "or SummarizedExperiment objetcs.") |>
            stop(call. = FALSE)
    }
    if (update_genes) {
        message("\"update_genes\" is TRUE, updating gene symbols")
        dat_exprs_match <- lapply(dat_exprs_match, update_gene_symbol)
    }
    ## Combine sample with common genes from a list of objects.
    ## Input data type should be data.frame
    dat_exprs_combine <- Reduce(function(x, y)
        merge(x, y, by = "id", all = FALSE),
        lapply(dat_exprs_match, function(x) {
           x$id <- row.names(x)
           x
    }))
    row_names <- dat_exprs_combine$id
    dat_exprs_count <- dat_exprs_combine |>
        dplyr::select(-.data$id) |>
        as.data.frame()
    row.names(dat_exprs_count) <- row_names
    ## Create combined column data information
    col_data <- lapply(seq_len(n), function(x) {
        col_data <- SummarizedExperiment::colData(object_list[[x]])
        col_data$Study <- names(object_list[x])
        as.data.frame(col_data)
    })
    ## Combine list into data frame with unequal columns
    ## fill in NA when columns from studies are not found
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
    col_info <- rbindx(col_data)
    Sample <- lapply(object_list, function(x)
        SummarizedExperiment::colData(x) |>
            row.names()) |>
            unlist(use.names = FALSE)
    row.names(col_info) <- Sample
    ## Remove samples that does not exist in the count
    index <- match(colnames(dat_exprs_count), Sample)  |>
        stats::na.omit()
    col_info <- col_info[index, ]
    ## Create output in the format of SummarizedExperiment
    result <- SummarizedExperiment::SummarizedExperiment(
        assays = list(assay1 = as.matrix(dat_exprs_count)),
        colData = col_info)
    return(result)
}

#' Select assay based on input list type
#'
#' @inheritParams combine_objects
#' @param Sobject Boolean. Indicate whether the input is a \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects.
#' @return A \code{list} of selected assays.
#'
.select_assay <- function(object_list, experiment_name, Sobject) {
    ## Merge code starts here
    n <- length(object_list)
    object_list_seq <- seq_len(n)
    object_list_names <- names(object_list)
    if (length(experiment_name) > 1L) {
        message("Found more than one \"experiment_name\".")
        ## experiment name for the list of object is different
        if (length(experiment_name) == n) {
            dat_exprs_match <- mapply(function(i, y) {
                x <- object_list[[i]]
                ## Avoid using ifelse, it deals with vectorized arguments.
                ## Returns same shape with the test.
                if (Sobject) {
                    dat_assay <- SummarizedExperiment::assays(x)[[y]]
                } else {
                    dat_assay <- MultiAssayExperiment::experiments(x)[[y]]
                }
                if (is.null(dat_assay)) {
                    sprintf("Object: %s with experiment name: %s",
                                  object_list_names[i], experiment_name) |>
                        paste("has assay NULL.") |>
                        stop(call. = FALSE)
                }
                as.data.frame(dat_assay)
            }, object_list_seq, experiment_name, SIMPLIFY = FALSE)
        } else {
            paste("Input list length",
                  "is different from the \"experiment_name\" vector.") |>
                stop(call. = FALSE)
        }
    } else {
        dat_exprs_match <- lapply(object_list_seq, function(i) {
            x <- object_list[[i]]
            ## Avoid using ifelse, it deals with vectorized arguments.
            ## Returns same shape with the test.
            if (Sobject) {
                dat_assay <- SummarizedExperiment::assays(x)[[experiment_name]]
            } else {
                dat_assay <- MultiAssayExperiment::experiments(x)[[experiment_name]]
            }
            if (is.null(dat_assay)) {
                sprintf("Object: %s with experiment name: %s",
                        object_list_names[i], experiment_name) |>
                    paste("has assay NULL.") |>
                    stop(call. = FALSE)
            }
            as.data.frame(dat_assay)
        })
    }
    names(dat_exprs_match) <- object_list_names
    return(dat_exprs_match)
}

#' Update gene names from input data
#' @name update_gene_symbol
#' @param dat_exprs A \code{data.frame} with row names as gene symbols to be updated.
#' @return A \code{data.frame} with updated gene symbol as row names.
#' @importFrom stats median na.pass
update_gene_symbol <- function(dat_exprs) {
    ## Function for updating gene names from HGNChelper::checkGeneSymbols
    update_genenames <- function(siglist) {
        newgenes <- HGNChelper::checkGeneSymbols(siglist,
                                                 unmapped.as.na = FALSE) |>
            suppressWarnings() |>
            suppressMessages()
        newgenes <- newgenes$Suggested.Symbol
        ind <- grep("//", newgenes)
        if (length(ind) != 0) {
            newgenes[ind] <- strsplit(newgenes[ind], " /// ")[[1]][1]
        }
        return(newgenes)
    }
    new_gene_names <- row.names(dat_exprs) |>
        update_genenames()
    new_gene_names_tab <- table(new_gene_names)
    ## Get genes with duplicates
    gene_names_dup <- names(new_gene_names_tab)[new_gene_names_tab > 1]
    if (!is.null(gene_names_dup)) {
        ## when we find duplicated gene names, we collapse gene symbol
        index <- which(new_gene_names %in% gene_names_dup)
        dat_exprs_no_duplicates <- dat_exprs[-index, ]
        row.names(dat_exprs_no_duplicates) <- new_gene_names[-index]
        dat_exprs_with_duplicates <- dat_exprs[index, ] |>
            as.data.frame() |>
            dplyr::mutate(SYMBOL = new_gene_names[index])
        exprs2 <- stats::aggregate(stats::as.formula(". ~ SYMBOL"),
                                   data = dat_exprs_with_duplicates,
                                   FUN = median, na.action = na.pass)
        row.names(exprs2) <- exprs2$SYMBOL
        dat_exprs_with_duplicates <- exprs2 |>
            dplyr::select(-.data$SYMBOL)
        dat_exprs <- rbind(dat_exprs_with_duplicates, dat_exprs_no_duplicates)
    } else {
        row.names(dat_exprs) <- new_gene_names
    }
    return(dat_exprs)
}
