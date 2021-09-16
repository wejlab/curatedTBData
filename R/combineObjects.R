#' Merge samples with common gene names from selected studies
#' @name combineObjects
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
#' combineObjects(data_list, experiment_name = "assay_curated")
#' @export
combineObjects <- function(object_list, experiment_name, update_genes = TRUE) {
    ## check the experiment_name argument
    if (base::missing(experiment_name)) {
        base::stop("Argument \"experiment_name\" is missing, with no default.")
    }
    ## check length of the list, should be greater than 1
    n <- base::length(object_list)
    if (n <= 1L) {
        base::sprintf("The length of the input list is %i,", n) %>%
            base::paste("expecting more than 1 elements from the list") %>%
            base::stop(call. = FALSE)
    }
    ## check names of the input object list
    obj_name <- base::names(object_list)
    if (base::is.null(obj_name)) {
        base::paste("Names of the input list should not be NULL.",
                    "Add unique name for each element in the list.") %>%
            base::stop(call. = FALSE)
    } else if (!base::is.na(base::match("", obj_name))) {
        base::paste("Names of the input contains \"\".",
                    "Replace  \"\" with unique character.") %>%
            base::stop(call. = FALSE)
    }
    ## Check whether it is a list of SummarizedExperiment or MultiAssayExperiment objects
    isSummarizedExperiment <- base::all(base::vapply(object_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE))
    isMultiAssayExperiement <- base::all(base::vapply(object_list, function(x)
        methods::is(x, "MultiAssayExperiment"), TRUE))
    if (isSummarizedExperiment) {
        dat_exprs_match <- .select_assay(object_list, experiment_name,
                                         Sobject = TRUE)
    } else if (isMultiAssayExperiement) {
        dat_exprs_match <- .select_assay(object_list, experiment_name,
                                         Sobject = FALSE)
    } else {
        base::paste("Input is not a list of MultiAssayExperiment",
                    "or SummarizedExperiment objetcs.") %>%
            base::stop(call. = FALSE)
    }
    if (update_genes) {
        base::message("\"update_genes\" is TRUE, updating gene symbols")
        dat_exprs_match <- base::lapply(dat_exprs_match, update_gene_symbol)
    }
    ## Combine sample with common genes from a list of objects.
    ## Input data type should be data.frame
    dat_exprs_combine <- base::Reduce(function(x, y)
        base::merge(x, y, by = "id", all = FALSE),
        base::lapply(dat_exprs_match, function(x) {
           x$id <- base::row.names(x)
           x
    }))
    row_names <- dat_exprs_combine$id
    dat_exprs_count <- dat_exprs_combine %>%
        dplyr::select(-.data$id) %>%
        base::as.data.frame()
    base::row.names(dat_exprs_count) <- row_names
    ## Create combined column data information
    col_data <- base::lapply(base::seq_len(n), function(x) {
        col_data <- SummarizedExperiment::colData(object_list[[x]])
        col_data$Study <- base::names(object_list[x])
        base::as.data.frame(col_data)
    })
    ## Combine list into data frame with unequal columns
    ## fill in NA when columns from studies are not found
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
    col_info <- rbindx(col_data)
    Sample <- base::lapply(object_list, function(x)
        SummarizedExperiment::colData(x) %>%
            base::row.names()) %>%
        base::unlist(use.names = FALSE)
    base::row.names(col_info) <- Sample
    ## Remove samples that does not exist in the count
    index <- base::match(base::colnames(dat_exprs_count), Sample) %>%
        stats::na.omit()
    col_info <- col_info[index, ]
    ## Create output in the format of SummarizedExperiment
    result <- SummarizedExperiment::SummarizedExperiment(
        assays = base::list(assay1 = base::as.matrix(dat_exprs_count)),
        colData = col_info)
    return(result)
}

#' Select assay based on input list type
#'
#' @inheritParams combineObjects
#' @param Sobject Boolean. Indicate whether the input is a \code{list} of
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   objects.
#' @return A \code{list} of selected assays.
#'
.select_assay <- function(object_list, experiment_name, Sobject) {
    ## Merge code starts here
    n <- base::length(object_list)
    object_list_seq <- base::seq_len(n)
    object_list_names <- base::names(object_list)
    if (base::length(experiment_name) > 1L) {
        base::message("Found more than one \"experiment_name\".")
        ## experiment name for the list of object is different
        if (base::length(experiment_name) == n) {
            dat_exprs_match <- base::mapply(function(i, y) {
                x <- object_list[[i]]
                ## Avoid using ifelse, it deals with vectorized arguments.
                ## Returns same shape with the test.
                if (Sobject) {
                    dat_assay <- SummarizedExperiment::assays(x)[[y]]
                } else {
                    dat_assay <- MultiAssayExperiment::experiments(x)[[y]]
                }
                if (base::is.null(dat_assay)) {
                    base::sprintf("Object: %s with experiment name: %s",
                                  object_list_names[i], experiment_name) %>%
                        base::paste("has assay NULL.") %>%
                        base::stop(call. = FALSE)
                }
                base::as.data.frame(dat_assay)
            }, object_list_seq, experiment_name)
        } else {
            base::paste("Input list length",
                        "is different from the \"experiment_name\" vector.") %>%
                base::stop(call. = FALSE)
        }
    } else {
        dat_exprs_match <- base::lapply(object_list_seq, function(i) {
            x <- object_list[[i]]
            ## Avoid using ifelse, it deals with vectorized arguments.
            ## Returns same shape with the test.
            if (Sobject) {
                dat_assay <- SummarizedExperiment::assays(x)[[experiment_name]]
            } else {
                dat_assay <- MultiAssayExperiment::experiments(x)[[experiment_name]]
            }
            if (base::is.null(dat_assay)) {
                base::sprintf("Object: %s with experiment name: %s",
                              object_list_names[i], experiment_name) %>%
                    base::paste("has assay NULL.") %>%
                    base::stop(call. = FALSE)
            }
            base::as.data.frame(dat_assay)
        })
    }
    base::names(dat_exprs_match) <- object_list_names
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
                                                 unmapped.as.na = FALSE) %>%
            base::suppressWarnings() %>%
            base::suppressMessages()
        newgenes <- newgenes$Suggested.Symbol
        ind <- base::grep("//", newgenes)
        if (base::length(ind) != 0) {
            newgenes[ind] <- base::strsplit(newgenes[ind], " /// ")[[1]][1]
        }
        return(newgenes)
    }
    new_gene_names <- base::row.names(dat_exprs) %>%
        update_genenames()
    new_gene_names_tab <- base::table(new_gene_names)
    ## Get genes with duplicates
    gene_names_dup <- base::names(new_gene_names_tab)[new_gene_names_tab > 1]
    if (!base::is.null(gene_names_dup)) {
        ## when we find duplicated gene names, we collapse gene symbol
        index <- base::which(new_gene_names %in% gene_names_dup)
        dat_exprs_no_duplicates <- dat_exprs[-index, ]
        base::row.names(dat_exprs_no_duplicates) <- new_gene_names[-index]
        dat_exprs_with_duplicates <- dat_exprs[index, ] %>%
            base::as.data.frame() %>%
            dplyr::mutate(SYMBOL = new_gene_names[index])
        exprs2 <- stats::aggregate(stats::as.formula(". ~ SYMBOL"),
                                   data = dat_exprs_with_duplicates,
                                   FUN = median, na.action = na.pass)
        base::row.names(exprs2) <- exprs2$SYMBOL
        dat_exprs_with_duplicates <- exprs2 %>%
            dplyr::select(-.data$SYMBOL)
        dat_exprs <- base::rbind(dat_exprs_with_duplicates,
                                 dat_exprs_no_duplicates)
    } else {
        base::row.names(dat_exprs) <- new_gene_names
    }
    return(dat_exprs)
}
