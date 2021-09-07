#' Merge samples with common genes from selected studies
#' @name combineObjects
#' @param object_list A \code{list} of \link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}
#' or \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} objects
#' The assays contains object's assay contain expression data with probes mapped to gene symbol
#' \code{names(object_list)} should not be \code{NULL}
#' @param experiment_name A character/vector of character to choose the name of the assay from the input
#' \code{list} of object
#' @return A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object contains combined data from the input
#' @examples
#' geo <-  c("GSE19435", "GSE19439")
#' data_list <-  curatedTBData(geo, dryrun = FALSE, curated.only = TRUE)
#' combineObjects(data_list, experiment_name = "assay_curated")
#' @export
combineObjects <- function(object_list, experiment_name) {
    ## check the experiment_name argument
    if (base::missing(experiment_name)) {
        base::stop("Argument \"experiment_name\" is missing, with no default.")
    }
    ## check length of the list, should be greater than 1
    if (base::length(object_list) <= 1L) {
        base::stop(sprintf("The length of the input list is %i, expecting more than 1 elments within list for combining objects.",
                           base::length(object_list)))
    }
    ## check names of the input object list
    obj_name <- base::names(object_list)
    if (base::is.null(obj_name)) {
        base::stop("names of the input list should not be NULL. Add unique name for each object within the list.")
    } else if (!base::is.na(base::match("", obj_name))) {
        base::stop("names of the input contains \"\". Replace \"\" with unique character.")
    }
    ## Check whether it is a list of SummarizedExperiment or MultiAssayExperiment objects
    isSummarizedExperiment <- base::all(base::vapply(object_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE))
    isMultiAssayExperiement <- base::all(base::vapply(object_list, function(x)
        methods::is(x, "MultiAssayExperiment"), TRUE))
    if (isSummarizedExperiment) {
        dat_exprs_match <- .select_assay(object_list, experiment_name, Sobject = TRUE)
    } else if (isMultiAssayExperiement) {
        dat_exprs_match <- .select_assay(object_list, experiment_name, Sobject = FALSE)
    } else {
        base::stop("Input is not a list of MultiAssayExperiment or SummarizedExperiment objetcs.")
    }
    ## Combine sample with common genes from a list of objects. Input data type should be data.frame
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
    col_data <- base::lapply(base::seq_len(base::length(object_list)), function(x) {
        col_data <- SummarizedExperiment::colData(object_list[[x]])
        col_data$Study <- base::names(object_list[x])
        base::as.data.frame(col_data)
    })
    ## Combine list into data frame with unequal columns, fill in NA when columns from studies are not found
    rbindx <- function(dfs) {
        ns <- base::unique(base::unlist(base::lapply(dfs, base::colnames)))
        base::do.call(base::rbind, base::lapply(dfs, function(x) {
            for (n in ns[!ns %in% base::colnames(x)]) {
                x[[n]] <- NA
            }
            x
        }))
    }
    col_info <- rbindx(col_data)
    Sample <- base::vapply(object_list, function(x)
        SummarizedExperiment::colData(x) %>%
            base::row.names(), base::character(1))
    base::row.names(col_info) <- Sample
    ## Remove samples that does not exist in the count
    index <- stats::na.omit(base::match(base::colnames(dat_exprs_count), Sample))
    col_info <- col_info[index, ]
    ## Create output in the format of SummarizedExperiment
    result <- SummarizedExperiment::SummarizedExperiment(assays = base::list(assay1 = base::as.matrix(dat_exprs_count)),
                                                         colData = col_info)
    return(result)
}

#' Select assay based on input list type
#'
#' @inheritParams combineObjects
#' @param Sobject Boolean. Indicate whether the input is a \code{list} of
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} objects.
#'
.select_assay <- function(object_list, experiment_name, Sobject) {
    ## Merge code starts here
    if (base::length(experiment_name) > 1L) {
        message("Found more than one \"experiment_name\".")
        ## experiment name for the list of object is different
        if (base::length(experiment_name) == base::length(object_list)) {
            dat_exprs_match <- base::mapply(function(i, y) {
                x <- object_list[[i]]
                if (Sobject) {
                    dat_assay <- SummarizedExperiment::assays(x)[[y]]
                } else {
                    dat_assay <- MultiAssayExperiment::experiments(x)[[y]]
                }
                if (base::is.null(dat_assay)) {
                    base::stop(base::sprintf("Object: %s with experiment name: %s has assay NULL.",
                                             base::names(object_list)[i], experiment_name))
                }
                base::as.data.frame(dat_assay)
            }, base::seq_len(base::length(object_list)), experiment_name)
            names(dat_exprs_match) <- names(object_list)
        } else {
            base::stop("The length of input list is not the same as the length of the \"experiment_name\" vector.")
        }
    } else {
        dat_exprs_match <- base::lapply(base::seq_len(base::length(object_list)), function(i) {
            x <- object_list[[i]]
            if (Sobject) {
                dat_assay <- SummarizedExperiment::assays(x)[[experiment_name]]
            } else {
                dat_assay <- MultiAssayExperiment::experiments(x)[[experiment_name]]
            }
            if (base::is.null(dat_assay)) {
                base::stop(base::sprintf("Object: %s with experiment name: %s has assay NULL.",
                                         base::names(object_list)[i], experiment_name))
            }
            base::as.data.frame(dat_assay)
        })
        base::names(dat_exprs_match) <- base::names(object_list)
    }
    return(dat_exprs_match)
}
