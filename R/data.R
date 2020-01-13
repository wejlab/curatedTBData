#' A list of published TB signatures.
#'
#' A set of 34 Tuberculosis gene signatures from various publications. This set
#' of signatures uses gene symbols. Attempts have been made to use updated gene
#' symbols and remove symbols that did not match the most recent annotation.
#' Additional sets for Entrez IDs and Ensembl IDs are forthcoming.
#'
#' Signature names are composed of the last name of the primary author,
#' followed by
#' a possible context for the signature, and ending with either the number of
#' gene transcripts or genes in the signature, with respect to however
#' it was described in the signature in the original publication.
#'
#' Possible signature contexts:
#' \itemize{
#' \item{OD: Other diseases}
#' \item{HIV: Human Immunodeficiency Virus}
#' \item{PNA: Pneumonia}
#' \item{RISK: Risk of developing active TB}
#' \item{RES: Response to TB treatment}
#' \item{FAIL: Failure of TB treatment}
#' }
#'
#' Note that in some cases signatures will be positive identifiers of TB
#' whereas others are negative identifiers; this should be taken into account
#' when creating ROC curves and computing any AUC estimates.
#'
