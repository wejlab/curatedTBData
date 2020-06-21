#' @keywords internal
#'
#' @importFrom stringr str_length
#' @importFrom stringr str_trunc
#' @importFrom magrittr %>%
#' @importFrom stringr str_pad
#' @importFrom stringr str_c
#' @importFrom utils download.file
#' @importFrom stringr str_subset
#' @importFrom stringr str_remove_all
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace_all
#' @importFrom tibble column_to_rownames
#' @importFrom tibble as_tibble
#' @importFrom tibble add_column
#' @importFrom dplyr select
#' @importFrom dplyr starts_with
#' @importFrom dplyr everything
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom stats sd
get_metadata <- function(series_accession) {
    accession_nchar <-
        stringr::str_length(series_accession)

    accession_short <-
        stringr::str_trunc(series_accession, 5, ellipsis = "") %>%
        stringr::str_pad(accession_nchar, side = "right", pad = "n")

    series_matrix_url <-
        stringr::str_c("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                       accession_short, "/", series_accession, "/matrix/",
                       series_accession, "_series_matrix.txt.gz")

    series_matrix_file <-
        base::tempfile()

    utils::download.file(series_matrix_url, series_matrix_file)

    series_matrix_connection <-
        base::gzfile(series_matrix_file)

    series_matrix_character <-
        base::readLines(series_matrix_connection)

    base::close(series_matrix_connection)

    series_metadata <-
        stringr::str_subset(series_matrix_character, "!Series_") %>%
        stringr::str_remove_all("!Series_") %>%
        readr::read_tsv(col_names = FALSE) %>%
        tidyr::pivot_wider(names_from = X1, names_prefix = "series_",
                           values_from = X2,
                           values_fn = base::list(X2 = base::list))

    sample_metadata <-
        stringr::str_subset(series_matrix_character, "!Sample_") %>%
        stringr::str_remove_all("!Sample_") %>%
        readr::read_tsv(col_names = FALSE) %>%
        dplyr::mutate(X1 = base::make.names(X1, unique = TRUE)) %>%
        dplyr::mutate(X1 = stringr::str_replace_all(X1, "\\.", "_")) %>%
        dplyr::mutate(X1 = stringr::str_c("sample_", X1)) %>%
        tibble::column_to_rownames(var = "X1") %>%
        base::t() %>%
        tibble::as_tibble()

    tibble::add_column(sample_metadata, series_metadata) %>%
        dplyr::select(dplyr::starts_with("series_"), dplyr::everything())
}
