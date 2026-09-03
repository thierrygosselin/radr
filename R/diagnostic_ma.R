# Minor Allele Diagnostic
#' @title MA diagnostic
#' @description Minor Allele diagnostic, help choose a filter threshold.
#' @param data A file in the working directory or object in the global environment
#' in wide or long (tidy) formats. To import, the function uses
#' \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' \code{\link[genometranslator]{read_genome}}.
#' \emph{See details of this function for more info}.
#' @param group.rank (Number) The number of group to class the MAF.
#' @param filename (optional) Name of the file written to the working directory.
#' @param strata Optional sample metadata containing `INDIVIDUALS` and
#'   `group.column`.
#' @param group.column Metadata column used to calculate local MAF.
#' @rdname diagnostic_maf
#' @keywords internal
#' @export
#' @details Highly recommended to look at the distribution of MAF

#' @examples
#' \dontrun{
#' problem <- radr::diagnostic_ma(
#' data = tidy.salmon.data, group.rank = 10)
#' }


#' @seealso \link{filter_ma}


diagnostic_ma <- function(
    data, group.rank, filename = NULL, strata = NULL,
    group.column = "STRATA"
){

  if (missing(data)) rlang::abort("Input file missing")
  if (missing(group.rank)) rlang::abort("group.rank argument missing")

  data <- summarise_genomic_data(
    data = data,
    strata = strata,
    group.column = group.column,
    by.strata = TRUE,
    write.files = FALSE,
    verbose = FALSE,
    internal = TRUE,
    path.folder = tempdir()
  )$stratum.statistics

  # Local
  test.local <- data %>%
    dplyr::filter(.data$GROUP != "OVERALL") %>%
    dplyr::group_by(.data$VARIANT_ID) %>%
    dplyr::summarise(
      MAF_L = mean(.data$MINOR_ALLELE_FREQUENCY, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(RANK = dplyr::ntile(MAF_L, group.rank)) %>%
    dplyr::summarise(
      LOCAL_MAF = mean(MAF_L, na.rm = TRUE),
      n = dplyr::n()
    ) %>%
    dplyr::select(-n)

  # Global

  test.global <- data %>%
    dplyr::filter(.data$GROUP == "OVERALL") %>%
    dplyr::transmute(
      VARIANT_ID = .data$VARIANT_ID,
      GLOBAL_MAF = .data$MINOR_ALLELE_FREQUENCY
    ) %>%
    dplyr::group_by(RANK = dplyr::ntile(GLOBAL_MAF, group.rank)) %>%
    dplyr::summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = TRUE),
      n = dplyr::n()
    ) %>%
    dplyr::select(GLOBAL_MAF, n)

  maf.diagnostic <- dplyr::bind_cols(test.local, test.global)

  if (is.null(filename)) {
    message("Saving file: not selected")
  } else {
    message("Saving file: selected")
    readr::write_tsv(maf.diagnostic, filename, append = FALSE, col_names = TRUE)
    message("filename in the working directory: ", filename)
  }
  return(maf.diagnostic)
}
