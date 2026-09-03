#' Diagnose biallelic and multiallelic variants
#'
#' Inspect active GDS variants and report the number of alleles declared by
#' each record and observed among active samples. The function is diagnostic
#' and does not modify GDS filters.
#'
#' `N_DECLARED_ALLELES` comes from the GDS `allele` node.
#' `N_OBSERVED_ALLELES` is calculated from called genotypes. A declared
#' multiallelic record can therefore have only one or two alleles represented
#' in the active sample.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param chunk.size Number of active variants read at a time.
#' @param write.files Write diagnostic tables.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return A `detect_biallelic_problems` object containing `statistics`,
#' `multiallelic`, `biallelic`, `summary`, and output information.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' variants <- radr::detect_biallelic_problems("study.gds")
#' variants$multiallelic
#' }
detect_biallelic_problems <- function(
    data, chunk.size = 2000L, write.files = TRUE, verbose = TRUE, ...
) {
  context <- .analysis_gds_start(
    data, "detect_biallelic_problems", write.files, verbose, ...
  )
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  alleles <- as.character(SeqArray::seqGetData(gds, "allele"))
  metadata <- .sex_marker_metadata(gds, variant.id)
  chunks <- split(seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size))
  observed <- unlist(lapply(chunks, function(index) {
    .biallelic_observed_counts(gds, variant.id[index])
  }), use.names = FALSE)
  declared <- vapply(strsplit(alleles, ",", fixed = TRUE), length, integer(1))
  statistics <- metadata |>
    dplyr::mutate(
      ALLELES = alleles,
      N_DECLARED_ALLELES = declared,
      N_ALT_ALLELES = pmax(declared - 1L, 0L),
      N_OBSERVED_ALLELES = observed,
      RECORD_TYPE = dplyr::case_when(
        declared < 2L ~ "INVALID_OR_MONOMORPHIC_RECORD",
        declared == 2L ~ "BIALLELIC",
        TRUE ~ "MULTIALLELIC"
      ),
      OBSERVED_TYPE = dplyr::case_when(
        observed == 0L ~ "ALL_MISSING",
        observed == 1L ~ "MONOMORPHIC_IN_ACTIVE_SAMPLES",
        observed == 2L ~ "BIALLELIC_IN_ACTIVE_SAMPLES",
        TRUE ~ "MULTIALLELIC_IN_ACTIVE_SAMPLES"
      )
    )
  multiallelic <- dplyr::filter(statistics, .data$N_DECLARED_ALLELES > 2L)
  biallelic <- dplyr::filter(statistics, .data$N_DECLARED_ALLELES == 2L)
  summary <- statistics |>
    dplyr::count(.data$RECORD_TYPE, .data$OBSERVED_TYPE, name = "N_VARIANTS")
  output.files <- character()
  if (write.files) {
    tables <- list(statistics, multiallelic, summary)
    output.files <- file.path(context$path.folder, c(
      "biallelic_variant_statistics.tsv",
      "multiallelic_variants.tsv",
      "biallelic_problem_summary.tsv"
    ))
    for (i in seq_along(output.files)) {
      readr::write_tsv(tables[[i]], output.files[[i]], na = "NA")
    }
  }
  if (verbose) {
    .summary_message(
      "Declared multiallelic variants: ", nrow(multiallelic), " of ",
      nrow(statistics), "."
    )
  }
  restored <- .analysis_gds_restore(context)
  out <- list(
    statistics = statistics, multiallelic = multiallelic,
    biallelic = biallelic, summary = summary,
    path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("detect_biallelic_problems", class(out))
  out
}

.biallelic_observed_counts <- function(gds, variant.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  genotype <- SeqArray::seqGetData(gds, "genotype")
  if (length(dim(genotype)) != 3L) {
    rlang::abort("The GDS genotype node is not a diploid variant array.")
  }
  vapply(seq_along(variant.id), function(i) {
    alleles <- genotype[, , i]
    length(unique(as.integer(alleles[alleles >= 0L])))
  }, integer(1))
}
