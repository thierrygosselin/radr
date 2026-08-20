#' Estimate population-specific beta
#'
#' @description
#' Estimates population-specific \eqn{\beta_i} following Weir and Goudet
#' (2017). The implementation uses diploid, biallelic alternate-allele dosage
#' (`ALT_DOSAGE`: 0, 1, 2, or `NA`) and works directly from a GDS file, an open
#' `SeqVarGDSClass` object, or a tidy genotype table.
#'
#' `beta_estimator()` does not filter individuals or markers. Missing genotypes
#' are omitted within each marker-population combination. A marker contributes
#' only when it has observed genotypes in at least two populations and its
#' between-population diversity is finite.
#'
#' @param data A GDS filepath, an open `SeqVarGDSClass` object, or a tidy data
#' frame containing `MARKERS`, `INDIVIDUALS`, `STRATA` (or `POP_ID`), and
#' `ALT_DOSAGE`.
#' @param strata Optional strata filepath or data frame used to add or replace
#' population assignments. It must contain `INDIVIDUALS` and `STRATA`.
#' Default: \code{strata = NULL}.
#' @param filename Optional output prefix. When supplied, three tab-delimited
#' files are written with suffixes `_beta.tsv`, `_within_population.tsv`, and
#' `_between_populations.tsv`. Default: \code{filename = NULL}.
#' @param parallel.core Number of processor cores passed to
#' [genometranslator::read_genome()] when file input must be imported.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param verbose Logical. Display progress and a population beta summary.
#' Default: \code{verbose = TRUE}.
#'
#' @details
#' For marker \eqn{l} and population \eqn{i}, within-population gene diversity
#' is estimated as
#' \deqn{H_{W,li} = \frac{n_{li}}{n_{li}-1}
#'   \left(1 - p_{li}^2 - (1-p_{li})^2\right),}
#' where \eqn{n_{li}} is the number of observed gene copies and \eqn{p_{li}}
#' is the alternate-allele frequency.
#'
#' Between-population diversity is calculated for every marker from the
#' populations with observed genotypes. Population-specific beta is then
#' \deqn{\beta_i = 1 - \frac{\sum_l H_{W,li}}{\sum_l H_{B,l}}.}
#' Consequently, populations can be compared using different numbers of loci
#' when their missingness patterns differ. Review `N_MARKERS` in the returned
#' summary before interpreting differences among populations.
#'
#' @return A named list containing:
#' \describe{
#'   \item{beta}{One row per population with `BETA`, the number of contributing
#'   markers, and the summed within- and between-population diversities.}
#'   \item{within_population}{Marker-population allele counts, frequencies, and
#'   `HW`.}
#'   \item{between_populations}{Marker-level number of represented populations
#'   and `HB`.}
#' }
#'
#' @references
#' Weir, B. S. and Goudet, J. (2017). A unified characterization of population
#' structure and relatedness. *Genetics*, 206, 2085–2103.
#' \doi{10.1534/genetics.116.198424}
#'
#' Goudet, J., Kay, T. and Weir, B. S. (2018). How to estimate kinship.
#' *Molecular Ecology*, 27, 4121–4135. \doi{10.1111/mec.14833}
#'
#' @examples
#' \dontrun{
#' genome <- genometranslator::read_genome(
#'   data = "individuals.vcf.gz",
#'   strata = "strata.tsv"
#' )
#'
#' beta <- beta_estimator(genome)
#' beta$beta
#' beta$within_population
#' beta$between_populations
#' }
#'
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
beta_estimator <- function(
    data,
    strata = NULL,
    filename = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {
  .start <- tgbase::startup(
    package = "radr",
    f.name = "beta_estimator",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (missing(data)) rlang::abort("Argument `data` is required.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    rlang::abort("`verbose` must be TRUE or FALSE.")
  }
  if (!is.numeric(parallel.core) || length(parallel.core) != 1L ||
      is.na(parallel.core) || parallel.core < 1) {
    rlang::abort("`parallel.core` must be one positive number.")
  }
  parallel.core <- as.integer(parallel.core)
  if (!is.null(filename) &&
      (!is.character(filename) || length(filename) != 1L ||
       is.na(filename) || !nzchar(filename))) {
    rlang::abort("`filename` must be NULL or one non-empty character value.")
  }

  gds <- NULL
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.data.frame(data)) {
    genotype.data <- tibble::as_tibble(data)
  } else {
    if (verbose) message("Importing genomic data with genometranslator...")
    gds <- genometranslator::read_genome(
      data = data,
      strata = strata,
      parallel.core = parallel.core,
      verbose = FALSE
    )
    on.exit(try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  }

  if (!is.null(gds)) {
    genotype.data <- genometranslator::extract_genotypes_metadata(
      gds = gds,
      genotypes.meta.select = c("MARKERS", "INDIVIDUALS", "ALT_DOSAGE"),
      whitelist = TRUE
    )
    if (is.null(strata)) {
      strata <- genometranslator::extract_individuals_metadata(
        gds = gds,
        ind.field.select = c("INDIVIDUALS", "STRATA"),
        whitelist = TRUE
      )
    }
  }

  if (!is.null(strata)) {
    genotype.data <- genometranslator::join_strata(
      data = genotype.data,
      strata = strata,
      verbose = FALSE
    )
  }
  if ("POP_ID" %in% names(genotype.data) &&
      !"STRATA" %in% names(genotype.data)) {
    genotype.data <- dplyr::rename(genotype.data, STRATA = POP_ID)
  }

  required.columns <- c("MARKERS", "INDIVIDUALS", "STRATA", "ALT_DOSAGE")
  missing.columns <- setdiff(required.columns, names(genotype.data))
  if (length(missing.columns) > 0L) {
    rlang::abort(paste0(
      "The input is missing required columns: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }

  genotype.data <- dplyr::select(
    genotype.data, MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE
  )
  if (anyNA(genotype.data$MARKERS) || anyNA(genotype.data$INDIVIDUALS) ||
      anyNA(genotype.data$STRATA)) {
    rlang::abort("`MARKERS`, `INDIVIDUALS`, and `STRATA` cannot contain missing values.")
  }

  observed.dosage <- genotype.data$ALT_DOSAGE[!is.na(genotype.data$ALT_DOSAGE)]
  if (length(observed.dosage) == 0L) {
    rlang::abort("No observed `ALT_DOSAGE` values are available.")
  }
  if (!is.numeric(observed.dosage) ||
      any(!observed.dosage %in% c(0, 1, 2))) {
    rlang::abort(
      "`ALT_DOSAGE` must contain only diploid biallelic dosages 0, 1, 2, or NA."
    )
  }

  duplicated.genotypes <- dplyr::count(
    genotype.data, MARKERS, INDIVIDUALS, name = "N"
  ) |>
    dplyr::filter(N > 1L)
  if (nrow(duplicated.genotypes) > 0L) {
    rlang::abort(
      "The input contains more than one genotype for at least one marker-individual pair."
    )
  }

  if (verbose) message("Computing population-specific beta...")
  within.population <- genotype.data |>
    dplyr::filter(!is.na(ALT_DOSAGE)) |>
    dplyr::group_by(MARKERS, STRATA) |>
    dplyr::summarise(
      N_INDIVIDUALS = dplyr::n(),
      N_CHROMOSOMES = 2L * N_INDIVIDUALS,
      ALT_COUNT = sum(ALT_DOSAGE),
      ALT_FREQ = ALT_COUNT / N_CHROMOSOMES,
      REF_FREQ = 1 - ALT_FREQ,
      HW = dplyr::if_else(
        N_CHROMOSOMES > 1L,
        (N_CHROMOSOMES / (N_CHROMOSOMES - 1)) *
          (1 - ALT_FREQ^2 - REF_FREQ^2),
        NA_real_
      ),
      .groups = "drop"
    )

  between.populations <- within.population |>
    dplyr::group_by(MARKERS) |>
    dplyr::summarise(
      N_POPULATIONS = dplyr::n(),
      ALT_PAIR_SIMILARITY = (sum(ALT_FREQ)^2 - sum(ALT_FREQ^2)) / 2,
      REF_PAIR_SIMILARITY = (sum(REF_FREQ)^2 - sum(REF_FREQ^2)) / 2,
      HB = dplyr::if_else(
        N_POPULATIONS > 1L,
        1 - (2 / (N_POPULATIONS * (N_POPULATIONS - 1))) *
          (ALT_PAIR_SIMILARITY + REF_PAIR_SIMILARITY),
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(N_POPULATIONS >= 2L, is.finite(HB)) |>
    dplyr::select(MARKERS, N_POPULATIONS, HB)

  if (nrow(between.populations) == 0L) {
    rlang::abort("No marker has observed genotypes in at least two populations.")
  }

  beta <- within.population |>
    dplyr::inner_join(between.populations, by = "MARKERS") |>
    dplyr::filter(is.finite(HW), is.finite(HB)) |>
    dplyr::group_by(STRATA) |>
    dplyr::summarise(
      N_MARKERS = dplyr::n(),
      SUM_HW = sum(HW),
      SUM_HB = sum(HB),
      BETA = dplyr::if_else(SUM_HB > 0, 1 - (SUM_HW / SUM_HB), NA_real_),
      .groups = "drop"
    )
  if (nrow(beta) < 2L) {
    rlang::abort("Beta estimation requires at least two represented populations.")
  }

  result <- list(
    beta = beta,
    within_population = within.population,
    between_populations = between.populations
  )

  if (!is.null(filename)) {
    extension <- tools::file_ext(filename)
    prefix <- if (tolower(extension) == "tsv") {
      substr(filename, 1L, nchar(filename) - 4L)
    } else {
      filename
    }
    output.files <- c(
      beta = paste0(prefix, "_beta.tsv"),
      within_population = paste0(prefix, "_within_population.tsv"),
      between_populations = paste0(prefix, "_between_populations.tsv")
    )
    readr::write_tsv(result$beta, output.files[["beta"]])
    readr::write_tsv(result$within_population, output.files[["within_population"]])
    readr::write_tsv(result$between_populations, output.files[["between_populations"]])
    result$files <- output.files
    if (verbose) {
      message("Files written: ", paste(basename(output.files), collapse = ", "))
    }
  }

  if (verbose) {
    message("Population-specific beta:")
    apply(beta, 1L, function(x) {
      message(
        "  ", x[["STRATA"]], " = ",
        format(round(as.numeric(x[["BETA"]]), 4), nsmall = 4),
        " (", x[["N_MARKERS"]], " markers)"
      )
    })
  }
  result
}


#' Legacy plural name for population-specific beta
#'
#' `betas_estimator()` is retained temporarily for compatibility. New code
#' should use [beta_estimator()].
#'
#' @inheritParams beta_estimator
#' @return The result returned by [beta_estimator()].
#' @export
betas_estimator <- function(
    data,
    strata = NULL,
    filename = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {
  .Deprecated("beta_estimator", package = "radr")
  beta_estimator(
    data = data,
    strata = strata,
    filename = filename,
    parallel.core = parallel.core,
    verbose = verbose
  )
}
