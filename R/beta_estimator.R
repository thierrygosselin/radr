#' Estimate population-specific allelic FST
#'
#' @description
#' Estimates the population-specific allelic FST of Weir and Goudet (2017),
#' \eqn{\beta^i_{WT}}. Each value measures the differentiation of one population
#' relative to the complete set of populations included in the analysis.
#'
#' Use a quality-controlled, filtered dataset. This function does not perform
#' quality filtering, LD pruning, or sample removal. It requires diploid,
#' biallelic alternate-allele dosage (`ALT_DOSAGE`: 0, 1, 2, or `NA`). For GDS
#' input, biallelic status is checked before genotype chunks are read.
#'
#' The estimator uses only markers with at least one observed genotype in every
#' included population. This common-marker requirement does not modify the GDS
#' file. It ensures that all reported population estimates use the
#' same loci and the same population reference.
#'
#' @param data A `.gds` filepath or an open `SeqVarGDSClass` object.
#' Import other formats first with [genometranslator::read_genome()]. Active
#' sample and marker filters are respected and restored. Files opened here are
#' closed on exit; supplied open GDS objects remain open.
#' @param strata Optional strata filepath or data frame used to add or replace
#' population assignments. It must contain `INDIVIDUALS` and `STRATA`.
#' Default: \code{strata = NULL}.
#' @param populations Optional character vector of populations to include. Use
#' exactly two names for a separate two-population analysis. Its estimates have
#' a pair-specific reference and should not be compared directly with estimates
#' from another pair or from the full population set.
#' Default: \code{populations = NULL}.
#' @param random.mating Logical. Random mating means that, within each `STRATA`,
#' gametes unite independently of the genotypes carried by individuals. Under
#' this assumption, genotype frequencies are expected to follow Hardy-Weinberg
#' proportions given the observed allele frequencies. If `FALSE`, use
#' dosage-based allelic matching among distinct individuals without making this
#' assumption. If `TRUE`, use the allele-frequency and gene-diversity
#' formulation implemented by `hierfstat::betas()`.
#' Default: \code{random.mating = FALSE}.
#' @param filename Optional output prefix. When supplied, three tab-delimited
#' files are written with suffixes `_beta.tsv`, `_within_population.tsv`, and
#' `_between_populations.tsv`. Default: \code{filename = NULL}.
#' @param chunk.size Number of markers processed at once from the GDS.
#' Smaller chunks reduce peak memory use but may increase computation time.
#' Default: \code{chunk.size = 10000L}.
#' @param verbose Logical. Display progress and a population-specific FST
#' summary. Default: \code{verbose = TRUE}.
#'
#' @details
#' With `random.mating = FALSE`, allelic matching for individuals \eqn{j} and
#' \eqn{k} is calculated from
#' their alternate-allele dosages across jointly observed markers. At a marker,
#' matching is 1 for the same homozygote, 0 for opposite homozygotes, and 0.5
#' when at least one individual is heterozygous. Missing genotypes are omitted
#' separately for each pair of individuals.
#'
#' Let \eqn{M_W^i} be the mean matching between distinct individuals within
#' population \eqn{i}, and \eqn{M_B} the equally weighted mean matching across
#' every distinct pair of populations. Population-specific allelic FST is
#' \deqn{\beta^i_{WT} = \frac{M_W^i-M_B}{1-M_B}.}
#' This dosage-based matching formulation follows `hierfstat::fs.dosage()` and
#' does not require random mating within populations.
#'
#' With `random.mating = TRUE`, within-population gene diversity is estimated
#' from sample allele frequencies with the finite-sample correction used by
#' `hierfstat::betas()`. Population-specific allelic FST is then
#' \deqn{\beta^i_{WT} = 1 - \frac{\sum_l H_{W,li}}{\sum_l H_{B,l}}.}
#' Here, random mating is an assumption about mating and the union of gametes
#' within each `STRATA`; it is not simply a description of how samples were
#' collected. The assumption can be violated by inbreeding, self-fertilization,
#' assortative mating, close-family sampling, recent admixture, or hidden
#' population subdivision. These processes can cause observed genotype
#' frequencies to differ from the Hardy-Weinberg proportions predicted by
#' allele frequencies. The allele-frequency and matching formulations therefore
#' answer closely related, but not identical, questions and can produce
#' different estimates.
#'
#' The reference is defined by the sampled population set. Adding, removing, or
#' pairing populations changes that reference. Negative estimates are valid and
#' indicate less within-population allele matching than the sampled reference.
#'
#' Filtering should address genotyping quality, problematic samples, batch
#' effects, duplicated or structurally unsuitable markers, and the intended
#' dependence among loci. However, Weir and Goudet (2017) recommend including
#' rare variants for this estimator; filtering loci solely by allele frequency
#' can alter the estimate.
#'
#' GDS genotypes are processed by marker chunks and are not expanded into one
#' complete individual-by-marker tidy table. With `random.mating = FALSE`, exact
#' pairwise matching requires two square sample-by-sample accumulator matrices,
#' so memory use grows with the square of the number of individuals. The
#' allele-frequency formulation is lighter because it accumulates population
#' summaries instead.
#'
#' @return A named list containing:
#' \describe{
#'   \item{beta}{One row per population with `BETA`, the shared number of
#'   contributing markers, estimator name, and estimator-specific statistics.}
#'   \item{within_population}{Within-population matching summaries when
#'   `random.mating = FALSE`, or summed within-population gene diversity when
#'   `random.mating = TRUE`.}
#'   \item{between_populations}{Between-population matching summaries when
#'   `random.mating = FALSE`, or summed between-population gene diversity when
#'   `random.mating = TRUE`.}
#' }
#'
#' @references
#' Weir, B. S. and Goudet, J. (2017). A unified characterization of population
#' structure and relatedness. *Genetics*, 206, 2085-2103.
#' \doi{10.1534/genetics.116.198424}
#'
#' Goudet, J., Kay, T. and Weir, B. S. (2018). How to estimate kinship.
#' *Molecular Ecology*, 27, 4121-4135. \doi{10.1111/mec.14833}
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
#'
#' beta_random_mating <- beta_estimator(
#'   genome,
#'   random.mating = TRUE
#' )
#'
#' # A separate analysis whose reference contains only NORTH and SOUTH.
#' north_south <- beta_estimator(
#'   genome,
#'   populations = c("NORTH", "SOUTH")
#' )
#' }
#'
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
beta_estimator <- function(
    data,
    strata = NULL,
    populations = NULL,
    random.mating = FALSE,
    filename = NULL,
    chunk.size = 10000L,
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
  if (!is.logical(random.mating) || length(random.mating) != 1L ||
      is.na(random.mating)) {
    rlang::abort("`random.mating` must be TRUE or FALSE.")
  }
  if (!is.numeric(chunk.size) || length(chunk.size) != 1L ||
      !is.finite(chunk.size) || chunk.size < 1 ||
      chunk.size > .Machine$integer.max || chunk.size != floor(chunk.size)) {
    rlang::abort("`chunk.size` must be one positive integer.")
  }
  chunk.size <- as.integer(chunk.size)
  if (!is.null(filename) &&
      (!is.character(filename) || length(filename) != 1L ||
       is.na(filename) || !nzchar(filename))) {
    rlang::abort("`filename` must be NULL or one non-empty character value.")
  }
  if (!is.null(populations) &&
      (!is.character(populations) || anyNA(populations) ||
       any(!nzchar(populations)))) {
    rlang::abort("`populations` must be NULL or a character vector of population names.")
  }
  populations <- unique(populations)

  common_in_chunk <- function(dosage, sample.strata, population.levels) {
    if (!is.matrix(dosage)) dosage <- as.matrix(dosage)
    if (ncol(dosage) != length(sample.strata)) {
      rlang::abort("Internal error: genotype columns and sample strata are not aligned.")
    }
    called <- do.call(cbind, lapply(population.levels, function(population) {
      idx <- which(sample.strata == population)
      rowSums(!is.na(dosage[, idx, drop = FALSE])) > 0L
    }))
    rowSums(called) == length(population.levels)
  }

  update_statistics <- function(dosage, sample.strata, population.levels) {
    if (nrow(dosage) == 0L) return(invisible(NULL))
    if (!random.mating) {
      observed <- !is.na(dosage)
      centred <- dosage - 1
      centred[!observed] <- 0
      matching.numerator <<- matching.numerator + crossprod(centred)
      matching.denominator <<- matching.denominator + crossprod(observed * 1)
    } else {
      within.chunk <- vapply(population.levels, function(population) {
        x <- dosage[, sample.strata == population, drop = FALSE]
        n.chromosomes <- 2L * rowSums(!is.na(x))
        alt.frequency <- rowSums(x, na.rm = TRUE) / n.chromosomes
        ref.frequency <- 1 - alt.frequency
        (n.chromosomes / (n.chromosomes - 1L)) *
          (1 - alt.frequency^2 - ref.frequency^2)
      }, numeric(nrow(dosage)))
      if (nrow(dosage) == 1L) {
        within.chunk <- matrix(within.chunk, nrow = 1L)
      }
      alt.frequency <- vapply(population.levels, function(population) {
        x <- dosage[, sample.strata == population, drop = FALSE]
        rowSums(x, na.rm = TRUE) / (2L * rowSums(!is.na(x)))
      }, numeric(nrow(dosage)))
      if (nrow(dosage) == 1L) {
        alt.frequency <- matrix(alt.frequency, nrow = 1L)
      }
      ref.frequency <- 1 - alt.frequency
      n.populations <- length(population.levels)
      alt.similarity <- (rowSums(alt.frequency)^2 -
        rowSums(alt.frequency^2)) / 2
      ref.similarity <- (rowSums(ref.frequency)^2 -
        rowSums(ref.frequency^2)) / 2
      between.chunk <- 1 -
        (2 / (n.populations * (n.populations - 1))) *
        (alt.similarity + ref.similarity)
      sum.within <<- sum.within + colSums(within.chunk)
      sum.between <<- sum.between + sum(between.chunk)
    }
    n.common.markers <<- n.common.markers + nrow(dosage)
    invisible(NULL)
  }

  close.gds <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L &&
             !is.na(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    if (!file.exists(data)) rlang::abort("The GDS file does not exist.")
    gds <- SeqArray::seqOpen(data, readonly = TRUE)
    close.gds <- TRUE
  } else {
    rlang::abort(
      "`data` must be a .gds filepath or an open SeqVarGDSClass object."
    )
  }
  on.exit({
    if (close.gds) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)
  SeqArray::seqSetFilter(gds, action = "push", verbose = FALSE)
  on.exit(try(SeqArray::seqSetFilter(
    gds, action = "pop", verbose = FALSE
  ), silent = TRUE), add = TRUE, after = FALSE)

  if (!isTRUE(genometranslator::detect_biallelic_markers(gds, verbose = FALSE))) {
    rlang::abort("`beta_estimator()` requires biallelic markers.")
  }

  active.samples <- SeqArray::seqGetData(gds, "sample.id")
  individuals <- if (is.null(strata)) {
    genometranslator::extract_individuals_metadata(
      gds = gds, ind.field.select = c("INDIVIDUALS", "STRATA"),
      whitelist = TRUE
    )
  } else {
    genometranslator::join_strata(
      data = tibble::tibble(INDIVIDUALS = active.samples),
      strata = strata, verbose = FALSE
    )
  }
  if (!all(c("INDIVIDUALS", "STRATA") %in% names(individuals))) {
    rlang::abort("Population assignments are missing. Supply GDS strata metadata or `strata`.")
  }
  individuals <- dplyr::select(individuals, INDIVIDUALS, STRATA)
  individuals <- individuals[individuals$INDIVIDUALS %in% active.samples, ]
  if (anyDuplicated(individuals$INDIVIDUALS) ||
      anyNA(individuals$STRATA) || any(!nzchar(individuals$STRATA))) {
    rlang::abort("Each individual must have one non-missing STRATA assignment.")
  }
  if (!is.null(populations)) {
    unknown <- setdiff(populations, unique(as.character(individuals$STRATA)))
    if (length(unknown) > 0L) {
      rlang::abort(paste0("Unknown populations: ", paste(unknown, collapse = ", "), "."))
    }
    individuals <- dplyr::filter(individuals, STRATA %in% populations)
  }
  population.levels <- if (is.null(populations)) {
    unique(as.character(individuals$STRATA))
  } else {
    populations
  }
  if (length(population.levels) < 2L) {
    rlang::abort("At least two populations are required.")
  }

  marker.ids <- SeqArray::seqGetData(gds, "variant.id")
  marker.chunks <- split(
    marker.ids,
    ceiling(seq_along(marker.ids) / chunk.size)
  )
  if (!random.mating) {
    matching.numerator <- matrix(0, nrow(individuals), nrow(individuals))
    matching.denominator <- matrix(0, nrow(individuals), nrow(individuals))
  } else {
    sum.within <- stats::setNames(numeric(length(population.levels)), population.levels)
    sum.between <- 0
  }
  n.common.markers <- 0L
  analysis.individuals <- NULL
  analysis.strata <- NULL
  for (i in seq_along(marker.chunks)) {
    variant.chunk <- marker.chunks[[i]]
    SeqArray::seqSetFilter(
      gds,
      sample.id = individuals$INDIVIDUALS,
      variant.id = variant.chunk,
      action = "set",
      verbose = FALSE
    )
    sample.order <- SeqArray::seqGetData(gds, "sample.id")
    dosage.list <- SeqArray::seqApply(
      gds,
      "$dosage_alt",
      FUN = function(x) as.numeric(x),
      margin = "by.variant",
      as.is = "list"
    )
    dosage <- do.call(rbind, dosage.list)
    sample.strata <- individuals$STRATA[
      match(sample.order, individuals$INDIVIDUALS)
    ]
    if (is.null(analysis.individuals)) {
      analysis.individuals <- sample.order
      analysis.strata <- as.character(sample.strata)
    }
    keep <- common_in_chunk(dosage, sample.strata, population.levels)
    update_statistics(
      dosage = dosage[keep, , drop = FALSE],
      sample.strata = sample.strata,
      population.levels = population.levels
    )
    if (verbose && length(marker.chunks) > 1L) {
      message("  Marker chunk ", i, " of ", length(marker.chunks), " completed")
    }
  }
  if (n.common.markers == 0L) {
    rlang::abort("No marker has observed genotypes in every included population.")
  }
  if (!random.mating) {
    if (verbose) message("Estimator: allelic matching without a random-mating assumption")
    matching <- 0.5 * (1 + matching.numerator / matching.denominator)
    diag(matching) <- NA_real_
    dimnames(matching) <- list(analysis.individuals, analysis.individuals)

    within.population <- dplyr::bind_rows(lapply(population.levels, function(population) {
      idx <- which(analysis.strata == population)
      block <- matching[idx, idx, drop = FALSE]
      tibble::tibble(
        STRATA = population,
        N_INDIVIDUALS = length(idx),
        N_INFORMATIVE_PAIRS = sum(is.finite(block)) / 2L,
        MATCHING_WITHIN = mean(block, na.rm = TRUE)
      )
    }))
    if (any(!is.finite(within.population$MATCHING_WITHIN))) {
      rlang::abort("Every population requires at least two individuals with jointly observed markers.")
    }

    population.pairs <- utils::combn(population.levels, 2L, simplify = FALSE)
    between.populations <- dplyr::bind_rows(lapply(population.pairs, function(pair) {
      idx.1 <- which(analysis.strata == pair[[1]])
      idx.2 <- which(analysis.strata == pair[[2]])
      block <- matching[idx.1, idx.2, drop = FALSE]
      tibble::tibble(
        STRATA_1 = pair[[1]],
        STRATA_2 = pair[[2]],
        N_INFORMATIVE_PAIRS = sum(is.finite(block)),
        MATCHING_BETWEEN = mean(block, na.rm = TRUE)
      )
    }))
    reference.matching <- mean(between.populations$MATCHING_BETWEEN)
    if (!is.finite(reference.matching) || reference.matching >= 1) {
      rlang::abort("Between-population matching does not provide a finite reference below 1.")
    }
    beta <- within.population |>
      dplyr::mutate(
        ESTIMATOR = "ALLELIC_MATCHING",
        N_MARKERS = n.common.markers,
        MATCHING_REFERENCE = reference.matching,
        BETA = (MATCHING_WITHIN - MATCHING_REFERENCE) /
          (1 - MATCHING_REFERENCE)
      ) |>
      dplyr::select(
        STRATA, ESTIMATOR, N_INDIVIDUALS, N_INFORMATIVE_PAIRS, N_MARKERS,
        MATCHING_WITHIN, MATCHING_REFERENCE, BETA
      )
  } else {
    if (verbose) message("Estimator: allele frequencies with random mating")
    within.population <- tibble::tibble(
      STRATA = population.levels,
      N_MARKERS = n.common.markers,
      SUM_HW = unname(sum.within)
    )
    between.populations <- tibble::tibble(
      N_POPULATIONS = length(population.levels),
      N_MARKERS = n.common.markers,
      SUM_HB = sum.between
    )
    beta <- within.population |>
      dplyr::mutate(
        ESTIMATOR = "ALLELE_FREQUENCY_RANDOM_MATING",
        SUM_HB = sum.between,
        BETA = 1 - (SUM_HW / SUM_HB)
      ) |>
      dplyr::select(STRATA, ESTIMATOR, N_MARKERS, SUM_HW, SUM_HB, BETA)
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
    message(
      "Markers shared by all ", length(population.levels),
      " populations: ", n.common.markers
    )
    message("Population-specific allelic FST:")
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
