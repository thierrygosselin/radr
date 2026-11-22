#' Explore and filter genomic data
#'
#' @description
#' Run a guided or scripted quality-control workflow on a Genomic Data
#' Structure (GDS) file. The workflow applies filters in a defined order,
#' records exclusions in the GDS, and writes diagnostics for the resulting
#' sample and marker panels. Choose filters and thresholds based on the
#' sequencing design, marker type, and analysis objective.
#'
#' @param data A GDS filepath or an open, writable \code{SeqVarGDSClass}
#'   connection. Use [genometranslator::read_vcf()], [genometranslator::read_dart()]
#'   or [genometranslator::read_plink()] to import other formats.
#' @param interactive.filter Logical. \code{TRUE} (default) offers diagnostic
#'   plots and prompts for thresholds. \code{FALSE} uses named settings
#'   supplied through \code{...} and the workflow defaults.
#' @param output Optional character vector of output formats passed to
#'   [genometranslator::genome_translator()]. The default, \code{NULL}, skips
#'   format conversion.
#' @param filename Optional filename or prefix for the translated output.
#'   In the current workflow the export uses the default basename
#'   \code{"filtered_genome"}. For a custom filename, export the returned
#'   GDS with [genometranslator::genome_translator()].
#' @param verbose Logical. Display progress messages. Default: \code{TRUE}.
#' @param parallel.core Number of processor cores passed to supported
#'   operations. Default: \code{parallel::detectCores() - 1}.
#' @param ... Named options for individual filters and workflow outputs. See
#'   the sections below for the available options and their batch defaults.
#'
#' @details
#' The workflow updates the supplied GDS in place. Use a writable working copy
#' for each comparison that starts from the same set of samples and markers.
#' At the end of the workflow, the active GDS selection is synchronized with
#' the saved sample and marker whitelists.
#'
#' @section Filter order and options in \code{...}:
#' The list gives the order of operations. Unless stated otherwise, optional
#' thresholds default to \code{NULL} and their filters are skipped in batch
#' mode when no threshold is supplied.
#' \enumerate{
#' \item \strong{Known sample exclusions.} Supply \code{blacklist.id} as a
#'   character vector of exact GDS sample IDs or a tab-separated file with one
#'   column named \code{INDIVIDUALS}. Use character IDs to preserve leading
#'   zeros. [filter_individuals()] applies the list without prompting,
#'   reports unknown IDs as errors before applying exclusions, and keeps
#'   earlier filtering reasons. When samples are newly excluded, it also
#'   checks for genotype-monomorphic markers before the next step.
#' \item \strong{Diversity Arrays Technology (DArT) reproducibility.} For
#'   DArT data with reproducibility scores, set \code{filter.reproducibility}
#'   to the minimum score or \code{"outliers"}. See
#'   [filter_dart_reproducibility()].
#' \item \strong{Marker panel.} The workflow checks for monomorphic markers
#'   and keeps markers present in every stratum by default
#'   (\code{filter.monomorphic = TRUE}, \code{filter.common.markers = TRUE}).
#'   Set either flag to \code{FALSE} to skip that dedicated step. Sample
#'   filters can still trigger a monomorphic-marker check. For an explicit
#'   marker panel, apply [filter_whitelist()] directly before exploration.
#'   The wrapper currently leaves \code{biallelic} undefined when a
#'   \code{whitelist.markers} list is supplied.
#' \item \strong{Sample quality.} \code{filter.individuals.missing} is the
#'   maximum proportion of missing genotypes per sample.
#'   \code{filter.individuals.heterozygosity} specifies lower and upper
#'   limits for observed heterozygosity. With depth data, use
#'   \code{filter.individuals.coverage.total},
#'   \code{filter.individuals.coverage.median}, or
#'   \code{filter.individuals.coverage.iqr} for per-sample depth filters.
#'   The last uses the interquartile range of depth. See
#'   [filter_individuals()] for numeric bounds and outlier rules.
#' \item \strong{Marker quality.} Set \code{filter.ma} to the minimum global
#'   minor allele count. Use \code{filter.coverage} for marker mean-depth
#'   bounds when depth data exist, and \code{filter.genotyping} for the
#'   maximum proportion of missing genotypes per marker. See
#'   [filter_ma()], [filter_coverage()], and [filter_genotyping()].
#' \item \strong{Marker position and number.}
#'   \code{filter.snp.position.read} screens single-nucleotide polymorphisms
#'   (SNPs) by position on the read when that information is available.
#'   \code{filter.snp.number} sets the maximum SNP count per locus and
#'   excludes loci above it. See [filter_snp_position_read()] and
#'   [filter_snp_number()].
#' \item \strong{Linkage disequilibrium (LD).} \code{filter.short.ld} selects
#'   a SNP per shared locus using \code{"mac"} (highest global minor allele
#'   count), \code{"random"}, \code{"first"}, \code{"middle"}, or
#'   \code{"last"}. Omit it in batch mode
#'   to skip within-locus pruning. Set \code{filter.long.ld} for longer-range
#'   pruning. The default \code{ld.method = "r2"} interprets the cutoff as
#'   squared correlation. \code{long.ld.missing = TRUE} requests the
#'   missingness-based selection in [filter_ld()].
#' \item \strong{Mixed and duplicate genomes.} In batch mode,
#'   \code{detect.mixed.genomes = TRUE} enables [detect_mixed_genomes()],
#'   with optional lower and upper \code{ind.heterozygosity.threshold}
#'   bounds. \code{detect.duplicate.genomes = TRUE} enables
#'   [detect_duplicate_genomes()]. Supply a positive \code{dup.threshold}
#'   to blacklist pairs with relative distance strictly below that cutoff.
#'   With no positive cutoff in batch mode, duplicate detection produces
#'   diagnostics. Both detection flags default to \code{FALSE} in batch mode.
#' \item \strong{Hardy-Weinberg equilibrium (HWE).} The interactive workflow
#'   displays evidence and asks for thresholds. In batch mode,
#'   \code{filter.hwe = TRUE} enables [filter_hwe()] with its defaults,
#'   including a Benjamini-Hochberg-adjusted mid-p cutoff of \code{1e-4}
#'   and a departure in at least one eligible stratum. For custom HWE
#'   settings, call [filter_hwe()] separately. The legacy
#'   \code{hw.pop.threshold} and \code{midp.threshold} options have no effect
#'   in this wrapper.
#' }
#'
#' @section Running interactively or in a script:
#' In interactive mode, inspect the diagnostic plots before choosing cutoffs.
#' Answers to prompts can replace thresholds supplied in the call. For a
#' repeatable script, set \code{interactive.filter = FALSE}, give named
#' thresholds, and set \code{random.seed} for steps that use randomness.
#' If omitted, the workflow generates a seed. It records the seed used.
#' \code{path.folder} sets the workflow's parent directory, which defaults to
#' the working directory. Call individual filters directly when you need a
#' different order or settings specific to a lower-level function.
#'
#' @section Outputs:
#' The dated workflow directory contains the recorded arguments, seed,
#' filter history and per-step audit tables, diagnostic figures, and a
#' \code{filtered} subdirectory with the final marker whitelist, marker and
#' sample blacklists when exclusions exist, filtered strata, and final
#' statistics. Use [genometranslator::list_filters()] with
#' \code{history = TRUE} to inspect filtering operations saved in the GDS.
#' Use [genometranslator::extract_individuals_metadata()] with
#' \code{blacklist = TRUE} to inspect excluded sample IDs and their reasons.
#'
#' @return A named list with \code{gds}, the filtered open GDS connection, and
#'   \code{output}, the conversion result. The latter is \code{"not selected"}
#'   when \code{output = NULL} and \code{"written"} when the conversion writes
#'   a file and returns no object. Close the GDS connection with
#'   \code{SeqArray::seqClose()} when finished.
#'
#' @examples
#' \dontrun{
#' # Work on a copy of the imported GDS. The blacklist is applied first.
#' guided <- radr::explore_genomes(
#'   data = "shark_guided.gds",
#'   blacklist.id = "blacklist.tsv",
#'   interactive.filter = TRUE
#' )
#' SeqArray::seqClose(guided$gds)
#'
#' # Example thresholds for a scripted run on another working copy.
#' screened <- radr::explore_genomes(
#'   data = "shark_batch.gds",
#'   interactive.filter = FALSE,
#'   blacklist.id = c("001", "002"),
#'   filter.individuals.missing = 0.30,
#'   filter.genotyping = 0.20,
#'   random.seed = 42L
#' )
#' genometranslator::extract_individuals_metadata(
#'   screened$gds, blacklist = TRUE
#' )
#' SeqArray::seqClose(screened$gds)
#' }
#'
#' @seealso
#' [filter_individuals()], [filter_monomorphic()], [filter_common_markers()],
#' [filter_ma()], [filter_genotyping()], [filter_ld()], [filter_hwe()],
#' [genometranslator::list_filters()].
#'
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
explore_genomes <- function(
  data,
  interactive.filter = TRUE,
  output = NULL,
  filename = NULL,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "explore_genomes",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  res <- list()

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c(
      "subsample.markers.stats",
      "filter.reproducibility", "blacklist.id", "filter.individuals.missing",
      "filter.individuals.heterozygosity", "filter.individuals.coverage.total",
      "filter.individuals.coverage.median", "filter.individuals.coverage.iqr",
      "filter.common.markers", "filter.monomorphic", "filter.ma",
      "filter.coverage", "filter.genotyping", "filter.snp.position.read",
      "filter.snp.number", "filter.short.ld", "filter.long.ld", "long.ld.missing",
      "ld.method", "detect.mixed.genomes", "ind.heterozygosity.threshold",
      "detect.duplicate.genomes", "dup.threshold",
      "filter.hwe", "hw.pop.threshold", "midp.threshold",
      "random.seed", "path.folder", "filename",
      "whitelist.markers",
      "internal"),
    deprecated = c("pop.select, pop.labels"),
    verbose = FALSE
  )
  filter.common.markers.bk <- filter.common.markers
  filter.monomorphic.bk <- filter.monomorphic
  filter.individuals.missing.bk <- filter.individuals.missing
  filter.individuals.heterozygosity.bk <- filter.individuals.heterozygosity
  filter.individuals.coverage.total.bk <- filter.individuals.coverage.total
  filter.individuals.coverage.median.bk <- filter.individuals.coverage.median
  filter.individuals.coverage.iqr.bk <- filter.individuals.coverage.iqr
  filter.ma.bk <- filter.ma
  filter.coverage.bk <- filter.coverage
  filter.genotyping.bk <- filter.genotyping
  filter.snp.position.read.bk <- filter.snp.position.read
  filter.snp.number.bk <- filter.snp.number
  filter.short.ld.bk <- filter.short.ld
  filter.long.ld.bk <- filter.long.ld
  long.ld.missing.bk <- long.ld.missing
  ld.method.bk <- ld.method
  detect.mixed.genomes.bk <- detect.mixed.genomes
  ind.heterozygosity.threshold.bk <- ind.heterozygosity.threshold
  detect.duplicate.genomes.bk <- detect.duplicate.genomes
  dup.threshold.bk <- dup.threshold
  filter.hwe.bk <- filter.hwe
  hw.pop.threshold.bk <- hw.pop.threshold
  midp.threshold.bk <- midp.threshold


  # Test
  # internal <- FALSE
  # path.folder <- NULL
  # random.seed <- NULL


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  data.type <- tryCatch(
    genometranslator::detect_genomic_format(data),
    error = function(cnd) "unsupported"
  )
  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort(
      c(
        "`explore_genomes()` requires a GDS file or object.",
        i = paste0(
          "The supplied input was detected as `", data.type, "`."
        ),
        i = paste0(
          "Import it first with an appropriate `genometranslator` ",
          "reader, then pass the resulting GDS to `explore_genomes()`."
        )
      )
    )
  }

  # Required packages ----------------------------------------------------------
  if (interactive.filter) {
    tgbase::check_package(package = "HardyWeinberg")
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)
    tgbase::check_package(package = "SNPRelate", cran = FALSE, bioc = TRUE)
  }

  # Folders---------------------------------------------------------------------
  # wf for working folder
  # radr.folder : the 01_radr folder to put most things...
  wf <- path.folder <- generate_folder(
    rad.folder = "explore_genomes",
    path.folder = path.folder,
    prefix.int = TRUE,
    append.date = TRUE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  radr.folder <- generate_folder(
    rad.folder = "radr",
    path.folder = path.folder,
    prefix.int = TRUE,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = radr.folder,
    filename = "radr_explore_genomes_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # Random seed ----------------------------------------------------------------
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }
  readr::write_lines(x = random.seed, file = file.path(radr.folder, "random.seed"))
  if (verbose) message("File written: random.seed (", random.seed,")")

  # Filter parameter file: generate --------------------------------------------
  filters.parameters <- filter_parameters(
    generate = TRUE,
    initiate = FALSE,
    update = FALSE,
    parameter.obj = NULL,
    path.folder = radr.folder,
    file.date = file.date,
    verbose = verbose)

  # GDS input ------------------------------------------------------------------
  if (data.type == "gds.file") {
    gds <- genometranslator::read_genome(data = data, verbose = verbose)
    completed <- FALSE
    on.exit(if (!completed) SeqArray::seqClose(gds), add = TRUE)
  } else {
    gds <- data
  }

  # Known exclusions precede QC and marker filters.
  if (length(blacklist.id)) {
    gds <- filter_individuals(
      data = gds,
      interactive.filter = FALSE,
      blacklist.id = blacklist.id,
      parallel.core = parallel.core,
      verbose = verbose,
      path.folder = wf,
      parameters = filters.parameters,
      internal = FALSE
    )
  }

  source <- genometranslator::extract_data_source(gds) # to know if dart data or not...

  # Filter reproducibility -----------------------------------------------------
  if ("dart" %in% source) {
    # Filter parameter file: initiate
    filters.parameters <- filter_parameters(
      generate = FALSE,
      initiate = TRUE,
      update = FALSE,
      parameter.obj = filters.parameters,
      data = gds,
      path.folder = radr.folder,
      file.date = file.date,
      verbose = verbose)

    gds <- filter_dart_reproducibility(
      data = gds,
      interactive.filter = interactive.filter,
      filter.reproducibility = filter.reproducibility,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = wf,
      internal = FALSE)
  }

  # Filter_monomorphic----------------------------------------------------------
  gds <- filter_monomorphic(
    data = gds,
    filter.monomorphic = filter.monomorphic.bk,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter common markers ------------------------------------------------------
  gds <- filter_common_markers(
    data = gds,
    filter.common.markers = filter.common.markers.bk,
    fig = TRUE,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter Whitelist------------------------------------------------------------
  gds <- filter_whitelist(
    data = gds,
    whitelist.markers = whitelist.markers,
    verbose = verbose,
    path.folder = wf,
    parameters = filters.parameters,
    biallelic = biallelic,
    internal = FALSE)

  # Genotypes metadata ---------------------------------------------------------
  # check for coverage information...
  # genotypes.meta <- extract_genotypes_metadata(gds, index.only = TRUE)
  # count.data <- "ALLELE_REF_DEPTH" %in% genotypes.meta
  # coverage.info <- "READ_DEPTH" %in% genotypes.meta
  # genotypes.meta <- NULL

  # Filter_individuals----------------------------------------------------------
  gds <- filter_individuals(
    data = gds,
    interactive.filter = interactive.filter,
    filter.individuals.missing = filter.individuals.missing.bk,
    filter.individuals.heterozygosity = filter.individuals.heterozygosity.bk,
    filter.individuals.coverage.total = filter.individuals.coverage.total.bk,
    filter.individuals.coverage.median = filter.individuals.coverage.median.bk,
    filter.individuals.coverage.iqr = filter.individuals.coverage.iqr.bk,
    parallel.core = parallel.core,
    verbose = verbose,
    path.folder = wf,
    parameters = filters.parameters,
    internal = FALSE
  )

  # Filter MAC------------------------------------------------------------------
  gds <- filter_ma(
    data = gds,
    interactive.filter = interactive.filter,
    filter.ma = filter.ma.bk,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)


  # Filter coverage-------------------------------------------------------------
  gds <- filter_coverage(
    data = gds,
    interactive.filter = interactive.filter,
    filter.coverage = filter.coverage.bk,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter genotyping---------------------------------------------------------
  gds <- filter_genotyping(
    data = gds,
    interactive.filter = interactive.filter,
    filter.genotyping = filter.genotyping.bk,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter SNP position on the read---------------------------------------------
  gds <- filter_snp_position_read(
    data = gds,
    interactive.filter = interactive.filter,
    filter.snp.position.read = filter.snp.position.read.bk,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)


  # Filter snp number per locus ------------------------------------------------
  gds <- filter_snp_number(
    data = gds,
    interactive.filter = interactive.filter,
    filter.snp.number = filter.snp.number.bk,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter Linkage disequilibrium --------------------------------------------
  continue.with.ld <- "y"

  if (interactive.filter) {
    continue.with.ld <- tgbase::question(
      x = "\n\nDo you want to filter the dataset using short distance ld  ? (y/n):",
      answer.opt = c("y", "n"))
  }

  if (continue.with.ld == "y") {
    gds <- filter_ld(
      data = gds,
      interactive.filter = interactive.filter,
      filter.short.ld = filter.short.ld.bk,
      filter.long.ld = filter.long.ld.bk,
      parallel.core = parallel.core,
      filename = NULL,
      verbose = verbose,
      long.ld.missing = long.ld.missing.bk,
      ld.method = ld.method.bk,
      parameters = filters.parameters,
      path.folder = wf,
      internal = FALSE)
  }

  # Detect mixed genomes -------------------------------------------------------
  gds <- detect_mixed_genomes(
    data = gds,
    interactive.filter = interactive.filter,
    detect.mixed.genomes = detect.mixed.genomes.bk,
    ind.heterozygosity.threshold = ind.heterozygosity.threshold.bk,
    parameters = filters.parameters,
    verbose = verbose,
    parallel.core = parallel.core,
    path.folder = wf,
    internal = FALSE)

  #Detect duplicate genomes ---------------------------------------------------
  gds <- detect_duplicate_genomes(
    data = gds,
    interactive.filter = interactive.filter,
    detect.duplicate.genomes = detect.duplicate.genomes.bk,
    dup.threshold = dup.threshold.bk,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    random.seed = random.seed,
    path.folder = wf,
    internal = FALSE)

  # Filter HWE -----------------------------------------------------------------
  if (interactive.filter || isTRUE(filter.hwe.bk)) {
    gds <- filter_hwe(
      data = gds,
      interactive.filter = interactive.filter,
      strata = NULL,
      path.folder = wf,
      verbose = verbose,
      internal = TRUE)
  }


  # Final Sync GDS -----------------------------------------------------------
  if (verbose) message("\nPreparing output files...")
  markers.meta <- genometranslator::extract_markers_metadata(
    gds = gds,
    whitelist = TRUE)

  strata <- genometranslator::extract_individuals_metadata(
    gds = gds,
    # ind.field.select = c("INDIVIDUALS", "STRATA"),
    whitelist = TRUE)

  genometranslator::sync_gds(
    gds = gds,
    samples = strata$INDIVIDUALS,
    variant.id = markers.meta$VARIANT_ID
  )

  # FINAL PREP
  # filtered data folder
  path.folder <- generate_folder(
    rad.folder = "filtered",
    path.folder = wf,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # Whitelist
  tgbase::write_tgbase_tsv(
    data = markers.meta,
    path.folder = path.folder,
    filename = "whitelist.markers",
    date = TRUE,
    internal = FALSE,
    write.message = "standard",
    verbose = verbose
  )

  # blacklist
  bl <- genometranslator::extract_markers_metadata(gds = gds, blacklist = TRUE)
  if (nrow(bl) > 0) {
    tgbase::write_tgbase_tsv(
      data = bl,
      path.folder = path.folder,
      filename = "blacklist.markers",
      date = TRUE,
      internal = FALSE,
      write.message = "standard",
      verbose = verbose
    )
  }


  # writing the blacklist of id
  blacklisted.individuals <- genometranslator::extract_individuals_metadata(
    gds = gds, blacklist = TRUE
  )
  if (nrow(blacklisted.individuals) > 0) {
    tgbase::write_tgbase_tsv(
      data = blacklisted.individuals,
      path.folder = path.folder,
      filename = "blacklist.id",
      date = TRUE,
      internal = FALSE,
      write.message = "standard",
      verbose = verbose
    )
  }

  # Generate new strata
  tgbase::write_tgbase_tsv(
    data = strata,
    path.folder = path.folder,
    filename = "strata.filtered",
    date = TRUE,
    internal = FALSE,
    write.message = "Writing the filtered strata: ",
    verbose = verbose
  )

  # Statistics after filtering -------------------------------------------------
  if (verbose) message("\nGenerating statistics after filtering")

  # Individuals and markers stats
  invisible(genometranslator::generate_stats(
    gds = gds,
    individuals = TRUE,
    markers = TRUE,
    missing = TRUE,
    heterozygosity = TRUE,
    coverage = TRUE,
    allele.coverage = TRUE,
    mac = TRUE,
    snp.position.read = TRUE,
    snp.per.locus = TRUE,
    path.folder = path.folder,
    file.date = file.date,
    parallel.core = parallel.core,
    verbose = verbose
  ))

  # Optional output translation -----------------------------------------------
  if (is.null(output)) {
    res$output <- "not selected"
  } else {
    output.filename <- filename
    if (is.null(output.filename)) output.filename <- "filtered_genome"
    output.filename <- file.path(path.folder, basename(output.filename))

    if (verbose) {
      message(
        "\nWriting requested genomic output: ",
        paste(output, collapse = ", ")
      )
    }
    res$output <- genometranslator::genome_translator(
      data = gds,
      strata = strata,
      output = output,
      parallel.core = parallel.core,
      filename = output.filename,
      verbose = verbose
    )
    if (is.null(res$output)) res$output <- "written"
  }
  genometranslator::summary_gds(gds = gds, verbose = TRUE)
  res$gds <- gds
  completed <- TRUE
  return(res)
} # End explore_genomes
