#' @name summarise_genomic_data
#' @title Summarise genomic data
#' @description Calculate summary statistics from tidy genomic data or a GDS
#' file. Statistics are calculated by population and marker and include REF and
#' ALT allele frequencies, observed and expected heterozygosity, and the
#' inbreeding coefficient (FIS).
#' @param data Tidy genomic data or a Genomic Data Structure (GDS) file or
#' object:
#' \itemize{
#' \item tidy data
#' \item Genomic Data Structure (GDS)
#' }
#'
#' \emph{How to get GDS and tidy data ?}
#' Use \code{\link[genometranslator]{read_genome}} to import supported formats
#' and \code{\link[genometranslator]{tidy_genome}} when a tidy table is needed.
# @param filename (optional) Name of the file written to the working directory.
#' @param path.folder (path, optional) By default will print results in the working directory.
#' Default: \code{path.folder = NULL}.
#' @param digits (integer, optional). Default: \code{digits = 4}.
#' @inheritParams radr_common_arguments
#' @rdname summarise_genomic_data
#' @export
# @keywords internal

summarise_genomic_data <- function(
  data,
  path.folder = NULL,
  digits = 4,
  verbose = TRUE
) {
  # Cleanup-------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("###################### summarise_genomic_data completed #######################\n"), add = TRUE)

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file missing")


  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    rad.folder = "summary_stats",
    path.folder = path.folder,
    prefix.int = FALSE,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # Detect format --------------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)
  if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)

    message("Importing gds file...")
    if (data.type == "gds.file") {
      data <- genometranslator::read_genome(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
    }

    tidy.data <- genometranslator::extract_genotypes_metadata(
      gds = data,
      genotypes.meta.select = c("MARKERS", "COL", "INDIVIDUALS", "POP_ID", "ALT_DOSAGE"))

    if (length(tidy.data) == 0) {
      tidy.data <- genometranslator::tidy_genome(data = data)
    }
    SeqArray::seqClose(data)
    data <- tidy.data
    tidy.data <- NULL

  } else {#tidy.data
    if (is.vector(data)) data <- genometranslator::read_genome(data = data, import.metadata = TRUE)
  }
  if (verbose) message("Summarizing...")
  # the function
  summarise_group <- function(x, maf = c("global", "local"), digits = 6) {
    maf <- match.arg(arg = maf, choices = c("global", "local"))
    x %<>%
      dplyr::summarise(
        N = as.numeric(length(unique(INDIVIDUALS))),
        PP = as.numeric(length(ALT_DOSAGE[ALT_DOSAGE == 0])),
        PQ = as.numeric(length(ALT_DOSAGE[ALT_DOSAGE == 1])),
        QQ = as.numeric(length(ALT_DOSAGE[ALT_DOSAGE == 2]))
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(
        FREQ_REF = ((PP*2) + PQ)/(2*N),
        FREQ_ALT = ((QQ*2) + PQ)/(2*N),
        HET_O = PQ/N,
        HET_E = 2 * FREQ_REF * FREQ_ALT,
        FIS = dplyr::if_else(HET_O == 0, 1, round(((HET_E - HET_O) / HET_E), digits)),
        PP = NULL, PQ = NULL, QQ = NULL
      )

    if (maf == "global") {
      x %<>% dplyr::rename(MAF_GLOBAL = FREQ_ALT)
    } else {
      x %<>% dplyr::rename(MAF_LOCAL = FREQ_ALT)
    }

    return(x)
  }#End summarise_group

  data  %<>% dplyr::filter(!is.na(ALT_DOSAGE))
  # ms: markers stats

  ms <- data %>%
    dplyr::group_by(MARKERS) %>%
    summarise_group(x = ., maf = "global", digits = digits)

  # mps: markers pop stats
  mps <- data %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    summarise_group(x = ., maf = "local", digits = digits)

  # ps: pop stats
  ps <- dplyr::bind_cols(
    dplyr::group_by(data, POP_ID) %>% dplyr::summarise(N = length(unique(INDIVIDUALS))),
    dplyr::group_by(mps, POP_ID) %>%
      dplyr::summarise(
        .data = .,
        dplyr::across(
          .cols = c(N, FREQ_REF, MAF_LOCAL, HET_O, HET_E),
          .fns = \(x) mean(x, na.rm = TRUE)
        ),
        .groups = "keep"
      ) %>%
      dplyr::rename(N_MEAN = N)
  ) %>%
    dplyr::mutate(
      POP_ID1 = NULL,
      FIS = dplyr::if_else(HET_O == 0, 1, round(((HET_E - HET_O) / HET_E), digits))
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(FREQ_REF, MAF_LOCAL, HET_O, HET_E, FIS),
        .fns  = \(x) round(x, digits = digits)
      )
    )

  # writting the results
  tgbase::write_tgbase_tsv(
    data = ms,
    path.folder = path.folder,
    filename = "summary.stats.markers",
    date = TRUE,
    internal = FALSE,
    write.message = "standard",
    verbose = verbose
  )

  tgbase::write_tgbase_tsv(
    data = mps,
    path.folder = path.folder,
    filename = "summary.stats.markers.pop",
    date = TRUE,
    internal = FALSE,
    write.message = "standard",
    verbose = verbose
  )

  tgbase::write_tgbase_tsv(
    data = ps,
    path.folder = path.folder,
    filename = "summary.stats.pop",
    date = TRUE,
    internal = FALSE,
    write.message = "standard",
    verbose = verbose
  )

  return(res = list(summary.stats.pop = ps,
                    summary.stats.markers.pop = mps,
                    summary.stats.markers = ms))
}#End summary_stats_vcf_tidy
