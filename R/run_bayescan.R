# classify_bayescan_selection --------------------------------------------------
# Central scientific classification used by single runs and subsampling
# summaries. The q-value establishes support for selection; ALPHA only
# determines its direction.
classify_bayescan_selection <- function(alpha, q.value, fdr = 0.05) {
  selection <- dplyr::case_when(
    !is.na(q.value) & q.value <= fdr & !is.na(alpha) & alpha > 0 ~
      "diversifying",
    !is.na(q.value) & q.value <= fdr & !is.na(alpha) & alpha < 0 ~
      "balancing",
    !is.na(q.value) & !is.na(alpha) ~ "neutral",
    .default = NA_character_
  )

  factor(
    selection,
    levels = c("neutral", "balancing", "diversifying")
  )
}

# run_bayescan -----------------------------------------------------------------
#' @name run_bayescan
#' @title Run a BayeScan genome scan
#' @description Run BayeScan, import its results, classify loci using a chosen
#' false-discovery-rate threshold, and generate diagnostic tables, figures,
#' whitelists, and blacklists.
#'
#' \strong{Function highlights:}
#'
#' \enumerate{
#' \item \strong{integrated and seamless pipeline: } generate \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' files within radr and run \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan} inside R!
#' \item \strong{unbalanced sampling sites impact: } measure and verify genome scan accurary in unbalanced sampling design with subsampling related arguments.
#' \item \strong{SNP linkage: } detect automatically the presence of multiple SNPs on the same locus and
#' measure/verify accuracy of genome scan within locus.
#' \item \strong{summary tables and visualization: } the function generate summary tables and plots of genome scan.
#' \item \strong{whitelists and blacklists} of markers under different selection identity are automatically generated !
#' }
#'
#' This function requires a working
#' \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan} program installed
#' on the computer (\href{http://cmpg.unibe.ch/software/BayeScan/download.html}{install instructions}).
#' For UNIX machines, please install the 64bits version.

#' @param data BayeScan input file, GDS file or open GDS object, or tidy genomic
#' data.
#' \enumerate{
#' \item Path to BayeScan input file.
#' Generate this with \code{\link[genometranslator]{write_bayescan}}.
#' \item GDS file or open GDS object. Population assignments are obtained from
#' the individual metadata stored in the GDS.
#' \item Tidy genomic data file or object. Tidy data and GDS input can be used
#' with \code{subsample} and \code{iteration.subsample}.
#' }
#' @param n Integer. Number of output iterations.
#' Default: \code{n = 5000}.
#' @param thin Integer. Thinning interval.
#' Default: \code{thin = 10}.
#' @param nbp Integer. Number of pilot runs.
#' Default: \code{nbp = 20}.
#' @param pilot Integer. Length of each pilot run.
#' Default: \code{pilot = 5000}.
#' @param burn Integer. Burn-in length.
#' Default: \code{burn = 50000}.
#' @param pr_odds Numeric. Prior odds for the neutral model. A \code{pr_odds = 10}
#' indicates that the neutral model is 10 times more likely than the
#' model with selection. Larger values make the analysis more conservative.
#' This argument is required and has no default.
#' @param fdr Numeric false-discovery-rate threshold used to classify loci from
#' BayeScan q-values. A locus is classified as diversifying or balancing only
#' when its q-value is at or below this threshold; the sign of \code{ALPHA}
#' determines the direction. Default: \code{fdr = 0.05}.

#' @param subsample Integer, proportion, or \code{"min"}.
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset. With \code{subsample = "min"}, the
#' minimum number of individual/population found in the data is used automatically.
#' Default: \code{subsample = NULL}.

#' @param iteration.subsample Integer. Number of repeated subsampling iterations.
#' subsampling.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.
#' Default: \code{iteration.subsample = 1}.

#' @param parallel.core Integer. Number of CPU cores available to BayeScan.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param bayescan.path Character. Full path to the BayeScan executable. When
#' \code{NULL}, \code{check_bayescan()} searches the system \code{PATH} and the
#' environment specified by \code{conda.env}.
#' Default: \code{bayescan.path = NULL}. See details.
#' @param conda.env Conda environment name or prefix used when discovering
#' BayeScan. Default: \code{conda.env = "genomics"}.
#' @param verbose Logical. Display progress messages.
#' Default: \code{verbose = TRUE}.
#' @param ... Additional arguments, including \code{path.folder} for the parent
#' results directory.


#' @rdname run_bayescan
#' @export
#' @return For specific \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' output files, see \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' documentation, please read the manual.
#'
#' radr::run_bayescan outputs without subsampling:
#'
#' \enumerate{
#' \item \code{bayescan}: dataframe with results of BayeScan analysis.
#' \item \code{selection.summary}: dataframe showing the number of markers in the different group of selections and model choice.
#' \item \code{whitelist.markers.positive.selection}: Whitelist of markers under diversifying selection and common in all iterations.
#' \item \code{whitelist.markers.neutral.selection}: Whitelist of neutral markers and common in all iterations.
#' \item \code{whitelist.markers.neutral.positive.selection}: Whitelist of neutral markers and markers under diversifying selection and common in all iterations.
#' \item \code{blacklist.markers.balancing.selection}: Blacklist of markers under balancing selection and common in all iterations.
#' \item \code{markers.dictionary}: BayeScan use integer for MARKERS info. In this dataframe, the corresponding values used inside the function.
#' \item \code{pop.dictionary}: BayeScan use integer for STRATA info. In this dataframe, the corresponding values used inside the function.
#' \item \code{bayescan.plot}: plot showing markers Fst and model choice.
#'
#' Additionnally, if multiple SNPs/locus are detected the object will also have:
#' \item \code{accurate.locus.summary}: dataframe with the number of accurate locus and the selection types.
#' \item \code{whitelist.accurate.locus}: whitelist of accurate locus.
#' \item \code{blacklist.not.accurate.locus}: blacklist of not accurate locus.
#' \item \code{accuracy.snp.number}: dataframe with the number of SNPs per locus and the count of accurate/not accurate locus.
#' \item \code{accuracy.snp.number.plot}: the plot showing the proportion of accurate/not accurate locus in relation to SNPs per locus.
#' \item \code{not.accurate.summary}: dataframe summarizing the number of not accurate locus with selection type found on locus.
#' }
#'
#' radr::run_bayescan outputs WITH subsampling:
#'
#' \enumerate{
#' \item \code{subsampling.individuals}: dataframe with indivuals subsample id and random seed number.
#' \item \code{bayescan.all.subsamples}: long dataframe with combined iterations of bayescan results.
#' \item \code{selection.accuracy}: dataframe with all markers with selection grouping and number of times observed throughout iterations.
#' \item \code{accurate.markers}: dataframe with markers attributed the same selection grouping in all iterations.
#' \item \code{accuracy.summary}: dataframe with a summary of accuracy of selection grouping.
#' \item \code{bayescan.summary}: dataframe with mean value, averaged accross iterations.
#' \item \code{bayescan.summary.plot}: plot showing markers Fst and model choice.
#' \item \code{selection.summary}: dataframe showing the number of markers in the different group of selections and model choice.
#' \item \code{whitelist.markers.positive.selection}: Whitelist of markers under diversifying selection and common in all iterations.
#' \item \code{whitelist.markers.neutral.selection}: Whitelist of neutral markers and common in all iterations.
#' \item \code{blacklist.markers.balancing.selection}: Blacklist of markers under balancing selection and common in all iterations.
#' \item \code{whitelist.markers.neutral.positive.selection}: Whitelist of neutral markers and markers under diversifying selection and common in all iterations.
#' \item \code{whitelist.markers.without.balancing.positive}:
#' Whitelist of all original markers with markers under balancing selection and directional selection removed.
#' The markers that remains are the ones to use in population structure analysis.
#' }
#'
#' Other files are present in the folder and subsampling folder.


#' @examples
#' \dontrun{
#' # library(radr)
#' # get a tidy data frame and a bayescan file with genometranslator::genome_translator:
#' # to run with a vcf haplotype file
#' data <- genometranslator::genome_translator(
#'     data = "batch_1.haplotypes.vcf",
#'     strata = "../../02_project_info/strata.stacks.TL.tsv",
#'     whitelist.markers = "whitelist.filtered.markers.tsv",
#'     blacklist.id = "blacklist.id.tsv",
#'     output = "bayescan",
#'     filename = "bayescan.haplotypes"
#'     )
#' # to run BayeScan:
#' scan.pops <- radr::run_bayescan(
#'     data = "bayescan.haplotypes.txt",
#'     pr_odds = 1000
#'     )
#'
#' # This will use the default values for argument: n, thin, nbp, pilot and burn.
#' # The number of CPUs will be the number available - 1 (the default).
#'
#' # To test the impact of unbalance sampling run BayeScan with subsampling,
#' # for this, you need to feed the function the tidy data frame generated above
#' # with genometranslator::genome_translator:
#' scan.pops.sub <- radr::run_bayescan(
#'     data = data$tidy.data,
#'     pr_odds = 1000,
#'     subsample = "min",
#'     iteration.subsample = 10
#'     )
#'
#' # This will run BayeScan 10 times, and for each iteration, the number of individuals
#' # sampled in each pop will be equal to the minimal number found in the pops
#' # (e.g. pop1 N = 36, pop2 N = 50 and pop3 N = 15, the subsampling will use 15
#' # individuals in each pop, taken randomly.
#' # You can also choose a specific subsample value with the argument.
#' }

#' @details
#' \strong{Selection classification:}
#' A q-value at or below \code{fdr} is required before a locus is classified as
#' a selection candidate. Among supported candidates, positive \code{ALPHA}
#' indicates diversifying selection and negative \code{ALPHA} indicates
#' balancing or purifying selection. Loci above the FDR threshold are neutral
#' regardless of the sign of \code{ALPHA}.
#'
#' \strong{subsampling:}
#' During subsampling the function will automatically remove monomorphic
#' markers that are generated by the removal of some individuals. Also, common markers
#' between all populations are also automatically detected. Consequently, the number of
#' markers will change throughout the iterations. The nice thing about the function
#' is that since everything is automated there is less chance of making an error...
#'
#' \strong{SNPs data set: }
#' You should not run BayeScan with SNPs data set that have multiple SNPs on the
#' same LOCUS. Instead, run genometranslator::genome_translator using the \code{snp.ld}
#' argument to keep only one SNP on the locus. Or run the function by first converting
#' an haplotype vcf or if your RAD dataset was produced by STACKS, use the
#' \code{batch_x.haplotypes.tsv} file! If the function detect multiple SNPs on
#' the same locus, accuracy will be measured automatically.
#'
#' \strong{UNIX install: } I like to transfer the \emph{BayeScan2.1_linux64bits}
#' (for Linux) or the \emph{BayeScan2.1_macos64bits} (for MACOs) in \code{/usr/local/bin}
#' and change it's name to \code{bayescan}. Too complicated ? and you've just
#' downloaded the last BayeScan version, I would try this :
#' \code{bayescan.path = "/Users/thierry/Downloads/BayeScan2.1/binaries/BayeScan2.1_macos64bits"}
#'
#' Make sure to give permission: \code{sudo chmod 777 /usr/local/bin/bayescan}


#' @seealso
#' \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}

#' @references Foll, M and OE Gaggiotti (2008) A genome scan method to identify
#' selected loci appropriate
#' for both dominant and codominant markers: A Bayesian perspective.
#' Genetics 180: 977-993

#' @references Foll M, Fischer MC, Heckel G and L Excoffier (2010)
#' Estimating population structure from
#' AFLP amplification intensity. Molecular Ecology 19: 4638-4647

#' @references Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP
#' genome scans detect
#' local adaptation in high-altitude populations of a small rodent (Microtus arvalis).
#' Molecular Ecology 20: 1450-1462

run_bayescan <- function(
  data,
  n = 5000,
  thin = 10,
  nbp = 20,
  pilot = 5000,
  burn = 50000,
  pr_odds,
  fdr = 0.05,
  subsample = NULL,
  iteration.subsample = 1,
  parallel.core = parallel::detectCores() - 1,
  bayescan.path = NULL,
  conda.env = "genomics",
  verbose = TRUE,
  ...
) {

  # test
  # n = 5000
  # thin = 10
  # nbp = 20
  # pilot = 5000
  # burn = 50000
  # subsample = NULL
  # iteration.subsample = 1
  # parallel.core = parallel::detectCores() - 1
  # bayescan.path = "/usr/local/bin/bayescan"

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "run_bayescan",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  res <- list() # return results in this list

  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "internal"),
    verbose = FALSE
  )

  # Check BayeScan installation -----------------------------------------------
  bayescan.path <- check_bayescan(
    bayescan.path = bayescan.path,
    conda.env = conda.env,
    verbose = verbose
  )

  if (missing(data)) rlang::abort("Input file missing")
  if (missing(pr_odds)) rlang::abort("Prior odds for the neutral model is missing.
                             No shortcut with default here, sorry.
                             Please read the BayeScan manual...")

  if (!is.numeric(fdr) || length(fdr) != 1L || is.na(fdr) ||
      !is.finite(fdr) || fdr <= 0 || fdr >= 1) {
    rlang::abort("fdr must be a single number greater than 0 and less than 1")
  }

  # Results folder and reproducibility log ------------------------------------
  path.folder <- tgbase::generate_folder(
    folder = "run_bayescan",
    path.folder = path.folder,
    internal = internal,
    file.date = file.date,
    prefix.int = TRUE,
    verbose = verbose
  )

  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_run_bayescan_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # Input format ---------------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data = data)

  # For a direct GDS run, generate the BayeScan input in the results folder.
  if (is.null(subsample) &&
      data.type %in% c("SeqVarGDSClass", "gds.file")) {
    bayescan.prefix <- file.path(path.folder, "radr_bayescan_input")
    genometranslator::write_bayescan(
      data = data,
      filename = bayescan.prefix,
      parallel.core = parallel.core,
      verbose = verbose
    )
    data <- paste0(bayescan.prefix, "_bayescan.txt")
  }

  # Subsampling ----------------------------------------------------------------
  if (!is.null(subsample)) {
    message("Subsampling: selected")
    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
      data <- genometranslator::tidy_genome(
        data = data,
        pop.id = FALSE,
        calibrate.alleles = FALSE,
        parallel.core = parallel.core,
        close.gds = identical(data.type, "gds.file")
      )
    } else if (data.type %in% c("fst.file", "parquet.file", "tbl_df", "data.frame")) {
      data <- genometranslator::read_genome(data = data, import.metadata = TRUE)
    } else if (is.vector(data)) {
      rlang::abort(
        "Subsampling requires a GDS file/object or tidy genomic data."
      )
    }

    { # Validate tidy data
      columns.tidy <- colnames(data)
      want <- c("GT_VCF", "GT_VCF_NUC", "GT", "ALT_DOSAGE")
      want.check <- TRUE %in% (unique(want %in% columns.tidy))
      want.more <- c("MARKERS", "INDIVIDUALS", "STRATA")
      want.more.check <- isTRUE(unique(want.more %in% columns.tidy))
      is.tidy <- isTRUE(unique(c(want.check, want.more.check)))
      if (!is.tidy) {
        rlang::abort(
          "Subsampling requires MARKERS, INDIVIDUALS, STRATA and genotypes."
        )
      }
    }

    ind.pop.df <- dplyr::distinct(.data = data, STRATA, INDIVIDUALS)

    # Print some statistics ----------------------------------------------------
    strata.stats <- ind.pop.df %>%
      dplyr::group_by(STRATA) %>%
      dplyr::tally(.) %>%
      dplyr::mutate(STRATA = stringi::stri_join(STRATA, n, sep = " = "))

    n.pop <- dplyr::n_distinct(ind.pop.df$STRATA)
    n.ind <- dplyr::n_distinct(ind.pop.df$INDIVIDUALS)
    message("Number of populations: ", n.pop)
    message("Number of individuals: ", n.ind)
    message("Number of ind/pop:\n", stringi::stri_join(strata.stats$STRATA, collapse = "\n"))
    message("Number of markers: ", dplyr::n_distinct(data$MARKERS))


    if (subsample == "min") {
      subsample <- ind.pop.df %>%
        dplyr::group_by(STRATA) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(n) %>%
        purrr::flatten_int(.)
      message("\nSubsample used: ", subsample)
    }
    subsample.list <- purrr::map(
      .x = 1:iteration.subsample,
      .f = subsampling_data,
      ind.pop.df = ind.pop.df,
      subsample = subsample
    )
    # keep track of subsampling individuals and write to directory
    subsampling.individuals <- dplyr::bind_rows(subsample.list)
    readr::write_tsv(
      x = subsampling.individuals,
      file = file.path(path.folder, "radr_bayescan_subsampling_individuals.tsv"),
      col_names = TRUE,
      append = FALSE
    )
    res$subsampling.individuals <- subsampling.individuals
  } else {
    iteration.subsample <- 1
  }

  # Run BayeScan iterations-----------------------------------------------------
  if (is.null(subsample)) {
    res <- bayescan_one(
      data = data,
      n = n,
      thin = thin,
      nbp = nbp,
      pilot = pilot,
      burn = burn,
      pr_odds = pr_odds,
      fdr = fdr,
      parallel.core = parallel.core,
      path.folder = path.folder,
      file.date = file.date,
      bayescan.path = bayescan.path
    )
  } else {# iterations
    subsample.bayescan <- purrr::map(
      .x = subsample.list,
      .f = bayescan_one,
      data = data,
      n = n,
      thin = thin,
      nbp = nbp,
      pilot = pilot,
      burn = burn,
      pr_odds = pr_odds,
      fdr = fdr,
      subsample = subsample,
      iteration.subsample = iteration.subsample,
      parallel.core = parallel.core,
      path.folder = path.folder,
      file.date = file.date,
      bayescan.path = bayescan.path
    )
    # Manage subsampling results -----------------------------------------------
    cat("\n\n#######################################################################\n")
    message("Summarizing subsampling results...")
    res$bayescan.all.subsamples <- purrr::map_df(subsample.bayescan, "bayescan") %>%
      dplyr::select(-BAYESCAN_MARKERS)
    readr::write_tsv(
      x = res$bayescan.all.subsamples,
      file = file.path(path.folder, "bayescan.all.subsamples.tsv"),
      col_names = TRUE,
      append = FALSE
    )

    iteration.number <- dplyr::n_distinct(res$bayescan.all.subsamples$ITERATIONS)

    # keep only markers present for all iterations
    markers.summary <- dplyr::ungroup(res$bayescan.all.subsamples) %>%
      dplyr::select(MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::tally(.)

    markers.whitelist <- dplyr::filter(markers.summary, n == iteration.number) %>%
      dplyr::distinct(MARKERS)

    markers.all.iterations <- nrow(markers.whitelist)
    total.unique.markers <- dplyr::n_distinct(markers.summary$MARKERS)
    proportion.keeper <- round(markers.all.iterations / total.unique.markers, 2)
    message("BayeScan subsampling summary: ")
    message("    number of unique markers: ", total.unique.markers)
    message("    keeping markers common in all iterations: ", markers.all.iterations, " (= ", proportion.keeper, ")")

    bayescan.all.subsamples.filtered <- dplyr::left_join(
      markers.whitelist, res$bayescan.all.subsamples, by = "MARKERS")

    res$selection.accuracy <- bayescan.all.subsamples.filtered %>%
      dplyr::group_by(MARKERS, SELECTION) %>%
      dplyr::tally(.)

    readr::write_tsv(
      x = res$selection.accuracy,
      file = file.path(path.folder, "selection.accuracy.tsv"),
      col_names = TRUE,
      append = FALSE
    )

    res$accurate.markers <- dplyr::ungroup(res$selection.accuracy) %>%
      dplyr::filter(n == iteration.number) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::select(-n)

    readr::write_tsv(
      x = res$accurate.markers,
      file = file.path(path.folder, "accurate.markers.tsv"),
      col_names = TRUE,
      append = FALSE
    )

    accurate.markers.summary <- res$accurate.markers %>%
      dplyr::group_by(SELECTION) %>%
      dplyr::tally(.)

    accurate.markers.number <- nrow(res$accurate.markers)

    selection_count <- function(selection) {
      value <- accurate.markers.summary$n[
        accurate.markers.summary$SELECTION == selection
      ]
      if (length(value) == 0L) 0L else value[[1L]]
    }

    res$accuracy.summary <- tibble::tibble(
      total = total.unique.markers,
      `found in all iterations` = markers.all.iterations,
      `not accurate` = markers.all.iterations - accurate.markers.number,
      `accurate` = accurate.markers.number,
      `accurate + neutral` = selection_count("neutral"),
      `accurate + balancing` = selection_count("balancing"),
      `accurate + diversifying` = selection_count("diversifying")
    ) %>%
      tidyr::pivot_longer(
        data = .,
        cols = tidyselect::everything(),
        names_to = "ACCURACY_MARKERS",
        values_to = "N"
      ) %>%
      dplyr::mutate(PROP = N / total.unique.markers)

    readr::write_tsv(
      x = res$accuracy.summary,
      file = file.path(path.folder, "accuracy.summary.tsv"),
      col_names = TRUE,
      append = FALSE
    )
    res$bayescan.summary <- dplyr::left_join(
      dplyr::select(res$accurate.markers, MARKERS), res$bayescan.all.subsamples, by = "MARKERS") %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise_if(.tbl = ., .predicate = is.numeric, .funs = mean) %>%
      dplyr::select(-ITERATIONS) %>%
      # Grouping des groupes LOG10_PO and Quantile of FST
      dplyr::mutate(
        SELECTION = classify_bayescan_selection(
          alpha = ALPHA,
          q.value = Q_VALUE,
          fdr = fdr
        ),
        PO_GROUP = factor(
          dplyr::if_else(LOG10_PO > 2, "decisive",
                         dplyr::if_else(LOG10_PO > 1.5, "very strong",
                                        dplyr::if_else(LOG10_PO > 1, "strong",
                                                       dplyr::if_else(LOG10_PO > 0.5, "substantial", "no evidence")))),
          levels = c("no evidence","substantial","strong","very strong","decisive"), ordered = TRUE)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(
        FST_GROUP = dplyr::ntile(FST, 5),
        FST_GROUP = dplyr::if_else(FST_GROUP == 1, "0-20%",
                                   dplyr::if_else(FST_GROUP == 2,  "20-40%",
                                                  dplyr::if_else(FST_GROUP == 3, "40-60%",
                                                                 dplyr::if_else(FST_GROUP == 4, "60-80%", "80-100%"))))
      ) %>%
      dplyr::arrange(FST)

    readr::write_tsv(
      x = res$bayescan.summary,
      file = file.path(path.folder, "bayescan.summary.tsv"),
      col_names = TRUE,
      append = FALSE
    )

    res$bayescan.summary.plot <- plot_bayescan(res$bayescan.summary)
    ggplot2::ggsave(
      filename = file.path(path.folder, "bayescan.summary.plot.pdf"),
      plot = res$bayescan.summary.plot,
      width = 30, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE)

    res$selection.summary <- res$bayescan.summary %>%
      dplyr::group_by(SELECTION, PO_GROUP) %>%
      dplyr::tally(.) %>%
      dplyr::rename(MARKERS = n)

    readr::write_tsv(
      x = res$selection.summary,
      file = file.path(path.folder, "selection.summary.tsv"),
      col_names = TRUE,
      append = FALSE
    )

    # Generating blacklists and whitelists of all iterations -------------------
    message("Generating blacklist and whitelists for all iterations")
    all.markers <- dplyr::distinct(markers.summary, MARKERS)

    # positive  ----------------------------------------------------------------
    res$whitelist.markers.positive.selection <- res$bayescan.summary %>%
      dplyr::filter(SELECTION == "diversifying" & PO_GROUP != "no evidence") %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::arrange (MARKERS)

    if (nrow(res$whitelist.markers.positive.selection) > 0) {
      readr::write_tsv(
        x = res$whitelist.markers.positive.selection,
        file = file.path(path.folder, "whitelist.markers.positive.selection.tsv"))
      positive <- TRUE
      message("    whitelist positive/directional selection: generated")
    } else {
      message("    whitelist positive/directional selection: not generated")
      positive <- FALSE
    }

    # neutral ------------------------------------------------------------------
    res$whitelist.markers.neutral.selection <- res$bayescan.summary %>%
      dplyr::filter(SELECTION == "neutral") %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::arrange (MARKERS)
    if (nrow(res$whitelist.markers.neutral.selection) > 0) {
      readr::write_tsv(
        x = res$whitelist.markers.neutral.selection,
        file = file.path(path.folder, "whitelist.markers.neutral.selection.tsv"))
      neutral <- TRUE
      message("    whitelist neutral selection: generated")
    } else {
      message("    whitelist neutral selection: not generated")
      neutral <- FALSE
    }

    # Whitelist neutral and positive -------------------------------------------
    if (neutral && positive) {
      res$whitelist.markers.neutral.positive.selection <- res$bayescan.summary %>%
        dplyr::filter(SELECTION == "neutral" | (SELECTION == "diversifying" & PO_GROUP != "no evidence")) %>%
        dplyr::distinct(MARKERS) %>%
        dplyr::arrange (MARKERS)
      readr::write_tsv(
        x = res$whitelist.markers.neutral.positive.selection,
        file = file.path(path.folder, "whitelist.markers.neutral.positive.selection.tsv"))
      message("    whitelist neutral and positive/directional selections: generated")
    } else {
      message("    whitelist neutral and positive/directional selections: not generated")
    }

    # blacklist of balancing selected markers-----------------------------------
    res$blacklist.markers.balancing.selection <- res$bayescan.summary %>%
      dplyr::filter(SELECTION == "balancing") %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::arrange (MARKERS)

    if (nrow(res$blacklist.markers.balancing.selection) > 0) {
      readr::write_tsv(
        x = res$blacklist.markers.balancing.selection,
        file = file.path(path.folder, "blacklist.markers.balancing.selection.tsv"))
      balancing <- TRUE
      message("    blacklist balancing selection: generated")
    } else {
      message("    blacklist balancing selection: not generated")
      balancing <- FALSE
    }


    # whitelist without balancing and positive ---------------------------------
    if (neutral && positive && balancing) {
      res$whitelist.markers.without.balancing.positive <- dplyr::anti_join(
        all.markers, res$blacklist.markers.balancing.selection, by = "MARKERS") %>%
        dplyr::anti_join(res$whitelist.markers.positive.selection, by = "MARKERS")
      readr::write_tsv(
        x = res$whitelist.markers.without.balancing.positive,
        file = file.path(path.folder, "whitelist.markers.without.balancing.positive.tsv"))
      message("    whitelist without balancing and positive selection: generated")
    }

    # if no positive

    if (neutral && balancing && !positive) {
      res$whitelist.markers.without.balancing.positive <- dplyr::anti_join(
        all.markers, res$blacklist.markers.balancing.selection, by = "MARKERS")
      readr::write_tsv(
        x = res$whitelist.markers.without.balancing.positive,
        file = file.path(path.folder, "whitelist.markers.without.balancing.positive.tsv"))
      message("    whitelist without balancing and positive selection: generated")
    }
  }# End

  return(res)
}# end bayescan


# internal function ------------------------------------------------------------

# subsampling_data --------------------------------------------------------------
#' @title subsampling data
#' @description subsampling data
#' @rdname subsampling_data
#' @export
#' @keywords internal


subsampling_data <- function(
  iteration.subsample = 1,
  ind.pop.df = NULL,
  subsample = NULL,
  random.seed = NULL
) {
  # message(paste0("Creating data subsample: ", iteration.subsample))
  if (is.null(subsample)) {
    subsample.select <- ind.pop.df %>%
      dplyr::mutate(SUBSAMPLE = rep(iteration.subsample, n()))
  } else {

    # Set seed for sampling reproducibility
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
      set.seed(random.seed)
    } else {
      set.seed(random.seed)
    }

    if (subsample > 1) {# integer
      subsample.select <- ind.pop.df %>%
        dplyr::group_by(STRATA) %>%
        dplyr::sample_n(tbl = ., size = subsample, replace = FALSE)# sampling individuals for each pop
    }
    if (subsample < 1) { # proportion
      subsample.select <- ind.pop.df %>%
        dplyr::group_by(STRATA) %>%
        dplyr::sample_frac(tbl = ., size = subsample, replace = FALSE)# sampling individuals for each pop
    }
    subsample.select <- subsample.select %>%
      dplyr::mutate(
        SUBSAMPLE = rep(iteration.subsample, n()),
        RANDOM_SEED_NUMBER = rep(random.seed, n())
      ) %>%
      dplyr::arrange(STRATA, INDIVIDUALS) %>%
      dplyr::ungroup(.)
  }
  return(subsample.select)
} # End subsampling function

# bayescan_one --------------------------------------------------------------
#' @title bayescan one iteration
#' @description bayescan_one
#' @rdname bayescan_one
#' @export
#' @keywords internal


bayescan_one <- function(
  x = NULL,
  data,
  n = 5000,
  thin = 10,
  nbp = 20,
  pilot = 5000,
  burn = 50000,
  pr_odds,
  fdr = 0.05,
  subsample = NULL,
  iteration.subsample = 1,
  parallel.core = parallel::detectCores() - 1,
  path.folder,
  file.date,
  bayescan.path = "/usr/local/bin/bayescan"
) {
  res <- list()
  if (!is.null(subsample)) {
    # x <- subsample.list[[1]] # test
    subsample.id <- unique(x$SUBSAMPLE)
    message("\nBayeScan, subsample: ", subsample.id, "\n")
    path.folder.subsample <- stringi::stri_join(path.folder, "/bayescan_subsample_", subsample.id)
    dir.create(file.path(path.folder.subsample))
    folder.message <- stringi::stri_join("bayescan_subsample_", subsample.id)
    message("Subsampling folder created: ", folder.message)
  } else {
    path.folder.subsample <- path.folder
  }
  output.folder <- stringi::stri_join("-od ", path.folder.subsample)
  # output.folder <- path.folder.subsample

  log.file <- stringi::stri_join(path.folder.subsample, "/radr_bayescan_", file.date,".log")
  message("For progress, look in the log file: radr_bayescan_", file.date,".log")

  # arguments -------------------------------------
  all.trace <- "-all_trace "
  parallel.core.bk <- parallel.core
  parallel.core <- stringi::stri_join("-threads ", parallel.core)
  n <- stringi::stri_join("-n ", n)
  thin <- stringi::stri_join("-thin ", thin)
  nbp <- stringi::stri_join("-nbp ", nbp)
  pilot <- stringi::stri_join("-pilot ", pilot)
  burn <- stringi::stri_join("-burn ", burn)
  pr.odds <- stringi::stri_join("-pr_odds ", pr_odds)

  if (!is.null(subsample)) {
    # Keep only the subsample
    bayescan.filename <- stringi::stri_join(
      "radr_bayescan_subsample_", subsample.id)
    bayescan.data <- dplyr::semi_join(
      data,
      x,
      by = c("STRATA", "INDIVIDUALS")
    )
    bayescan.sub <- genometranslator::write_bayescan(
      data = bayescan.data,
      parallel.core = parallel.core.bk,
      filename = bayescan.filename)
    x <- NULL #unused object
    data <- paste0(bayescan.filename, "_bayescan.txt")
  }

  # Moving input file in folder
  data <- normalizePath(data, mustWork = TRUE)
  new.data <- file.path(path.folder.subsample, basename(data))
  if (!identical(
    normalizePath(data, mustWork = FALSE),
    normalizePath(new.data, mustWork = FALSE)
  )) {
    message("Copying input BayeScan file in folder")
    file.copy(from = data, to = new.data, overwrite = TRUE)
  }

  # moving dictionary files ----------------------------------------------------
  if (!is.null(subsample)) {
    pop.dictionary <- bayescan.sub$pop.dictionary
    markers.dictionary <- bayescan.sub$markers.dictionary
  } else {
    pop.dic.file <- sub("\\.txt$", "_pop_dictionary.tsv", data)
    markers.dic.file <- sub("\\.txt$", "_markers_dictionary.tsv", data)

    if (file.exists(pop.dic.file) && file.exists(markers.dic.file)) {
      pop.dictionary <- readr::read_tsv(
        file = pop.dic.file,
        col_types = "ci")
      markers.dictionary <- readr::read_tsv(
        file = markers.dic.file,
        col_types = "ci")
      pop.dic.destination <- file.path(
        path.folder.subsample,
        basename(pop.dic.file)
      )
      markers.dic.destination <- file.path(
        path.folder.subsample,
        basename(markers.dic.file)
      )
      if (!identical(pop.dic.file, pop.dic.destination)) {
        file.copy(pop.dic.file, pop.dic.destination, overwrite = TRUE)
      }
      if (!identical(markers.dic.file, markers.dic.destination)) {
        file.copy(markers.dic.file, markers.dic.destination, overwrite = TRUE)
      }
    } else {
      pop.dictionary <- markers.dictionary <- NULL
    }
  }
  pop.dic.file <- markers.dic.file <- bayescan.sub <- NULL
  # command --------------------------------------------------------------------
  command.arguments <- paste(new.data, output.folder, all.trace, parallel.core, n, thin, nbp, pilot, burn, pr.odds)

  os <- Sys.info()[['sysname']]

  if (os == "Windows") {
    system(command = paste0(bayescan.path, command.arguments))
  } else {
    system2(
      command = bayescan.path,
      args = command.arguments,
      stderr = log.file, stdout = log.file
    )
  }


  # Importing BayeScan file  ---------------------------------------------------
  message("Importing BayeScan results")
  res$bayescan <- suppressWarnings(readr::read_table2(
    file = list.files(path = path.folder.subsample, pattern = "_fst.txt", full.names = TRUE),
    skip = 1,
    col_names = c("BAYESCAN_MARKERS", "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST"),
    col_types = c("iddddd"))) %>%
    dplyr::mutate(
      SELECTION = classify_bayescan_selection(
        alpha = ALPHA,
        q.value = Q_VALUE,
        fdr = fdr
      ),
      Q_VALUE = dplyr::if_else(Q_VALUE <= 0.0001, 0.0001, Q_VALUE),
      Q_VALUE = round(Q_VALUE, 4),
      POST_PROB = round(POST_PROB, 4),
      LOG10_PO = round(LOG10_PO, 4),
      ALPHA = round(ALPHA, 4),
      FST = round(FST, 6),
      LOG10_Q = log10(Q_VALUE)
    )

  if (!is.null(markers.dictionary)) {
    res$bayescan <- dplyr::right_join(markers.dictionary, res$bayescan, by = "BAYESCAN_MARKERS")
  } else {
    res$bayescan <- dplyr::mutate(res$bayescan, MARKERS = BAYESCAN_MARKERS)
  }

  res$bayescan <- res$bayescan %>%
    dplyr::mutate(# Grouping des groupes LOG10_PO & #Quantile of FST
      PO_GROUP = factor(
        dplyr::if_else(LOG10_PO > 2, "decisive",
                       dplyr::if_else(LOG10_PO > 1.5, "very strong",
                                      dplyr::if_else(LOG10_PO > 1, "strong",
                                                     dplyr::if_else(LOG10_PO > 0.5, "substantial", "no evidence")))),
        levels = c("no evidence","substantial","strong","very strong","decisive"), ordered = TRUE)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      FST_GROUP = dplyr::ntile(FST, 5),
      FST_GROUP = dplyr::if_else(FST_GROUP == 1, "0-20%",
                                 dplyr::if_else(FST_GROUP == 2,  "20-40%",
                                                dplyr::if_else(FST_GROUP == 3, "40-60%",
                                                               dplyr::if_else(FST_GROUP == 4, "60-80%", "80-100%"))))
    ) %>%
    dplyr::arrange(FST)

  if (!is.null(subsample)) {
    res$bayescan <- dplyr::mutate(res$bayescan, ITERATIONS = rep(subsample.id, n()))
  }
  # Accuracy within LOCUS ------------------------------------------------------
  # special concern when > 1 SNP / LOCUS...
  radr.markers <- dplyr::distinct(res$bayescan, MARKERS)
  radr.markers <- dplyr::filter(radr.markers, !is.na(MARKERS))

  radr.markers <- unique(stringi::stri_detect_fixed(
    str = radr.markers$MARKERS, pattern = "__"))
  if (radr.markers) {
    message("Detected SNP and LOCUS information in markers")
    res$bayescan <- res$bayescan %>%
      tidyr::separate(
        data = .,
        col = MARKERS,
        into = c("CHROM", "LOCUS", "POS"),
        sep = "__",
        remove = FALSE,
        extra = "warn"
      )

    whitelist.multiple.snp <- dplyr::distinct(res$bayescan, LOCUS, POS) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n > 1) %>%
      dplyr::select(LOCUS)
    markers.more.snp <- nrow(whitelist.multiple.snp)
    n.markers <- dplyr::n_distinct(res$bayescan$MARKERS)
    if (markers.more.snp > 0) {
      message("Detected markers > 1 SNP per LOCUS...")
      message("    total number of markers: ", n.markers)
      message("    markers with >1 SNPs/LOCUS: ", markers.more.snp, " (", round(markers.more.snp/n.markers, 2), ")")
      message("\nCalculating accuracy within LOCUS...")

      locus.accuracy <- dplyr::left_join(whitelist.multiple.snp, res$bayescan, by = "LOCUS") %>%
        dplyr::select(-BAYESCAN_MARKERS, -MARKERS) %>%
        dplyr::group_by(LOCUS, SELECTION) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::mutate(
          SNP_NUMBER = sum(n),
          ACCURACY = dplyr::if_else(n == SNP_NUMBER, "accurate", "not accurate"))

      accurate.locus <- locus.accuracy %>%
        dplyr::filter(ACCURACY == "accurate") %>%
        dplyr::select(LOCUS, SELECTION, SNP_NUMBER)
      n.accurate.locus <- nrow(accurate.locus)
      n.not.accurate.locus <- markers.more.snp - n.accurate.locus
      message("Number of locus accurate: ", n.accurate.locus, " (", round(n.accurate.locus/markers.more.snp, 2), ")")
      message("Number of locus NOT accurate: ", n.not.accurate.locus, " (", round(n.not.accurate.locus/markers.more.snp, 2), ")")
      res$accurate.locus.summary <- accurate.locus %>%
        dplyr::group_by(SELECTION) %>%
        dplyr::tally(.)

      res$whitelist.accurate.locus <- dplyr::distinct(accurate.locus, LOCUS) %>%
        dplyr::arrange(LOCUS)
      readr::write_tsv(
        x = res$whitelist.accurate.locus,
        file = file.path(path.folder.subsample, "whitelist.accurate.locus.tsv"))

      res$blacklist.not.accurate.locus <- locus.accuracy %>%
        dplyr::filter(ACCURACY == "not accurate") %>%
        dplyr::distinct(LOCUS) %>%
        dplyr::arrange(LOCUS)
      readr::write_tsv(
        x = res$blacklist.not.accurate.locus,
        file = file.path(path.folder.subsample, "blacklist.not.accurate.locus.tsv"))


      # correlation between number of snps and accuracy... ?
      res$accuracy.snp.number <- locus.accuracy %>%
        dplyr::distinct(LOCUS, SNP_NUMBER, ACCURACY) %>%
        dplyr::group_by(SNP_NUMBER, ACCURACY) %>%
        dplyr::tally(.)
      readr::write_tsv(
        x = res$accuracy.snp.number,
        file = file.path(path.folder.subsample, "accuracy.snp.number.tsv"))

      res$accuracy.snp.number.plot <- ggplot2::ggplot(res$accuracy.snp.number, ggplot2::aes(y = n, x = SNP_NUMBER, fill = ACCURACY)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(y = "Number of locus") +
        ggplot2::labs(x = "Number of SNPs per locus") +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
          legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
          legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
          strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
        )

      ggplot2::ggsave(
        filename = file.path(path.folder.subsample, "accuracy.snp.number.plot.pdf"),
        plot = res$accuracy.snp.number.plot,
        width = 20, height = 15,
        dpi = 600, units = "cm", useDingbats = FALSE)

      res$not.accurate.summary <- locus.accuracy %>%
        dplyr::filter(ACCURACY == "not accurate") %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(
          SELECTION_TYPE_ON_LOCUS = stringi::stri_join(SELECTION, collapse = " <-> ")
        ) %>%
        dplyr::group_by(SELECTION_TYPE_ON_LOCUS) %>%
        dplyr::tally(.) %>%
        dplyr::rename(LOCUS_NUMBER = n) %>%
        dplyr::mutate(
          PROP = round(LOCUS_NUMBER/sum(LOCUS_NUMBER), 4),
          PROP_TOTAL_MARKERS = round(LOCUS_NUMBER/markers.more.snp, 4))
      readr::write_tsv(
        x = res$not.accurate.summary,
        file = file.path(path.folder.subsample, "not.accurate.summary.tsv"))
    }

  }
  radr.markers <- NULL

  # selection ------------------------------------------------------------------
  res$whitelist.markers.positive.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "diversifying" & PO_GROUP != "no evidence") %>%
    # dplyr::filter(SELECTION == "diversifying") %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange (MARKERS)

  if (!is.null(subsample)) {
    res$whitelist.markers.positive.selection <- dplyr::mutate(res$whitelist.markers.positive.selection, ITERATIONS = rep(subsample.id, n()))
  }

  readr::write_tsv(
    x = res$whitelist.markers.positive.selection,
    file = file.path(path.folder.subsample, "whitelist.markers.positive.selection.tsv"))

  res$whitelist.markers.neutral.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "neutral") %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange (MARKERS)

  if (!is.null(subsample)) {
    res$whitelist.markers.neutral.selection <- dplyr::mutate(res$whitelist.markers.neutral.selection, ITERATIONS = rep(subsample.id, n()))
  }
  readr::write_tsv(
    x = res$whitelist.markers.neutral.selection,
    file = file.path(path.folder.subsample, "whitelist.markers.neutral.selection.tsv"))

  # neutral and positive
  res$whitelist.markers.neutral.positive.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "neutral" | (SELECTION == "diversifying" & PO_GROUP != "no evidence")) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange (MARKERS)
  if (!is.null(subsample)) {
    res$whitelist.markers.neutral.positive.selection <- dplyr::mutate(res$whitelist.markers.neutral.positive.selection, ITERATIONS = rep(subsample.id, n()))
  }
  readr::write_tsv(
    x = res$whitelist.markers.neutral.positive.selection,
    file = file.path(path.folder.subsample, "whitelist.markers.neutral.positive.selection.tsv"))

  res$blacklist.markers.balancing.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "balancing") %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange (MARKERS)
  if (!is.null(subsample)) {
    res$blacklist.markers.balancing.selection <- dplyr::mutate(res$blacklist.markers.balancing.selection, ITERATIONS = rep(subsample.id, n()))
  }
  readr::write_tsv(
    x = res$blacklist.markers.balancing.selection,
    file = file.path(path.folder.subsample, "blacklist.markers.balancing.selection.tsv"))

  # Get the numbers of LOCI under various evolutionary forces
  # Get the numbers for markers under directional selection
  selection <- dplyr::group_by(res$bayescan, SELECTION, PO_GROUP) %>%
    dplyr::tally(.) %>%
    dplyr::rename(MARKERS = n)
  if (!is.null(subsample)) {
    selection <- dplyr::mutate(selection, ITERATIONS = rep(subsample.id, n()))
  }
  readr::write_tsv(
    x = selection,
    file = file.path(path.folder.subsample, "selection.summary.tsv"))



  # Generating plot ------------------------------------------------------------
  message("Generating plot")
  res$bayescan.plot <- plot_bayescan(res$bayescan)

  if (!is.null(subsample)) {
    temp.name <-  stringi::stri_join("bayescan_plot_", subsample.id, ".pdf")
    ggplot2::ggsave(
      filename = file.path(path.folder.subsample, temp.name),
      plot = res$bayescan.plot,
      width = 30, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE)
  } else {
    ggplot2::ggsave(
      filename = file.path(path.folder.subsample, "bayescan_plot.pdf"),
      plot = res$bayescan.plot,
      width = 30, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE)
  }

  # Saving bayescan data frame--------------------------------------------------
  if (!is.null(subsample)) {
    temp.name <-  stringi::stri_join("bayescan_", subsample.id, ".tsv")
    readr::write_tsv(
      x = res$bayescan,
      file = file.path(path.folder.subsample, temp.name))
  } else {
    readr::write_tsv(
      x = res$bayescan,
      file = file.path(path.folder.subsample, "bayescan.tsv"))
  }
  # Update results list --------------------------------------------------------
  res$selection.summary <- selection
  if (!is.null(markers.dictionary)) res$markers.dictionary <- markers.dictionary
  if (!is.null(pop.dictionary)) res$pop.dictionary <- pop.dictionary
  return(res)
} #End bayescan_one


# BayeScan plot function -------------------------------------------------------
#' @title plot_bayescan
#' @description plot_bayescan
#' @rdname plot_bayescan
#' @export
#' @keywords internal

plot_bayescan <- function(data){
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = LOG10_Q, y = FST)) +
    ggplot2::geom_point(ggplot2::aes(colour = PO_GROUP,shape = FST_GROUP)) +
    ggplot2::scale_shape_manual(name = "FST quantile group",values = c(5,2,3,4,1)) +
    ggplot2::scale_colour_manual(name = "Model choice",values = c("darkred","yellow","orange","green","forestgreen")) +
    ggplot2::labs(x = "Log10(Q_VALUE)") +
    ggplot2::labs(y = "Fst") +
    ggplot2::geom_vline(xintercept = c(log10(0.05)),color = "black") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16, family = "Helvetica",face = "bold"),
      legend.title = ggplot2::element_text(size = 16,family = "Helvetica",face = "bold"),
      legend.text = ggplot2::element_text(size = 16,family = "Helvetica",face = "bold"),
      legend.position = "right")
  return(plot)
}#End plot_bayescan
