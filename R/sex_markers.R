
#' @name sexy_markers
#' @title Identify sex-linked markers and reassign genetic sex
#'
#' @description Identifies markers putatively located on heterogametic or
#' homogametic sex chromosomes and can reassign genetic sex using candidate
#' Y- or W-linked markers. The function works best in: DArT silico (counts) >
#' DArT counts or RADseq with allele read depth > DArT silico (genotypes) >
#' RADseq (genotypes) and DArT (1-row, 2-rows genotypes).

#' @param data A genomic input supported by
#'   \code{\link[genometranslator]{read_genome}} that resolves to a GDS object.
#'   An open GDS object is preferred for repeated analyses. VCF, DArT, PLINK,
#'   and other supported inputs can be translated automatically. For
#'   format-specific import control, call the corresponding
#'   \code{genometranslator::read_*()} function first.

#' @param silicodata (optional, file) A silico (dominant marker) DArT file \code{.csv or .tsv}.\cr
#' This can be count or genotyped data. Note that both \code{data} and
#' \code{silicodata} can be used at the same time.\cr
#' Default: \code{silicodata = NULL}.

#' @param strata (file) A tab delimited file with a minimum of
#' 2 columns (\code{INDIVIDUALS, STRATA}) for VCF files and 3 columns for DArT files
#' (\code{TARGET_ID, INDIVIDUALS, STRATA}).
#' \itemize{
#' \item \code{TARGET_ID:} it's the column header of the DArT file.
#' \item \code{STRATA:} is the grouping column, here the sex info.
#' 3 values: \code{M} for male, \code{F} for female and \code{U} for unknown.
#' Anything else is converted to \code{U}.
#' \item \code{INDIVIDUALS:} Is how you want your samples to be named.
#' }
#' Default: \code{strata = NULL}, the function will look for sex info in the
#' tidy data or GDS data (individuals.meta section).
#' You can easily build the strata file by starting with the output of these
#' functions: \code{\link[genometranslator]{extract_dart_target_id}} and
#' \code{\link[genometranslator]{extract_individuals_vcf}}


#' @param boost.analysis (optional, logical) This method uses machine learning
#' approaches to find sex markers and re-assign samples in sex group.\cr
#' The approach is currently under construction.
#' Default: \code{boost.analysis = FALSE}.

#' @param coverage.thresholds (optional, integer) The minimum coverage required
#' to call a marker absent. For silico genotype data this must be < 1.\cr
#' Default: \code{coverage.thresholds = 1}.
#'
#'
#' @param filters (optional, logical) When \code{filters = TRUE}, the data will
#' be filtered for monomorphic loci, missingness of individuals and heterozygosity of individuals.\cr
#' CAUTION: we advice to use these filter, since not filtering or filtering too
#' stringently will results in false positive or false negative detections.\cr
#' Default: \code{filters = TRUE}.
#'
#' @param interactive.filter (optional, logical) When \code{interactive.filter = TRUE}
#' the function will ask for your input to define thresholds. If \code{interactive.filter = FALSE}
#' the function expects additional arguments (see Advance mode).\cr
#' Default: \code{interactive.filter = TRUE}.
#'
#' @param folder.name (optional,character) Name of the folder to store the results.
#' Default: \code{folder.name = NULL}. The name sexy_markers_datetime will be generated.
#'
#' @param verbose Logical. Display progress messages.
#' Default: \code{verbose = TRUE}.
#'
#' @inheritParams radr_common_arguments
#'
#'
#' @details
#' This function takes DArT and RAD-type data to find markers that have a specific
#' pattern that is linked to sex.\cr
#' The function hypothesizes the presence of sex-chromosomes in your
#' species/population. The tests are designed to identify markers that are
#' located on putative heterogametic (Y or W) or homogametic (X or Z) chromosomes.\cr
#' \emph{Note:} Violating Assumptions or Prerequisites (see below) can lead to
#' false positive or the absence of detection of sex-linked markers.



#' @section Assumptions:
#' \enumerate{
#' \item \strong{Genetic Sex Determination System} over a e.g. environmental-sex-determination system.
#' \item \strong{Genome coverage:} restriction sites randomly spread throughout the whole genome.
#' \item \strong{Mutations:} Processes such as sex-specific mutations
#' in the restriction sites could lead to false positive results.
#' \item \strong{Deletions/duplications:} Processes such as
#' sex-specific deletions or duplications could lead to false positive results.
#' \item \strong{Homology:} The existence of homologous sequences between the
#' homogametic and heterogametic chromosomes could lead to false negative
#' results.
#' \item \strong{Absence of population signal}
#' }

#' @section Prerequisites:
#' \enumerate{
#' \item \strong{Sample size:} Ideally, the data must have enough individuals (n > 100).
#' \item \strong{Batch effect:} Sex should be randomized on lanes/chips during sequencing.
#' \item \strong{Sex ratio:} Dataset with equal ratio work best.
#' \item \strong{Genotyping rate}: for DArT data, if the minimum call rate is
#' greater than 0.5, ask DArT to lower its filtering threshold.
#' RADseq data, lower markers missingness thresholds during filtering
#' (e.g. stacks \code{r} and \code{p}).
#' \item \strong{Identity-by-Missingness:} Absence of artifactual pattern
#' of missingness (\href{https://github.com/thierrygosselin/grur}{see missing visualization})
#' \item \strong{Low genotyping error rate:} see \code{\link{detect_het_outliers}}
#' and \href{https://github.com/eriqande/whoa}{whoa}.
#' \item \strong{Low heterozygosity miscall rate:} see \code{\link{detect_het_outliers}} and
#' \href{https://github.com/eriqande/whoa}{whoa}.
#' \item \strong{Absence of pattern of heterozygosity driven by missingness:}
#' see \code{\link[radr]{detect_mixed_genomes}}.
#' \item \strong{Absence of paralogous sequences in the data}.
#' }


#' @section Sex methods:
#' \strong{Heterogametic sex-markers:}
#' \itemize{
#' \item \strong{Presence/Absence method}: To identify markers on Y or W chromosomes,
#' we look at the presence or
#' absence of a marker between females and males. More specifically, if a marker
#' is always present in males but never in females, they are putatively located
#' in the Y-chromosome; and vice versa for the W-linked markers.
#' }
#' \strong{Homogametic sex-markers:}
#' We have two different methods to identify markers on X or Z chromosomes:
#' \itemize{
#' \item \strong{Heterozygosity method:} By looking at the heterozygosity of a
#' marker between sexes, we can
#' identify markers that are always homozygous in one sex (e.g. males for an
#' XY system), while exhibiting an intermediate range of heterozygosity in the
#' other sex (0.1 - 0.5).
#' \item \strong{Coverage method:} If the data includes count information, this
#' function will look for
#' markers that have double the number of counts for either of the sexes.
#' For example if an XY/XX system is present, females are expected to have
#' double the number of counts for markers on the X chromosome.
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \itemize{
#' \item \code{species}: To give your figures some meanings.
#' Default \code{species = NULL}.
#' \item \code{population}: To give your figures some meanings.
#' Default \code{species = NULL}.
#' \item \code{tau}: The quantile used in regression to distinguish homogametic markers
#' with the \strong{heterozygosity method}. See \code{\link[quantreg]{rq}}.\cr
#' Default \code{tau = 0.03}.

#' \item \code{mis.threshold.data}: Threshold to filter the SNP data on missingness.
#' Only if \code{interactive.filter = FALSE}.
#' \item \code{mis.threshold.silicodata}: Threshold to filter the silico data on
#' missingness. No default. Only if \code{interactive.filter = FALSE}.
#' \item \code{threshold.y.markers}: Threshold to select heterogametic sex-linked
#' markers from the SNP data with the \strong{presence/absence method}.No default.
#' Only if \code{interactive.filter = FALSE}.
#' \item \code{threshold.x.markers.qr}: Threshold to select homogametic sex-linked
#' markers from the SNP data with the \strong{heterozygosity method}. No default.
#' Only if \code{interactive.filter = FALSE}.
#' \item \code{zoom.data}: Threshold to subset the F/M ratio of mean SNP coverage.
#' Used to improve the histogram resolution to select a better \code{threshold.x.markers.RD}
#' threshold. No default. Only if \code{interactive.filter = FALSE}.
#' \item \code{threshold.x.markers.RD}: Threshold to select homogametic sex-linked
#' markers from the SNP data with the \strong{coverage method}.No default.
#' Only if \code{interactive.filter = FALSE}.
#' \item \code{zoom.silicodata}: Threshold to subset the F/M ratio of mean silico coverage.
#' Used to improve the histogram resolution to select a better \code{threshold.x.markers.RD.silico}
#' threshold. No default. Only if \code{interactive.filter = FALSE}.
#' \item \code{threshold.x.markers.RD.silico}: Threshold to select homogametic sex-linked
#' markers from the silico data with the \strong{coverage method}. No default.
#' Only if \code{interactive.filter = FALSE}.
#' \item \code{sex.id.input}: (integer) \code{sex.id.input = c(1, 2 or 3)}
#' to recalculate the sex based on (1) 'visual', (2) 'genetic SNP' or (3) 'genetic SILICO' sexID.
#' No default. Only if \code{interactive.filter = FALSE}.
#' \item \code{het.qr.input}: (integer) \code{het.qr.input = c(1 or 2)}
#' to plot the heterozygosity residual plot for (1) X-linked markers (heterozygous for females),
#' or (2) Z-linked markers (heterozygous for males). No default.
#' Only if \code{interactive.filter = FALSE}.
#' }
#'
#'
#' @section Life cycle:
#'
#' Machine Learning approaches (Random Forest and Extreme Gradient Boosting Trees)
#' are currently been tested. They usually show a lower discovery rate but tend to
#' perform better with new samples.
#'
#'

#' @seealso Eric Anderson's \href{https://github.com/eriqande/whoa}{whoa} package.

#' @export
#' @rdname sexy_markers

#' @return The created object contains:
#' \enumerate{
#' \item A list with (1) the summarised SNP data per sex and
#' (2) the summarised silico data per sex. This should allow you to re-create the various plots.
#' \item A vector with the names of the sex-linked marker. One vector for each sex method.
#' \item A dataframe with a summary of the sex-linked markers and their sequence (if available).
#' }

#' @examples
#' \dontrun{

#' # The minimum
#' sex.markers <- radr::sexy_markers(
#'     data = "shark.dart.data.csv",
#'     strata = "shark.strata.tsv")

#' # This will use the default: interactive version and a list is created and to view the sex markers
#' }


#' @author Floriaan Devloo-Delva \email{Floriaan.Devloo-Delva@@csiro.au} and with help from
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}

sexy_markers <- function(data,
                         silicodata = NULL,
                         strata = NULL,
                         boost.analysis = FALSE,
                         coverage.thresholds = 1,
                         filters = TRUE,
                         interactive.filter = TRUE,
                         folder.name = NULL,
                         parallel.core = parallel::detectCores() - 1,
                         verbose = TRUE,
                         ...
) {

  # Test
  # silicodata <- "../1.Data/G.galeus/SchoolShark_silico_counts.csv"
  # boost.analysis = FALSE
  # coverage.thresholds = 1 ##NEEDS TO BE SET TO 1 if genotype data
  # filters = TRUE
  # interactive.filter = TRUE
  # parallel.core = 1
  # species = "shark"
  # population = "Australia"
  # tau = 0.03
  # threshold.y.markers = NULL
  # threshold.y.silico.markers = NULL
  # sex.id.input = 2
  # het.qr.input = 1
  # threshold.x.markers.qr = NULL
  # data = "../1.Data/G.galeus/SchoolShark_SNP_counts.csv"
  # strata = "../1.Data/G.galeus/SchoolShark_strata.tsv"

  # parallel.core = parallel::detectCores() - 1

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "sexy_markers",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  res <- list()

  # Validate formal arguments --------------------------------------------------
  if (missing(data)) rlang::abort("Argument `data` is required.")
  logical.args <- c("boost.analysis", "filters", "interactive.filter", "verbose")
  invalid.logical <- logical.args[!vapply(
    logical.args,
    function(x) is.logical(get(x, inherits = FALSE)) &&
      length(get(x, inherits = FALSE)) == 1L &&
      !is.na(get(x, inherits = FALSE)),
    logical(1)
  )]
  if (length(invalid.logical) > 0L) {
    rlang::abort(paste0(
      "Arguments must be single TRUE/FALSE values: ",
      paste(invalid.logical, collapse = ", ")
    ))
  }
  if (!is.numeric(coverage.thresholds) || length(coverage.thresholds) != 1L ||
      is.na(coverage.thresholds) || coverage.thresholds < 0) {
    rlang::abort("`coverage.thresholds` must be one non-negative number.")
  }

  if (boost.analysis) {
    rlang::warn("`boost.analysis` is under development and is not run yet.")
  }

  # required packages ----------------------------------------------------------
  tgbase::check_package(package = "quantreg")

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c(
      "species",
      "population",
      "tau",
      "threshold.y.markers",
      "threshold.y.silico.markers",
      "sex.id.input",
      "threshold.x.markers.qr",
      "threshold.x.markers.RD",
      "threshold.x.markers.RD.silico",
      "mis.threshold.data",
      "mis.threshold.silicodata",
      "zoom.data",
      "zoom.silicodata",
      "sex.id.input",
      "het.qr.input"
    ),
    verbose = FALSE
  )

  if (is.null(tau)) tau <- 0.03

  # Folders---------------------------------------------------------------------
  if (is.null(folder.name)) folder.name <- "sexy_markers"
  wd <- path.folder <- tgbase::generate_folder(
    folder = folder.name,
    path.folder = NULL,
    internal = FALSE,
    prefix.int = FALSE,
    file.date = file.date,
    verbose = verbose
  )

  # write the dots file
  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_sexy_markers_args",
    date = TRUE,
    internal = FALSE,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )


  # Import primary genomic data ------------------------------------------------
  count.data <- FALSE
  input.type <- genometranslator::detect_genomic_format(data)
  if (!identical(input.type, "SeqVarGDSClass")) {
    data <- genometranslator::read_genome(
      data = data,
      strata = strata,
      parallel.core = parallel.core,
      verbose = FALSE
    )
  }

  data.type <- genometranslator::detect_genomic_format(data)
  if (!identical(data.type, "SeqVarGDSClass")) {
    rlang::abort(paste0(
      "`data` must resolve to an open GDS object; `read_genome()` returned format `",
      data.type, "`. Use a format-specific genometranslator reader to create GDS input."
    ))
  }

  # Detect source --------------------------------------------------------------
  data.source <- genometranslator::extract_data_source(gds = data)
  if ("counts" %in% data.source) count.data <- TRUE

  # Filters----------------------------------------------------------------------
  if (filters) {
    data <- radr::filter_monomorphic(
      data = data,
      parallel.core = parallel.core,
      verbose = FALSE,
      internal = FALSE,
      path.folder = path.folder
    ) |>
      radr::filter_individuals(
        interactive.filter = interactive.filter,
        filter.individuals.missing = "outliers",
        filter.individuals.heterozygosity = "outliers",
        filter.individuals.coverage.total = "outliers",
        dp = count.data,
        internal = FALSE,
        path.folder = path.folder,
        parallel.core = parallel.core,
        verbose = FALSE
      ) |>
      radr::filter_ld(
        interactive.filter = interactive.filter,
        filter.short.ld = "mac",
        long.ld.missing = FALSE,
        parallel.core = parallel.core,
        verbose = FALSE,
        internal = FALSE,
        path.folder = path.folder
      )
  }


  # clean the strata -----------------------------------------------------------
  strata <- genometranslator::extract_individuals_metadata(
    gds = data,
    ind.field.select = c("TARGET_ID", "INDIVIDUALS", "STRATA"),
    whitelist = TRUE
  )

  # clean strata----------------------------------------------------------------
  strata %<>%
    dplyr::mutate(
      STRATA = stringi::stri_trans_toupper(str = STRATA),
      STRATA = stringi::stri_sub(str = STRATA, from = 1, to = 1),
      STRATA = replace(x = STRATA, !STRATA %in% c("F", "M"), "U")
    )

  #checks
  strata.groups <- unique(sort(strata$STRATA))
  if (length(strata.groups) > 3 || length(strata.groups) < 2) {
    rlang::abort("The strata requires a minimum of 2 groups: F and M")
  }
  if (length(strata.groups) == 2 &&
      all(c("F", "M") %in% strata.groups) == FALSE) {
    rlang::abort("The strata requires a minimum of 2 groups: F and M")
  }

  # Generate new strata and write to disk
  readr::write_tsv(
    x = strata,
    file = file.path(path.folder, "sexy_markers_strata-original_cleaned.tsv")
  )


  # check Sex Ratio ------------------------------------------------------------
  sex.ratio <- dplyr::filter(strata, STRATA != "U") %>%
    dplyr::count(STRATA, name = "SEX_RATIO")

  message("\n\nSex-ratio (F/M): ", round((sex.ratio$SEX_RATIO[sex.ratio$STRATA == "F"] / sex.ratio$SEX_RATIO[sex.ratio$STRATA == "M"]), 2))


  gds.bk <- data
  data <- genometranslator::extract_genotypes_metadata(
    gds = data,
    genotypes.meta.select = c(
      "MARKERS", "INDIVIDUALS", "ALT_DOSAGE", "READ_DEPTH"
    ),
    whitelist = TRUE
  ) |>
    genometranslator::join_strata(
      strata = strata,
      verbose = FALSE
    )

  if (!rlang::has_name(data, "ALT_DOSAGE")) {
    rlang::abort(
      "ALT_DOSAGE is unavailable in the GDS genotype metadata; recreate the GDS with current genometranslator readers."
    )
  }



  # SILICO files ----------------------------------------------------------------

  if (!is.null(silicodata)) {
    if (is.character(silicodata)) {
      data.type <- genometranslator::detect_genomic_format(silicodata)
      silicodata <- genometranslator::read_dart(
        data = silicodata,
        strata = strata,
        tidy.dart = FALSE,
        parallel.core = parallel.core,
        path.folder = path.folder,
        internal = TRUE,
        verbose = verbose
      ) %>%
        dplyr::rename(., MARKERS = CLONE_ID)
      if (max(silicodata$VALUE, na.rm = TRUE) > 1) {
        data.source <- c("counts", data.type, data.source)
      } else{
        data.source <- c("genotype", data.type, data.source)
      }
    } else if (is.data.frame(silicodata)) {
      data.type <- "silico.dart"
      if (max(silicodata$VALUE, na.rm = TRUE) > 1) {
        data.source <- c("counts", data.type, data.source)
      } else{
        data.source <- c("genotype", data.type, data.source)
      }
    } else {
      rlang::abort(paste0(
        "`silicodata` must be a DArT filepath or a tidy data frame returned ",
        "by `genometranslator::read_dart()`."
      ))
    }
  }

  ### add warning about not using count data
  if (!(rlang::has_name(data, "READ_DEPTH"))) {
    if (verbose) {
      message(
        "Read depth is unavailable; coverage-based sex-marker methods will be skipped."
      )
    }
  }

  # Validate non-interactive inputs -------------------------------------------
  if (!interactive.filter) {
    required <- c("mis.threshold.data", "sex.id.input", "het.qr.input")
    if (!is.null(silicodata)) required <- c(required, "mis.threshold.silicodata")
    if (rlang::has_name(data, "READ_DEPTH")) required <- c(required, "zoom.data")
    if (all(c("silico.dart", "counts") %in% data.source)) {
      required <- c(required, "zoom.silicodata")
    }

    missing.required <- required[vapply(
      required,
      function(x) !exists(x, inherits = FALSE) || is.null(get(x, inherits = FALSE)),
      logical(1)
    )]
    if (length(missing.required) > 0L) {
      rlang::abort(paste0(
        "With `interactive.filter = FALSE`, supply these arguments in `...`: ",
        paste(missing.required, collapse = ", "), "."
      ))
    }

    missing.thresholds <- c(mis.threshold.data, mis.threshold.silicodata)
    missing.thresholds <- missing.thresholds[!vapply(
      missing.thresholds,
      is.null,
      logical(1)
    )]
    if (length(missing.thresholds) > 0L &&
        any(!is.finite(missing.thresholds) |
            missing.thresholds < 0 | missing.thresholds > 1)) {
      rlang::abort("Missingness thresholds must be between 0 and 1.")
    }
    if (!sex.id.input %in% 1:3) {
      rlang::abort("`sex.id.input` must be 1, 2, or 3.")
    }
    if (!het.qr.input %in% 1:2) {
      rlang::abort("`het.qr.input` must be 1 or 2.")
    }
  }


  # START Sex markers-----------------------------------------------------------------
  # 1. presence/absence per strata and markers: yw
  # 2. heterozygosity per strata and markers: xz
  # 3. read depth per strata and markers markers: xz

  {
    cat("\n\n################################################################################\n")
    cat("######################## Start finding sex-linked markers ######################\n")
    cat("################################################################################\n")
  }


  # Summarise DArT and VCFs---------------------------------------
  if (is.null(tau)) tau <- 0.03

  res$sum <- summarize_sex(
    data = data,
    silicodata = silicodata,
    coverage.thresholds = coverage.thresholds,
    data.source = data.source,
    tau = tau
  )
  data.sum <- res$sum$data.sum
  silico.sum <- res$sum$silico.sum

  # FIGURES -------------------------------------------------------------------

  # Here we start to generate the figures for each sex determination mechanism
  # The function for figures was moved outside the function below
  # you need to load the function first so that it's accessible

  SexID <- "visually"

  #### P/A #####

  # TODO perhaps for count data, set coverage threshold, based on plot?

  ###* For data.sum ####
  ### SCATTER ----
  plot.filename <- stringi::stri_join("1A.sexy_markers_PA_scatter_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
    stringi::stri_join(".pdf")

  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "PRESENCE_ABSENCE_F",
    y = "PRESENCE_ABSENCE_M",
    scat = TRUE,
    qreg = FALSE,
    tuckey = FALSE,
    x.title = "Proportion of females",
    y.title = "Proportion of males",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Absence of each SNP marker between females and males",
    plot.filename = plot.filename,
    path = path.folder
  )
  print(scat.fig)

  message("Figure written: ", plot.filename)

  ### TUCKEY ----
  plot.filename <- stringi::stri_join("1B.sexy_markers_PA_tuckey_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
    stringi::stri_join(".pdf")

  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "MEAN",
    y = "DIFF",
    scat = TRUE,
    qreg = FALSE,
    tuckey = TRUE,
    x.title = "Mean of males and females",
    y.title = "Difference between males and females (<- W/Y ->)",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Tukey mean-difference plot of each SNP marker between females \nand males",
    plot.filename = plot.filename,
    path = path.folder
  )
  print(scat.fig)
  message("Figure written: ", plot.filename)



  # Interacive selection of threshold
  if (interactive.filter) {
    filter.y.markers <- tgbase::question(x = "P/A method of SNPs: \nLook at the figures: Do you want to select Y/W-linked markers (y/n): ", answer.opt = c("y", "n"))

    if (filter.y.markers == "y") {
      threshold.y.markers <-
        tgbase::question(x = "Choose the threshold for Y/W-linked markers, note that the Y-axis is inverted (-1 to 1): ",
                                    minmax = c(-1.5, 1.5))
    } else {
      threshold.y.markers <- NULL
    }
  }

  # Select and print sex markers
  if (!is.null(threshold.y.markers)) {
    if (threshold.y.markers < 0) {
      y.markers <- dplyr::filter(data.sum, DIFF < threshold.y.markers)$MARKERS
    } else if (threshold.y.markers > 0) {
      y.markers <- dplyr::filter(data.sum, DIFF > threshold.y.markers)$MARKERS
    }
    res$heterogametic.markers <- y.markers
  } else{
    y.markers <- NULL
    res$heterogametic.markers <- NULL
  }

  if (length(y.markers) == 0) {
    res$heterogametic.markers <- NULL
  }


  ###* For silico.sum ####
  if ("silico.dart" %in% data.source) {
    ### SCATTER

    plot.filename <- stringi::stri_join("2A.sexy_markers_SILICO_PA_scatter_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "PRESENCE_ABSENCE_F",
      y = "PRESENCE_ABSENCE_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      x.title = "Proportion of females",
      y.title = "Proportion of males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Absence of each SILICO marker between females and males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(scat.fig)
    message("Figure written: ", plot.filename)

    ### TUCKEY
    plot.filename <- stringi::stri_join("2B.sexy_markers_SILICO_PA_tuckey_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
    stringi::stri_join(".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "MEAN",
      y = "DIFF",
      scat = TRUE,
      qreg = FALSE,
      tuckey = TRUE,
      x.title = "Mean of males and females",
      y.title = "Difference between males and females (<- W/Y ->)",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Tukey mean-difference plot of each SILICO marker between females \nand males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(scat.fig)

    message("Figure written: ", plot.filename)


    # Interacive selection of threshold
    if (interactive.filter) {
      filter.y.markers <- tgbase::question(x = "P/A method of SILICOs:\nLook at the figures: Do you want to select Y/W-linked markers (y/n): ", answer.opt = c("y", "n"))

      if (filter.y.markers == "y") {
        threshold.y.silico.markers <-
          tgbase::question(x = "Choose the threshold for Y/W-linked SILICO markers, note that the Y-axis is inverted (-1 to 1): ",
                                      minmax = c(-1.5, 1.5))
      } else {
        threshold.y.silico.markers <- NULL
        y.silico.markers <- NULL
      }
    }

    # Select and print sex markers
    if (!is.null(threshold.y.silico.markers)) {
      if (threshold.y.silico.markers < 0) {
        y.silico.markers <- dplyr::filter(silico.sum, DIFF < threshold.y.silico.markers)$MARKERS
      } else if (threshold.y.silico.markers > 0) {
        y.silico.markers <- dplyr::filter(silico.sum, DIFF > threshold.y.silico.markers)$MARKERS
      }
      res$heterogametic.silico.markers <- y.silico.markers
    } else{
      y.silico.markers <- NULL
      res$heterogametic.silico.markers <- NULL
    }
    # if (length(y.silico.markers) == 0){
    #   res$heterogametic.silico.markers <- NULL
  } else{
    y.silico.markers <- NULL
    res$heterogametic.silico.markers <- NULL
  }



  ####* SET GENETIC SEX STRATA ####
  if (exists("y.markers") &&
      !rlang::is_empty(y.markers)) {
    ####** Y markers ####
    if (threshold.y.markers < 0) {
      if (rlang::has_name(data, "READ_DEPTH")) {
        y.data <-
          dplyr::filter(data, MARKERS %in% y.markers) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          dplyr::summarise(MEAN_RD = mean(READ_DEPTH, na.rm = TRUE)) %>%
          dplyr::ungroup(.) %>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::mutate(VISUAL_SEX = STRATA) %>%
          dplyr::mutate(GENETIC_SEX =
                          dplyr::case_when(
                            is.na(MEAN_RD) |
                              MEAN_RD <= coverage.thresholds ~ "F",
                            !(is.na(MEAN_RD) |
                                MEAN_RD <= coverage.thresholds) ~ "M"
                          )) %>%
          dplyr::mutate(STRATA = GENETIC_SEX) %>%
          dplyr::select(-TARGET_ID)
      }
      else if (rlang::has_name(data, "ALT_DOSAGE")) {
        if (length(y.markers) > 1) {
          y.data <- dplyr::filter(data, MARKERS %in% y.markers) %>%
            dplyr::group_by(INDIVIDUALS) %>%
            dplyr::summarise(MEAN_GT = mean(ALT_DOSAGE),
                             NA_COUNT = sum(is.na(ALT_DOSAGE))) %>%
            dplyr::ungroup(.) %>%
            dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(GENETIC_SEX =
                            dplyr::case_when(
                              # is.na(MEAN_GT) &
                              NA_COUNT >= length(y.markers)/2 ~ "F", #majority rule
                              !(
                                # is.na(MEAN_GT) &
                                NA_COUNT >= length(y.markers)/2) ~ "M"
                              # VISUAL_SEX == "M" & (
                              #   # is.na(MEAN_GT) &
                              #   NA_COUNT >= length(y.markers)/2) ~ "U"
                            )) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(-TARGET_ID)
        } else {
          message("You have only 1 Y-linked marker: 'M' for which this marker is absent are reassigned as 'U' since this marker could be absent due to an error.")
          y.data <- dplyr::filter(data, MARKERS %in% y.markers) %>%
            dplyr::mutate(MEAN_GT = ALT_DOSAGE) %>%
            # dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(GENETIC_SEX =
                            dplyr::case_when(
                              is.na(MEAN_GT) ~ "F",
                              !(is.na(MEAN_GT)) ~ "M",
                              VISUAL_SEX == "M" &
                                is.na(MEAN_GT) ~ "U"
                            )) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(INDIVIDUALS, MEAN_GT, STRATA, VISUAL_SEX, GENETIC_SEX)
        }
      } else {
        message("READ_DEPTH and ALT_DOSAGE not availabe in the data.")
      }
      ####** W markers ####
    } else if (threshold.y.markers > 0) {
      if (rlang::has_name(data, "READ_DEPTH")) {
        y.data <-
          dplyr::filter(data, MARKERS %in% y.markers) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          dplyr::summarise(MEAN_RD = mean(READ_DEPTH, na.rm = TRUE)) %>%
          dplyr::ungroup(.) %>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::mutate(VISUAL_SEX = STRATA) %>%
          dplyr::mutate(GENETIC_SEX =
                          dplyr::case_when(
                            is.na(MEAN_RD) |
                              MEAN_RD <= coverage.thresholds ~ "M",
                            !(is.na(MEAN_RD) |
                                MEAN_RD < coverage.thresholds) ~ "F"
                          )) %>%
          dplyr::mutate(STRATA = GENETIC_SEX) %>%
          dplyr::select(-TARGET_ID)

      } else if (rlang::has_name(data, "ALT_DOSAGE")) {
        if (length(y.markers) > 1) {
          y.data <- dplyr::filter(data, MARKERS %in% y.markers) %>%
            dplyr::group_by(INDIVIDUALS) %>%
            dplyr::summarise(MEAN_GT = mean(ALT_DOSAGE),
                             NA_COUNT = sum(is.na(ALT_DOSAGE))) %>%
            dplyr::ungroup(.) %>%
            dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(GENETIC_SEX =
                            dplyr::case_when(
                              # is.na(MEAN_GT) &
                              NA_COUNT >= length(y.markers)/2  ~ "M", #majority rule
                              !(
                                # is.na(MEAN_GT) &
                                NA_COUNT >= length(y.markers)/2) ~ "F"
                              # VISUAL_SEX == "F" & (
                              #   # is.na(MEAN_GT) &
                              #   NA_COUNT >= length(y.markers)/2) ~ "U"
                            )) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(-TARGET_ID)
        } else {
          message("You have only 1 W-linked marker: 'F' for which this marker is absent are reassigned as 'U'.
This marker could be absent due to an error.")
          y.data <- dplyr::filter(data, MARKERS %in% y.markers) %>%
            dplyr::mutate(MEAN_GT = ALT_DOSAGE) %>%
            # dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(GENETIC_SEX =
                            dplyr::case_when(
                              is.na(MEAN_GT) ~ "M",
                              !(is.na(MEAN_GT)) ~ "F",
                              VISUAL_SEX == "F" &
                                is.na(MEAN_GT) ~ "U"
                            )) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(INDIVIDUALS, MEAN_GT, STRATA, VISUAL_SEX, GENETIC_SEX)
        }
      } else {
        rlang::abort("READ_DEPTH and ALT_DOSAGE not availabe in the data.")
      }
    }
    # print visual and gentic sex table
    SumTable <-
      y.data %>%
      dplyr::count(VISUAL_SEX, GENETIC_SEX) %>%
      dplyr::rename(Visual_Sex = VISUAL_SEX, Genetic_Sex_SNP = GENETIC_SEX)
    print(SumTable)
    readr::write_tsv(
      x = SumTable,
      file = file.path(
        path.folder,
        "sexy_markers_summary_table_visual.VS.geneticSNP_sex.tsv"
      )
    )
  }

  if (exists("y.silico.markers") &&
      !rlang::is_empty(y.silico.markers)) {
    ###* Y SILICO markers ####
    if (threshold.y.silico.markers < 0) {
      if (max(silicodata$VALUE, na.rm = TRUE) > 1) {
        y.silico.data <-
          dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          dplyr::summarise(MEAN_RD = mean(VALUE, na.rm = TRUE)) %>%
          dplyr::ungroup(.) %>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::mutate(VISUAL_SEX = STRATA) %>%
          dplyr::mutate(
            GENETIC_SEX =
              dplyr::case_when(
                MEAN_RD <= coverage.thresholds ~ "F",
                is.na(MEAN_RD) ~ "U",
                MEAN_RD > coverage.thresholds &
                  !is.na(MEAN_RD) ~ "M"
              )
          ) %>%
          dplyr::mutate(STRATA = GENETIC_SEX) %>%
          dplyr::select(-TARGET_ID)
      } else { #GT
        message("If all Y-linked markers are assigned NA (i.e. uncertain genotype), this individuals' sex is assigned 'U'.")
        if (length(y.silico.markers) > 1) {
          y.silico.data <-
            dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
            dplyr::group_by(INDIVIDUALS) %>%
            dplyr::summarise(MEAN_GT = mean(VALUE, na.rm = TRUE)) %>%
            dplyr::ungroup(.) %>%
            dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(
              GENETIC_SEX =
                dplyr::case_when(
                  MEAN_GT < coverage.thresholds / length(y.silico.markers) ~ "F",
                  is.na(MEAN_GT) ~ "U",
                  MEAN_GT >= coverage.thresholds / length(y.silico.markers) &
                    !is.na(MEAN_GT) ~ "M"
                )
            ) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(-TARGET_ID)
        } else {
          y.silico.data <-
            dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
            dplyr::mutate(MEAN_GT = VALUE) %>%
            # dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(
              GENETIC_SEX =
                dplyr::case_when(
                  MEAN_GT < coverage.thresholds ~ "F",
                  is.na(MEAN_GT) ~ "U",
                  MEAN_GT >= coverage.thresholds &
                    !is.na(MEAN_GT) ~ "M"
                )
            ) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(INDIVIDUALS, MEAN_GT, STRATA, VISUAL_SEX, GENETIC_SEX)
        }
      }
    }
    ####** W SILICO markers ####
    else if (threshold.y.silico.markers > 0) {
      if (max(silicodata$VALUE, na.rm = TRUE) > 1) {
        y.silico.data <-
          dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          dplyr::summarise(MEAN_RD = mean(VALUE, na.rm = TRUE)) %>%
          dplyr::ungroup(.) %>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::mutate(VISUAL_SEX = STRATA) %>%
          dplyr::mutate(
            GENETIC_SEX =
              dplyr::case_when(
                MEAN_RD <= coverage.thresholds ~ "M",
                is.na(MEAN_RD) ~ "U",
                MEAN_RD > coverage.thresholds &
                  !is.na(MEAN_RD) ~ "F"
              )
          ) %>%
          dplyr::mutate(STRATA = GENETIC_SEX) %>%
          dplyr::select(-TARGET_ID)
      } else { #GT
        message("If all W-linked markers are assigned NA (i.e. uncertain genotype), this individuals' sex is assigned 'U'.")
        if (length(y.silico.markers) > 1) {
          y.silico.data <-
            dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
            dplyr::group_by(INDIVIDUALS) %>%
            dplyr::summarise(MEAN_GT = mean(VALUE, na.rm = TRUE)) %>%
            dplyr::ungroup(.) %>%
            dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(
              GENETIC_SEX =
                dplyr::case_when(
                  MEAN_GT < coverage.thresholds / length(y.silico.markers) ~ "M",
                  is.na(MEAN_GT) ~ "U",
                  MEAN_GT >= coverage.thresholds / length(y.silico.markers) &
                    !is.na(MEAN_GT) ~ "F"
                )
            ) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(-TARGET_ID)
        } else {
          y.silico.data <-
            dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
            dplyr::mutate(MEAN_GT = VALUE) %>%
            # dplyr::left_join(strata, by = "INDIVIDUALS") %>%
            dplyr::mutate(VISUAL_SEX = STRATA) %>%
            dplyr::mutate(
              GENETIC_SEX =
                dplyr::case_when(
                  MEAN_GT < coverage.thresholds ~ "M",
                  is.na(MEAN_GT) ~ "U",
                  MEAN_GT >= coverage.thresholds &
                    !is.na(MEAN_GT) ~ "F"
                )
            ) %>%
            dplyr::mutate(STRATA = GENETIC_SEX) %>%
            dplyr::select(INDIVIDUALS, MEAN_GT, STRATA, VISUAL_SEX, GENETIC_SEX)
        }
      }
    }
    # print visual and gentic sex table
    SumTable <-
      y.silico.data %>%
      dplyr::count(VISUAL_SEX, GENETIC_SEX) %>%
      dplyr::rename(Visual_Sex = VISUAL_SEX, Genetic_Sex_SILICO = GENETIC_SEX)
    print(SumTable)
    readr::write_tsv(
      x = SumTable,
      file = file.path(
        path.folder,
        "sexy_markers_summary_table_visual.VS.geneticSILICO_sex.tsv"
      )
    )
  }

  if (exists("y.markers") && !rlang::is_empty(y.markers) &&
      exists("y.silico.markers") && !rlang::is_empty(y.silico.markers)) {
    # print genetic SNP and gentic SILICO sex table
    SumTable <- dplyr::left_join(y.data, y.silico.data, by = "INDIVIDUALS") %>%
      dplyr::count(GENETIC_SEX.x, GENETIC_SEX.y) %>%
      dplyr::rename(Genetic_Sex_SNP = GENETIC_SEX.x, Genetic_Sex_SILICO = GENETIC_SEX.y)
    print(SumTable)
    readr::write_tsv(
      x = SumTable,
      file = file.path(
        path.folder,
        "sexy_markers_summary_table_geneticSNP.VS.geneticSILICO_sex.tsv"
      )
    )
  }

  if (exists("y.markers") && !rlang::is_empty(y.markers) ||
      exists("y.silico.markers") && !rlang::is_empty(y.silico.markers)) {
    # Interacive selection which sex info
    if (is.null(sex.id.input)) sex.id.input <- "1"
    if (!rlang::is_empty(y.markers) &&
        !rlang::is_empty(y.silico.markers)) {
      if (interactive.filter) {
        sex.id.input <-
          tgbase::question(
            x = "For further analysis, do you want to continue based on (1) visual, (2) genetic SNP or (3) genetic SILICO sex?\nWe advise (3) for better results",
            answer.opt = c("1", "2", "3"))
        sex.id.input <- as.integer(sex.id.input)
      } else{
        # sex.id.input <- as.integer(3)
        sex.id.input <- as.integer(sex.id.input)
      }
    } else if (!rlang::is_empty(y.markers)) {
      # Interacive selection which sex info
      if (is.null(sex.id.input)) sex.id.input <- "1"
      if (interactive.filter) {
        sex.id.input <-
          tgbase::question(
            x = "For further analysis, do you want to continue based on (1) visual or (2) genetic SNP sex?\nWe advise (2) for better results",
            answer.opt = c("1", "2"))
        sex.id.input <- as.integer(sex.id.input)
      } else {
        # sex.id.input <- as.integer(2)
        sex.id.input <- as.integer(sex.id.input)
      }
    } else if (!rlang::is_empty(y.silico.markers)) {
      # Interacive selection which sex info
      if (is.null(sex.id.input)) sex.id.input <- "1"
      if (interactive.filter) {
        sex.id.input <-
          tgbase::question(x = "For further analysis, do you want to continue based on (1) visual or (3) genetic SILICO sex?\nWe advise (3) for better results",
                                      answer.opt = c("1", "3"))
        sex.id.input <- as.integer(sex.id.input)
      } else {
        # sex.id.input <- as.integer(3)
        sex.id.input <- as.integer(sex.id.input)
      }
    }

    #set new sexID for Het analysis
    if (sex.id.input == 2) {
      readr::write_tsv(
        x = y.data,
        file = file.path(
          path.folder,
          "sexy_markers_strata_new_with_genetic_sex_according_to_SNPdata.tsv"
        )
      )
      message("New strata file with genetic sex is written.")
      ### recalculate data based on new sexID
      data <-
        dplyr::select(data,-STRATA) %>%
        # dplyr::rename(data, VISUAL_STRATA = STRATA) %>%
        dplyr::left_join(y.data, by = "INDIVIDUALS")
      # dplyr::rename(MARKERS = MARKERS.x, ALT_DOSAGE = ALT_DOSAGE.x, TARGET_ID = TARGET_ID.x) %>%
      # dplyr::rename(TARGET_ID = TARGET_ID.x) %>%
      # dplyr::select(-c(TARGET_ID.y))
      # dplyr::mutate(GENETIC_STRATA = STRATA)
      genometranslator::write_genome(
        data = data,
        filename = file.path(wd, "sexy_markers_genetic_SNP_sex_ID.arrow.parquet")
      )

      if ("silico.dart" %in% data.source) {
        silicodata <-
          dplyr::select(silicodata,-c(STRATA)) %>%
          # dplyr::rename(silicodata, VISUAL_STRATA = STRATA) %>%
          dplyr::left_join(y.data, by = "INDIVIDUALS")
        # dplyr::rename(TARGET_ID = TARGET_ID.x) %>%
        # dplyr::select(-c(TARGET_ID.y))
        # dplyr::mutate(GENETIC_STRATA = STRATA)
        genometranslator::write_genome(
          data = silicodata,
          filename = file.path(wd, "sexy_markers_silicodata_genetic_SNP_sex_ID.arrow.parquet")
          )
      } else{
        silicodata <- NULL
      }

    } else if (sex.id.input == 3) {
      readr::write_tsv(
        x = y.silico.data,
        file = file.path(
          path.folder,
          "sexy_markers_strata_new_with_genetic_sex_according_to_SILICOdata.tsv"
        )
      )
      message("New strata file with genetic sex is written.")
      SexID <- "genetically (SILICO)"
      ### recalculate data based on new sexID
      data <-
        dplyr::select(data,-c(STRATA)) %>%
        # dplyr::rename(data, VISUAL_STRATA = STRATA) %>%
        dplyr::left_join(y.silico.data, by = "INDIVIDUALS")
      # dplyr::rename(TARGET_ID = TARGET_ID.x) %>%
      # dplyr::select(-c(TARGET_ID.y))
      # dplyr::mutate(GENETIC_STRATA = STRATA)
      genometranslator::write_genome(
        data = data,
        filename = file.path(wd, "sexy_markers_data_genetic_SILCIO_sex_ID.arrow.parquet")
        )

      if ("silico.dart" %in% data.source) {
        silicodata <-
          dplyr::select(silicodata,-c(STRATA)) %>%
          # dplyr::rename(silicodata, VISUAL_STRATA = STRATA) %>%
          dplyr::left_join(y.silico.data, by = "INDIVIDUALS")
        # dplyr::rename(TARGET_ID = TARGET_ID.x) %>%
        # dplyr::select(-c(TARGET_ID.y))
        # dplyr::mutate(GENETIC_STRATA = STRATA)
        genometranslator::write_genome(
          data = silicodata,
          filename = file.path(wd, "sexy_markers_silicodata_genetic_SILICO_sex_ID.arrow.parquet")
          )
      } else{
        silicodata <- NULL
      }
    }


    if (sex.id.input != 1) {
      message("Sex and summary statistics will be calculated according to: ",
              sex.id.input)
      res$sum <- summarize_sex(
        data = data,
        silicodata = silicodata,
        coverage.thresholds = coverage.thresholds,
        data.source = data.source,
        tau = tau
      )
      data.sum <- res$sum$data.sum
      silico.sum <- res$sum$silico.sum

      # Add new summary for number F and M and ratio, because Unknown are now included
      sex.ratio <-
        dplyr::distinct(data, INDIVIDUALS, .keep_all = TRUE) %>%
        dplyr::count(GENETIC_SEX, name = "N")
      message(
        "\n\nIndividuals with unknown sex ID are now included in the analysis.
        The new sex-ratio (F/M) is: ",
        round(sex.ratio$N[sex.ratio$GENETIC_SEX == "F"] / sex.ratio$N[sex.ratio$GENETIC_SEX == "M"], 2)
      )
      print(sex.ratio)
    } else {message("Sex ID has not changed. Sex and summary statistics are not recalculated.")}
  }


  # Remove markers that are already extracted and have high missingness
  print(ggplot2::qplot(data.sum$MISSINGNESS, xlab = "Missingness per SNP marker"))
  if (interactive.filter) {
    mis.threshold.data <- tgbase::question(
      x = "Have a look at the plot: Choose the upper threshold for missingness per SNP marker (e.g. 0.2):",
      minmax = c(0, 1)
      )

    message(
      "For detecting heterogametic markers, the SILICO data is filtered on a missingness >: ",
      mis.threshold.data
    )
    data.sum <-
      dplyr::filter(data.sum,
                    !(MARKERS %in% y.markers | MISSINGNESS > mis.threshold.data))
  } else {
    message(
      "For detecting heterogametic markers, the SNP data is filtered on a missingness >: ",
      mis.threshold.data)
    data.sum <-
      dplyr::filter(data.sum, !(MARKERS %in% y.markers | MISSINGNESS > mis.threshold.data))
  }

  if (!is.null(silicodata)) {
    print(ggplot2::qplot(silico.sum$MISSINGNESS, xlab = "Missingness per SILICO marker"))
    if (interactive.filter) {
      print(ggplot2::qplot(silico.sum$MISSINGNESS, xlab = "Missingness per SILICO marker"))
      mis.threshold.silicodata <-
        tgbase::question(x = "Have a look at the plot: Choose the upper threshold for missingness per SILICO marker (e.g. 0.2):", minmax = c(0, 1))
      message(
        "For detecting homogametic markers, the SILICO data is filtered on a missingness > ",
        mis.threshold.silicodata
      )
      silico.sum <-
        dplyr::filter(silico.sum, !(MARKERS %in% y.silico.markers |
                                      MISSINGNESS > mis.threshold.silicodata
        ))
    } else {
      message(
        "For detecting homogametic markers, the SILICO data is filtered on a missingness > ",
        mis.threshold.silicodata)
      silico.sum <-
        dplyr::filter(silico.sum, !(MARKERS %in% y.silico.markers | MISSINGNESS > mis.threshold.silicodata))
    }
  }



  #### HET ####

  ### SCATTER
  message("Start finding X- or Z-linked markers. CAUTION: this step is largely affected by structure in your data (e.g. families, populations or species).")

  plot.filename <- stringi::stri_join("3A.sexy_markers_HET_scat_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
    stringi::stri_join(".pdf")



  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "MEAN_HET_F",
    y = "MEAN_HET_M",
    scat = TRUE,
    qreg = FALSE,
    x.title = "Proportion of females",
    y.title = "Proportion of males",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Heterozygosity of each SNP marker between females and males",
    plot.filename = plot.filename,
    path = path.folder
  )
  print(scat.fig)

  if (interactive.filter) {
    het.qr.input <-
      tgbase::question(
        x = "\nHeterozygosity method of SNPs:\nDo you want to plot the residuals for \n(1) X-linked markers: i.e. heterozygous for females\n(2) Z-linked markers: i.e. heterozygous for males",
        answer.opt = c("1", "2"))
    het.qr.input <- as.integer(het.qr.input)
  } else{
    het.qr.input <- as.integer(het.qr.input)
  }

  # Quantilte regression
  if (het.qr.input == 1) {
    plot.filename <- stringi::stri_join("3B.sexy_markers_HET_qr_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    qr.fig <- sex_markers_plot(
      data = data.sum,
      x = "QR_RESIDUALS",
      y = "MEAN_HET_M",
      qreg = TRUE,
      tuckey = FALSE,
      x.title = "Quantile residuals for M~F (<- X/Z ->)",
      y.title = "Proportion of heterozygosity in males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned; tau = ", tau)
      } else {
        paste0(population, ": Sex is ", SexID, " assigned; tau = ", tau)
      },
      title = "Quantile residual plot of each SNP marker between females and males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(qr.fig)

  } else if (het.qr.input == 2) {
    plot.filename <- stringi::stri_join("3B.sexy_markers_HET_qr_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    qr.fig <- sex_markers_plot(
      data = data.sum,
      x = "QR_RESIDUALS",
      y = "MEAN_HET_F",
      qreg = TRUE,
      tuckey = FALSE,
      x.title = "Quantile residuals for M~F (<- X/Z ->)",
      y.title = "Proportion of heterozygosity in females",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned; tau = ", tau)
      } else {
        paste0(population, ": Sex is ", SexID, " assigned; tau = ", tau)
      },
      title = "Quantile residual plot of each SNP marker between females and males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(qr.fig)
  } else {
    message("het.qr.input must be '1' or '2'.")
  }


  message(
    "Files written: '3A.sexy_markers_HET_scat_plot.pdf' & '3B.sexy_markers_HET_qr_plot.pdf'"
  )

  # Interacive selection of threshold
  if (interactive.filter) {
    filter.x.markers <- tgbase::question(x = "Look at the figures: Do you want to select X/Z-linked markers (y/n): ", answer.opt = c("y", "n"))

    if (filter.x.markers == "y") {
      threshold.x.markers.qr <-
        tgbase::question(x = "Choose the threshold for X/Z-linked markers (-1 to 1): ",
                                    minmax = c(-1, 1))
    } else {
      threshold.x.markers.qr <- NULL
    }
  }

  if (!is.null(threshold.x.markers.qr)) {
    if (threshold.x.markers.qr < 0) {
      x.markers <- dplyr::filter(data.sum, QR_RESIDUALS < threshold.x.markers.qr)$MARKERS
    } else{
      x.markers <- dplyr::filter(data.sum, QR_RESIDUALS > threshold.x.markers.qr)$MARKERS
    }
    res$homogametic.het.markers <- x.markers
  } else {
    x.markers <- NULL
    res$homogametic.het.markers <- NULL
  }
  if (length(x.markers) == 0) {
    res$homogametic.het.markers <- NULL
  }



  #### READ DEPTH ####
  ####* For data.sum ####
  if (rlang::has_name(data, "READ_DEPTH")) {
    ### SCATTER
    plot.filename <- stringi::stri_join("4A.sexy_markers_RD_scat_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")


    scat.fig <- sex_markers_plot(
      data = data.sum,
      x = "MEAN_READ_DEPTH_F",
      y = "MEAN_READ_DEPTH_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      RD = TRUE,
      x.title = "Mean coverage for females",
      y.title = "Mean coverage for males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Average coverage of each marker between females and males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(scat.fig)


    ##HIST
    plot.filename <- stringi::stri_join("4B.sexy_markers_RD_hist_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    hist1.fig <- sex_markers_plot(
      data = data.sum,
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male coverage ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females over males coverage for each marker",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(hist1.fig)

    ## INTERACTIVE ZOOM
    if (interactive.filter) {
      zoom.data <-
        tgbase::question(x = "Choose the lower OR upper threshold to subset the histogram (e.g. scaling the plot to a ratio < 0.8 or ratio > 1.2)",minmax = c(-2,2))
    } else {
      zoom.data <- zoom.data
    }
    ##HIST2
    plot.filename <- stringi::stri_join("4C.sexy_markers_RD_hist_subsetted_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    hist2.fig <- sex_markers_plot(
      data =
        if (zoom.data > 1) {
          dplyr::filter(data.sum, RATIO > zoom.data)
        } else if (zoom.data < 1) {
          dplyr::filter(data.sum, RATIO < zoom.data)
        },
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male count ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females counts over males counts for each marker",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(hist2.fig)

    message(
      "Files written: '4A.sexy_markers_RD_scat_plot.pdf' & '4B.sexy_markers_RD_hist_plot.pdf' & '4C.sexy_markers_RD_hist_subsetted_plot.pdf'"
    )


    # Interacive selection of threshold
    if (interactive.filter) {
      filter.x.markers <- tgbase::question(x = "Coverage method of SNPs:\nLook at the figures: Do you want to select X/Z-linked markers (y/n): ", answer.opt = c("y", "n"))

      if (filter.x.markers == "y") {
        threshold.x.markers.RD <-
          tgbase::question(x = "Choose the RATIO threshold for X/Z-linked markers: ", minmax = c(-Inf,Inf) )
      } else {
        threshold.x.markers.RD <- NULL
      }
    }

    if (!is.null(threshold.x.markers.RD)) {
      if (threshold.x.markers.RD > 1) {
        x.markers <- dplyr::filter(data.sum, RATIO > threshold.x.markers.RD)$MARKERS
      } else{
        x.markers <- dplyr::filter(data.sum, RATIO < threshold.x.markers.RD)$MARKERS
      }
      res$homogametic.RD.markers <- x.markers
    } else {
      x.markers <- NULL
      res$homogametic.RD.markers <- NULL
    }
    if (length(x.markers) == 0) {
      res$homogametic.RD.markers <- NULL
    }
  }





  ####* For silico.sum ####
  if (all(c("silico.dart", "counts") %in% data.source)) {
    ### SCATTER
    plot.filename <- stringi::stri_join("5A.sexy_markers_SILICO_RD_scat_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "MEAN_READ_DEPTH_F",
      y = "MEAN_READ_DEPTH_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      RD = TRUE,
      x.title = "Mean coverage for females",
      y.title = "Mean coverage for males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Average coverage of each silico marker between females and males",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(scat.fig)


    ##HIST
    plot.filename <- stringi::stri_join("5B.sexy_markers_SILICO_RD_hist_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    hist1.fig <- sex_markers_plot(
      data = silico.sum,
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male coverage ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females coverage over males coverage for each SILICO",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(hist1.fig)


    ## INTERACTIVE ZOOM
    if (interactive.filter) {
      zoom.silicodata <-
        tgbase::question(x = "Choose the lower OR upper threshold to subset the histogram (e.g. scaling the plot to a ratio < 0.8 or ratio > 1.2)",minmax = c(-2,2))
    } else {
      zoom.silicodata <- zoom.silicodata
    }
    ##HIST2
    plot.filename <- stringi::stri_join("5C.sexy_markers_SILICO_RD_hist_subsetted_plot", species, population, file.date, sep = "_", ignore_null = TRUE) %>%
      stringi::stri_join(".pdf")

    hist2.fig <- sex_markers_plot(
      data =
        if (zoom.silicodata > 1) {
          dplyr::filter(silico.sum, RATIO > zoom.silicodata)
        } else if (zoom.silicodata < 1) {
          dplyr::filter(silico.sum, RATIO < zoom.silicodata)
        },
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male coverage ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females coverage over males coverage for each SILICO",
      plot.filename = plot.filename,
      path = path.folder
    )
    print(hist2.fig)

    message(
      "Files written: '5A.sexy_markers_SILICO_RD_scat_plot.pdf' & '5B.sexy_markers_SILICO_RD_hist_plot.pdf' & '5C.sexy_markers_SILICO_RD_hist_subsetted_plot.pdf'"
    )

    # Interacive selection of threshold
    if (interactive.filter) {
      filter.x.markers <-
        tgbase::question(x = "Coverage method of SILICOs:\nLook at the figures: Do you want to select X/Z-linked markers (y/n): ", answer.opt = c("y", "n"))

      if (filter.x.markers == "y") {
        threshold.x.markers.RD.silico <-
          tgbase::question(x = "Choose the RATIO threshold for X/Z-linked SILICO markers: ", minmax = c(-Inf, Inf))
      } else {
        threshold.x.markers.RD.silico <- NULL
      }
    }

    if (!is.null(threshold.x.markers.RD.silico)) {
      if (threshold.x.markers.RD.silico > 1) {
        x.markers <-
          dplyr::filter(silico.sum, RATIO > threshold.x.markers.RD.silico)$MARKERS
      } else{
        x.markers <-
          dplyr::filter(silico.sum, RATIO < threshold.x.markers.RD.silico)$MARKERS
      }
      res$homogametic.RD.silico.markers <- x.markers
    } else {
      x.markers <- NULL
      res$homogametic.RD.silico.markers <- NULL
    }
    if (length(x.markers) == 0) {
      res$homogametic.RD.silico.markers <- NULL
    }
  }


  # Export -------------------------------------------------------------------

  ## SEX MARKER SUMMARY ##

  ### Get the sequence metadata from the gds object
  if (class(gds.bk)[1] == "SeqVarGDSClass") {
    # Set sex-marker to whitelist and allign the sex-marker method with the markers
    meta <- genometranslator::extract_markers_metadata(
      gds = gds.bk,
      markers.meta.select = c("MARKERS","SEQUENCE", "LOCUS"), #LOCUS is CLONE_ID
      metadata.node = TRUE,
      whitelist = FALSE,
      blacklist = FALSE,
      verbose = FALSE
    ) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)
  }
  if (!is.null(silicodata)) {
    silico <- dplyr::select(silicodata, MARKERS, SEQUENCE) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)
  } else{
    silico <- NULL
  }

  if (!is.null(
    c(
      res$heterogametic.markers,
      res$heterogametic.silico.markers,
      res$homogametic.het.markers,
      res$homogametic.RD.markers,
      res$homogametic.RD.silico.markers
    )
  )) {
    res$sexy.summary <-
      tibble::tibble(
        SEX_MARKERS =
          c(
            res$heterogametic.markers,
            res$heterogametic.silico.markers,
            res$homogametic.het.markers,
            res$homogametic.RD.markers,
            res$homogametic.RD.silico.markers
          ),
        CLONE_ID =
          c(
            meta$LOCUS[meta$MARKERS %in% res$heterogametic.markers],
            res$heterogametic.silico.markers,
            meta$LOCUS[meta$MARKERS %in% res$homogametic.het.markers],
            meta$LOCUS[meta$MARKERS %in% res$homogametic.RD.markers],
            res$homogametic.RD.silico.markers
          ),
        METHOD =
          c(
            rep(
              "presence/absence_method-SNP",
              length(res$heterogametic.markers)
            ),
            rep(
              "presence/absence_method-SILICO",
              length(res$heterogametic.silico.markers)
            ),
            rep(
              "heterozygosity_method-SNP",
              length(res$homogametic.het.markers)
            ),
            rep("Coverage_method-SNP", length(res$homogametic.RD.markers)),
            rep(
              "Coverage_method-SILICO",
              length(res$homogametic.RD.silico.markers)
            )
          ),
        MARKER_TYPE =
          c(
            rep("Heterogametic_sex-marker", length(
              c(
                res$heterogametic.markers,
                res$heterogametic.silico.markers
              )
            )),
            rep("Homogametic_sex-marker", length(
              c(
                res$homogametic.het.markers,
                res$homogametic.RD.markers,
                res$homogametic.RD.silico.markers
              )
            ))
          )
      )
    # TODO add check for when SEQUENCE data is not available
    if (rlang::has_name(meta, "SEQUENCE")) {
      res$sexy.summary %<>% dplyr::mutate(
        SEQUENCE =
          c(
            meta$SEQUENCE[meta$MARKERS %in% SEX_MARKERS[METHOD == "presence/absence_method-SNP"]],
            silico$SEQUENCE[silico$MARKERS %in% SEX_MARKERS[METHOD ==
                                                              "presence/absence_method-SILICO"]],
            meta$SEQUENCE[meta$MARKERS %in% SEX_MARKERS[METHOD ==
                                                          "heterozygosity_method-SNP"]],
            meta$SEQUENCE[meta$MARKERS %in% SEX_MARKERS[METHOD ==
                                                          "Coverage_method-SNP"]],
            silico$SEQUENCE[silico$MARKERS %in% SEX_MARKERS[METHOD ==
                                                              "Coverage_method-SILICO"]]
          )
      )
    }
    message("Summary table of sex-linked markers by method of discovery:")
    print(summary(as.factor(res$sexy.summary$METHOD)))
  } else {
    res$sexy.summary <- NULL
  }


  ## Method-overlap plot
  if (sum(
    !is.null(res$heterogametic.markers),
    !is.null(res$heterogametic.silico.markers),
    !is.null(res$homogametic.het.markers),
    !is.null(res$homogametic.RD.markers),
    !is.null(res$homogametic.RD.silico.markers)) > 1) {
    # Common markers plot
    plot.data <- res$sexy.summary %>%
      dplyr::distinct(SEX_MARKERS, METHOD) %>%
      dplyr::mutate(
        n = rep(1, n()),
        METHOD = stringi::stri_join("METHOD_", METHOD)
        ) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = SEX_MARKERS ~ METHOD,
        value.var = "n",
        fill = 0
      ) %>%
      tibble::as_tibble(.)

    method.cols <- grep("^METHOD_", names(plot.data), value = TRUE)
    plot.data <- dplyr::mutate(
      plot.data,
      dplyr::across(tidyselect::all_of(method.cols), as.logical)
    )

    overlap.plot <- ComplexUpset::upset(
      data = plot.data,
      intersect = method.cols,
      name = "Detection method"
    )
    print(overlap.plot)

    plot.filename <- stringi::stri_join(
      "6.sexy_markers_method_overlap_", file.date, ".pdf"
    )
    plot.filename <- file.path(path.folder, plot.filename)

    ggplot2::ggsave(
      filename = plot.filename,
      plot = overlap.plot,
      width = 10,
      height = 7,
      units = "in"
    )
    if (verbose) message("File written: ", basename(plot.filename))
  } else{
    if (verbose) {
      message("At least two sex-marker methods are required for an overlap plot.")
    }
  }

  ## FASTA file with sex markers for all methods ##
  # TODO add check for when SEQUENCE data is not available

  if (!is.null(
    c(res$heterogametic.markers,
      res$heterogametic.silico.markers,
      res$homogametic.het.markers,
      res$homogametic.RD.markers,
      res$homogametic.RD.silico.markers
    )) & (rlang::has_name(meta, "SEQUENCE"))) {
    afile <-
      file(file.path(path.folder, "7.sexy_markers_sequences.fasta"),
           open = 'w')
    for (i in 1:length(res$sexy.summary$SEX_MARKERS)) {
      cat(
        paste(
          '>',
          res$sexy.summary$MARKER_TYPE[i],
          "|",
          res$sexy.summary$METHOD[i],
          "|",
          res$sexy.summary$SEX_MARKERS[i],
          '\n',
          res$sexy.summary$SEQUENCE[i],
          '\n',
          sep = ''
        ),
        sep = '',
        file = afile
      )
    }
    close(afile)
    message("File written:'7.sexy_markers_sequences.fasta'")
  }

  return(res)
}#End sexy_markers


# INTERNAL NESTED FUNCTION------------------------------------------------------
#' @title sex_markers_plot
#' @description Function to generate the different figures required for 'sexy_markers'.
#' @keywords internal
#' @export

sex_markers_plot <- function(
  data, x, y, x.title = "x title", y.title = "y title", subtitle = "subtitle", title = "Big title",
  width = 15, height = 15, scat = FALSE, tuckey = FALSE, qreg = FALSE, hist = FALSE, RD = FALSE, plot.filename = NULL, path = NULL) {

  x <- dplyr::sym(x)

  element.text <- ggplot2::element_text(size = 10, face = "bold")

  if (qreg) data %<>% dplyr::filter(!is.na(QR_RESIDUALS))
  if (hist) {
    sex.plot <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = !!x)
    ) +
      ggplot2::geom_histogram(
        col = "black",
        fill = "grey",
        alpha = .2,
        binwidth = 0.1)
  } else{
    y <- dplyr::sym(y)
    sex.plot <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = !!x, y = !!y)) +
      ggplot2::geom_point(size = 2, alpha = 0.3)

    if (tuckey) {
      sex.plot <- sex.plot + ggplot2::scale_y_reverse()
    } else if (qreg) {
      sex.plot <- sex.plot
    } else if (RD) {
      sex.plot <- sex.plot + ggplot2::geom_abline(slope = c(2,1,1/2))
    } else if (scat) {
      sex.plot <- sex.plot + ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::scale_y_continuous(limits = c(0, 1))
    } else {
      sex.plot <- sex.plot + ggplot2::geom_abline()
    }
  }

  sex.plot <- sex.plot +
    ggplot2::labs(x = x.title, y = y.title, subtitle = subtitle, title = title) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      axis.title.y = element.text,
      axis.title.x = element.text,
      axis.text.x = element.text)
  # print(sex.plot)

  ggplot2::ggsave(
    filename = plot.filename,
    path = path,
    plot = sex.plot,
    width = width,
    height = height,
    dpi = 300,
    units = "cm",
    useDingbats = FALSE
  )
  return(sex.plot)
}

#' @title summarize_sex
#' @description Function to summarize the presence-absence, heterozygosity and read depth per sex.
#' @keywords internal
#' @export

summarize_sex <- function(data, silicodata, data.source, coverage.thresholds = NULL, tau = 0.3) {
  if (rlang::has_name(data, "READ_DEPTH")) {
    # With DArT count data and VCFs

    mis <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, READ_DEPTH, ALT_DOSAGE) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[is.na(ALT_DOSAGE)]) / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    data.sum <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, READ_DEPTH, ALT_DOSAGE) %>%
      dplyr::group_by(STRATA, MARKERS) %>%
      dplyr::summarise(
        PRESENCE_ABSENCE = length(INDIVIDUALS[READ_DEPTH < coverage.thresholds | is.na(READ_DEPTH)]) / length(INDIVIDUALS), #added | is.na(READ_DEPTH)
        MEAN_HET = length(INDIVIDUALS[ALT_DOSAGE == 1]) / length(INDIVIDUALS),
        MEAN_READ_DEPTH = mean(READ_DEPTH, na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ STRATA,
        value.var = c("PRESENCE_ABSENCE", "MEAN_READ_DEPTH", "MEAN_HET")
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
        DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F,
        RATIO = MEAN_READ_DEPTH_F / MEAN_READ_DEPTH_M
      ) %>%
      dplyr::filter(RATIO != 0) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(mis, by = "MARKERS")


    data.sum[is.na(data.sum)] <- NA
    mis <- NULL

    if (is.null(tau)) tau <- 0.03
    #Quantile regression plot
    data.sum <- dplyr::bind_rows(
      dplyr::filter(data.sum, is.na(MEAN_HET_M) | is.na(MEAN_HET_F)),
      dplyr::filter(data.sum, !is.na(MEAN_HET_M) &
                      !is.na(MEAN_HET_F)) %>%
        dplyr::mutate(
          QR_RESIDUALS = quantreg::rq(MEAN_HET_M ~ MEAN_HET_F,
                                      tau = tau,
                                      data = .) %$%
            residuals
        )
    )
  } else if (rlang::has_name(data, "ALT_DOSAGE") ) {  #genotype data
    mis <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[is.na(ALT_DOSAGE)]) / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    data.sum <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE) %>%
      dplyr::group_by(STRATA, MARKERS) %>%
      dplyr::summarise(
        PRESENCE_ABSENCE = length(INDIVIDUALS[is.na(ALT_DOSAGE)]) / length(INDIVIDUALS),
        MEAN_HET = length(INDIVIDUALS[ALT_DOSAGE == 1]) / length(INDIVIDUALS)
      ) %>%
      dplyr::ungroup(.) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ STRATA,
        value.var = c("PRESENCE_ABSENCE", "MEAN_HET")
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
        DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(mis, by = "MARKERS")


    data.sum[is.na(data.sum)] <- NA
    mis <- NULL

    #Quantile regression plot
    data.sum <- dplyr::bind_rows(
      dplyr::filter(data.sum, is.na(MEAN_HET_M) | is.na(MEAN_HET_F)),
      dplyr::filter(data.sum, !is.na(MEAN_HET_M) &
                      !is.na(MEAN_HET_F)) %>%
        dplyr::mutate(
          QR_RESIDUALS = quantreg::rq(MEAN_HET_M ~ MEAN_HET_F,
                                      tau = tau,
                                      data = .) %$%
            residuals
        )
    )
  } else {
    data.sum <- NULL
  }


  if ("silico.dart" %in% data.source) {
    mis <-
      dplyr::select(silicodata, MARKERS, INDIVIDUALS, STRATA, VALUE) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[VALUE < coverage.thresholds])
                       / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    ### silico counts
    if ("counts" %in% data.source) {
      silico.sum <- silicodata %>%
        dplyr::group_by(STRATA, MARKERS) %>%
        dplyr::summarise(
          PRESENCE_ABSENCE = length(INDIVIDUALS[VALUE < coverage.thresholds])
          / length(INDIVIDUALS[!is.na(VALUE)]),
          MEAN_READ_DEPTH = mean(VALUE, na.rm = TRUE)
        ) %>%
        dplyr::ungroup(.) %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula = MARKERS ~ STRATA ,
          value.var = c("PRESENCE_ABSENCE", "MEAN_READ_DEPTH"),
          fill = NA
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(
          RATIO = MEAN_READ_DEPTH_F / MEAN_READ_DEPTH_M,
          MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
          DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
        ) %>%
        dplyr::left_join(mis, by = "MARKERS")
    } else {
      silico.sum <- silicodata %>%
        dplyr::group_by(STRATA, MARKERS) %>%
        # dplyr::summarise(PRESENCE_ABSENCE = length(INDIVIDUALS[VALUE < coverage.thresholds])
        #                  / length(INDIVIDUALS[!is.na(VALUE)])) %>%
        dplyr::summarise(PRESENCE_ABSENCE = length(INDIVIDUALS[VALUE[!is.na(VALUE)] < coverage.thresholds])
                         / length(INDIVIDUALS[!is.na(VALUE)])) %>% #na.rm = TRUE
        dplyr::ungroup(.) %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula = MARKERS ~ STRATA,
          value.var = "PRESENCE_ABSENCE",
          fill = NA
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::rename(
          PRESENCE_ABSENCE_M = M,
          PRESENCE_ABSENCE_F = F
        ) %>%
        dplyr::mutate(
          MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
          DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
        ) %>%
        dplyr::left_join(mis, by = "MARKERS")
    }
  } else {
    silico.sum <- NULL
  }
  return(list(data.sum = data.sum, silico.sum = silico.sum))
}


#' @title Export FASTA
#' @description Function to write fasta files for sex markers
#' @keywords internal
#' @export
# FASTA from gds


# FASTA file (different for silico)
write_fasta <- function(sexmarkdf, filename) {
  afile <- file(filename, open = 'w')
  if (!is.null(species) && !is.null(population)) {
    for (i in 1:length(sexmarkdf$MARKERS)) {
      cat(
        paste(
          '>',Species,"|",population,'|',sexmarkdf$method[i],
          sexmarkdf$MARKERS[i],'\n',
          sexmarkdf$SEQUENCE[i],'\n',
          sep = ''
        ),
        sep = '',
        file = afile
      )
    }
  } else {
    for (i in 1:length(sexmarkdf$MARKERS)) {
      cat(
        paste(
          '>', sexmarkdf$method[i],"|",sexmarkdf$MARKERS[i],'\n',
          sexmarkdf$SEQUENCE[i],'\n',
          sep = ''
        ),
        sep = '',
        file = afile
      )
    }
  }
  close(afile)
} #write_fasta
