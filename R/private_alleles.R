# private_alleles

#' @name private_alleles
#' @title Find private alleles
#' @description The function highlight private alleles by strata, using a GDS or tidy file or object.

#' @inheritParams genometranslator::read_strata
#' @inheritParams radr_common_arguments

#' @return A list with an object highlighting private alleles by markers and strata and
#' a second object with private alleles summarized by strata.

#' @examples
#' \dontrun{
#' corals.private.alleles.by.pop <- radr::private_alleles(data = tidy, strata = strata.pop)
#' }



#' @rdname private_alleles
#' @export

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

private_alleles <- function(data, strata = NULL, verbose = TRUE) {

  # Cleanup---------------------------------------------------------------------
  tgbase::function_header(f.name = "private_alleles", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- tgbase::tic()
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(tgbase::toc(timing), add = TRUE)
  on.exit(tgbase::function_header(f.name = "private_alleles", start = FALSE, verbose = verbose), add = TRUE)


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Detect format --------------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)
  if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data ...")
  if (data.type %in% c("tbl_df", "fst.file")) {
    if (is.vector(data)) data <- genometranslator::read_genome(data, import.metadata = TRUE)
    data.type <- "tbl_df"

    if (is.null(strata)) {
      if (rlang::has_name(data, "POP_ID") || rlang::has_name(data, "STRATA")) {
        if (rlang::has_name(data, "POP_ID")) data %<>% dplyr::rename(STRATA = POP_ID)
        strata <- genometranslator::generate_strata(data = data)
      } else {
        rlang::abort("strata required and if not provided a STRATA or POP_ID column in the data...")
      }
    } else {
      strata <- genometranslator::read_strata(strata = strata) %$% strata
      data %<>% genometranslator::join_strata(data = ., strata = strata)
    }

    if (rlang::has_name(data, "GT_VCF_NUC")) {
      private <- dplyr::filter(data, GT_VCF_NUC != "./.") %>%
        dplyr::distinct(STRATA, MARKERS, GT_VCF_NUC) %>%
        genometranslator::separate_gt(x = ., gt = "GT_VCF_NUC", exclude = c("MARKERS", "STRATA"), split.chunks = 1L) %>%
        dplyr::select(-ALLELES_GROUP) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELES)
    } else if (rlang::has_name(data, "GT")) {
      private <- dplyr::filter(data, GT != "000000") %>%
        dplyr::distinct(STRATA, MARKERS, GT) %>%
        genometranslator::separate_gt(x = ., gt = "GT",exclude = c("MARKERS", "STRATA"), split.chunks = 1L) %>%
        dplyr::select(-ALLELES_GROUP) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELES)
    } else {
      # work with ALT_DOSAGE
      private <- dplyr::filter(data, !is.na(ALT_DOSAGE)) %>%
        dplyr::distinct(STRATA, MARKERS, ALT_DOSAGE) %>%
        genometranslator::separate_gt(x = ., gt = "ALT_DOSAGE",exclude = c("MARKERS", "STRATA"), split.chunks = 1L) %>%
        dplyr::select(-ALLELES_GROUP) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELES)

        #
        # dplyr::mutate(
        #   ALLELE = dplyr::case_when(
        #     ALT_DOSAGE == 0 ~ "REFREF",
        #     ALT_DOSAGE == 1 ~ "REFALT",
        #     ALT_DOSAGE == 2 ~ "ALTALT"
        #   ),
        #   ALT_DOSAGE = NULL
        # ) %>%
        # genometranslator::separate_gt(x = ., gt = "ALLELE", exclude = c("MARKERS", "STRATA")) %>%
        # dplyr::select(-ALLELE_GROUP, ALLELE = HAPLOTYPES) %>%
        # dplyr::distinct(MARKERS, STRATA, ALLELE)
    }
  }

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- genometranslator::read_genome(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
    }
    strata <- genometranslator::extract_individuals_metadata(
      gds = data,
      ind.field.select = c("INDIVIDUALS", "STRATA"),
      whitelist = TRUE
    )
    private <- genometranslator::generate_gt_vcf_nuc(data) %>%
      magrittr::set_rownames(x = ., value = SeqArray::seqGetData(gdsobj = data, "sample.id")) %>%
      magrittr::set_colnames(
        x = .,
        value = genometranslator::extract_markers_metadata(
          gds = data,
          markers.meta.select = "MARKERS",
          whitelist = TRUE
          ) %$%
          MARKERS
      ) %>%
      tibble::as_tibble(x = ., rownames = "INDIVIDUALS") %>%
      tgbase::trans_long(
        x = .,
        cols = "INDIVIDUALS",
        names_to = "MARKERS",
        values_to = "GT_VCF_NUC",
        variable_factor = FALSE
      ) %>%
      genometranslator::separate_gt(x = ., gt = "GT_VCF_NUC", exclude = c("MARKERS", "STRATA"), split.chunks = 1L) %>%
      dplyr::select(-ALLELES_GROUP)

    private %<>%
      genometranslator::join_strata(data = ., strata = strata, verbose = FALSE) %>%
      dplyr::distinct(MARKERS, STRATA, ALLELES)
  }

  data <- NULL
  private.search <- private %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::tally(.) %>%
    dplyr::filter(n == 1) %>%
    dplyr::distinct(MARKERS, ALLELES) %>%
    dplyr::left_join(private, by = c("MARKERS", "ALLELES")) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(STRATA, MARKERS, ALLELES) %>%
    dplyr::arrange(STRATA, MARKERS, ALLELES) %>%
    readr::write_tsv(x = ., file = "private.alleles.tsv")

  private.summary <- private.search %>%
    dplyr::group_by(STRATA) %>%
    dplyr::tally(.) %>%
    dplyr::ungroup(.) %>%
    tibble::add_row(.data = ., STRATA = "OVERALL", n = sum(.$n)) %>%
    dplyr::rename(PRIVATE_ALLELES = n) %>%
    readr::write_tsv(x = ., file = "private.alleles.summary.tsv")

  if (nrow(private.summary) > 0) {
    message("Number of private alleles per strata:")
    priv.message <- dplyr::mutate(private.summary, PRIVATE = stringi::stri_join(STRATA, PRIVATE_ALLELES, sep = " = "))
    message(stringi::stri_join(priv.message$PRIVATE, collapse = "\n"))
    res <- list(private.alleles = private.search, private.alleles.summary = private.summary)
  } else {
    message("Private allele(s) found: 0")
    res <- "0 private allele"
  }

  return(res)
}#End private_alleles
