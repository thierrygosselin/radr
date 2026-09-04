# strip_rad --------------------------------------------------------------------
#' @rdname strip_rad
#' @title strip_rad
#' @description Strip a tidy data set of it's strata and markers meta. Used internally.
#' @param x The data
#' @param m (character, string) The variables part of the markers metadata.
#' @param env.arg You want to redirect \code{rlang::current_env()} to this argument.
#' @param keep.strata (logical) Keep the strata in the dataset or remove the info
#' and keep only the sample ids.
#' @param verbose (logical) The function will chat more when allowed.
# @noRd
#' @keywords internal
#' @export
strip_rad <- function(
    x,
    m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
    env.arg = NULL,
    keep.strata = TRUE,
    verbose = TRUE
) {
  objs <- utils::object.size(x)

  # If ID_SEQ, STRATA_SEQ and M_SEQ already present... remove...

  if (rlang::has_name(x, "ID_SEQ") && rlang::has_name(x, "INDIVIDUALS")) {
    x %<>% dplyr::select(-ID_SEQ)
  }

  if (rlang::has_name(x, "STRATA_SEQ")) {
    strata.n <- intersect(colnames(x), c("STRATA", "POP_ID"))
    if (length(strata.n) > 0 && strata.n %in% c("STRATA", "POP_ID")) {
      x %<>% dplyr::select(-STRATA_SEQ)
    }
  }

  if (rlang::has_name(x, "M_SEQ")) {
    markers.n <- intersect(colnames(x), c("VARIANT_ID", "CHROM", "LOCUS", "POS", "MARKERS"))
    if (length(markers.n) > 0) {
      x %<>% dplyr::select(-M_SEQ)
    }
  }
  strata.n <- markers.n <- NULL


  # STRATA ----------
  if (rlang::has_name(x, "POP_ID")) {
    strata <- dplyr::ungroup(x) |>
      dplyr::transmute(POP_ID = .data[[if ("POP_ID" %in% names(x))
        "POP_ID" else "STRATA"]], INDIVIDUALS = .data$INDIVIDUALS) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$POP_ID, .data$INDIVIDUALS) %>%
      dplyr::mutate(
        ID_SEQ = seq_len(length.out = dplyr::n()),
        STRATA_SEQ = as.integer(factor(x = POP_ID, levels = unique(POP_ID)))
      )
    pop.levels <- as.character(unique(strata$POP_ID))
  } else {#STRATA
    strata <- dplyr::ungroup(x) |>
      dplyr::transmute(STRATA = .data[[if ("POP_ID" %in% names(x))
        "POP_ID" else "STRATA"]], INDIVIDUALS = .data$INDIVIDUALS) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$STRATA, .data$INDIVIDUALS) %>%
      dplyr::mutate(
        ID_SEQ = seq_len(length.out = dplyr::n()),
        STRATA_SEQ = as.integer(factor(x = STRATA, levels = unique(STRATA)))
      )
    pop.levels <- as.character(unique(strata$STRATA))
  }


  cm <- intersect(colnames(strata), colnames(x))
  x %<>%
    dplyr::left_join(strata, by = cm) %>%
    dplyr::select(-tidyselect::any_of(cm))

  if (!keep.strata) x %<>% dplyr::select(-STRATA_SEQ)

  # Note to myself: maybe write using vroom and read back when necessary ?
  assign(
    x = "strata.bk",
    value = strata,
    pos = env.arg,
    envir = env.arg
  )
  cm <- keep.strata <- pop.id <- strata <- NULL


  assign(
    x = "pop.levels.bk",
    value = pop.levels,
    pos = env.arg,
    envir = env.arg
  )

  # MARKERS ---------
  x %<>%
    dplyr::mutate(
      M_SEQ = as.integer(factor(x = MARKERS, levels = unique(MARKERS)))
    )
  m <- c("M_SEQ", m)
  markers.meta <- x %>%
    dplyr::select(tidyselect::any_of(m)) %>%
    dplyr::distinct(M_SEQ, .keep_all = TRUE)

  cm <- intersect(colnames(markers.meta), colnames(x)) %>%
    purrr::discard(.x = ., .p = . %in% "M_SEQ")

  # remove markers meta
  x %<>% dplyr::select(-tidyselect::any_of(cm))

  # Note to myself: maybe write using vroom and read back when necessary ?
  assign(
    x = "markers.meta.bk",
    value = markers.meta,
    pos = env.arg,
    envir = env.arg
  )
  markers.meta <- cm <- m <- NULL

  # GENOTYPE META ---------

   # g <- c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH",
  # "GL", "CATG", "PL", "HQ", "GQ", "GOF","NR", "NV")
  want <- c("GT", "GT_VCF_NUC", "GT_VCF", "ALT_DOSAGE", "REF", "ALT", "MARKERS",
            "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID"
  )
  g <- purrr::keep(.x = colnames(x), .p = !colnames(x) %in% want) %>%
    purrr::discard(.x = ., .p = . %in% c("STRATA_SEQ"))
  # purrr::discard(.x = ., .p = . %in% c("ID_SEQ", "STRATA_SEQ", "M_SEQ"))
  # g
  genotypes.meta <- NULL # default

  if (length(g) != 0L) {
    genotypes.meta <- x %>% dplyr::select(tidyselect::any_of(g))
    want <- c("ID_SEQ", "STRATA_SEQ", "M_SEQ", "GT", "GT_VCF_NUC", "GT_VCF", "ALT_DOSAGE", "REF", "ALT")
    x %<>% dplyr::select(tidyselect::any_of(want))
  }

  assign(
    x = "genotypes.meta.bk",
    value = genotypes.meta,
    pos = env.arg,
    envir = env.arg
  )
  g <- want <- genotypes.meta <- NULL

  if (verbose) message("Proportion of size reduction: ", round(1 - (utils::object.size(x) / objs), 2))
  return(x)
  # return(list(data = x, markers.meta = markers.meta, genotypes.meta = genotypes.meta))
}# End strip_rad


# join_rad ---------------------------------------------------------------------
#' @rdname join_rad
#' @title join_rad
#' @description Join back the parts stripped in strip_rad. Used internally.
#' @param x The data
#' @param s The strata metadata.
#' @param m The markers metadata.
#' @param g The genotypes metadata.
#' @param env.arg You want to redirect \code{rlang::current_env()} to this argument.
# @noRd
#' @keywords internal
#' @export
join_rad <- function(x, s, m, g, env.arg = NULL) {
  if (!is.null(g)) {
    g.by <- intersect(colnames(x), colnames(g))
    x %<>% dplyr::left_join(g, by = g.by)
    env.arg$genotypes.meta.bk <- NULL
  }

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
            "COL", "REF", "ALT", "INDIVIDUALS", "POP_ID", "STRATA",
            "GT_VCF", "GT_VCF_NUC", "GT", "ALT_DOSAGE")
  s.by <- intersect(colnames(x), colnames(s))
  m.by <- intersect(colnames(x), colnames(m))
  x %<>%
    dplyr::left_join(s, by = s.by) %>%
    dplyr::select(-tidyselect::any_of(c("ID_SEQ", "STRATA_SEQ"))) %>%
    dplyr::left_join(m, by = m.by) %>%
    dplyr::select(-M_SEQ) %>%
    dplyr::select(tidyselect::any_of(want), tidyselect::everything())

  env.arg$strata.bk <- env.arg$markers.meta.bk <- NULL
  return(x)
}#End join_rad

