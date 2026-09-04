# Internal helpers retained for the miscall model and legacy HWE workflow.
summarise_genotypes <- function(data, path.folder = NULL) {
  if (is.null(path.folder)) path.folder <- getwd()

  # data.bk <- data
  # data <- data.bk
  data %<>% dplyr::rename(POP_ID = tidyselect::any_of("STRATA"))
  want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "ALT_DOSAGE", "READ_DEPTH")
  data %<>% dplyr::select(tidyselect::any_of(want))

  n.pop <- length(unique(data$POP_ID))
  if (is.factor(data$POP_ID)) {
    pop.levels <- c(levels(data$POP_ID), "OVERALL")
  } else {
    pop.levels <- c(sort(unique(data$POP_ID)), "OVERALL")
  }
  data$POP_ID <- as.character(data$POP_ID)
  replace_zero <- function(x) replace(x = x, list = which(is.na(x)), 0)

  # sum of read depth per pop and overall
  # scaled separately for pop and overall before merging
  if (rlang::has_name(data, "READ_DEPTH")) {
    rd.pop <- dplyr::select(data, MARKERS, POP_ID, READ_DEPTH) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(READ_DEPTH_SCALED = READ_DEPTH / max(READ_DEPTH, na.rm = TRUE))

    rd <- dplyr::bind_rows(
      dplyr::select(rd.pop, MARKERS, POP_ID, READ_DEPTH = READ_DEPTH_SCALED),
      dplyr::mutate(rd.pop, POP_ID = "OVERALL") %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(READ_DEPTH = READ_DEPTH / max(READ_DEPTH, na.rm = TRUE))
    ) %>%
      dplyr::arrange(MARKERS, POP_ID) %>%
      dplyr::select(-MARKERS, -POP_ID)
    rd.pop <- NULL
  } else {
    rd <- NULL
  }

  pop <- data %>%
    dplyr::mutate(
      ALT_DOSAGE = dplyr::case_when(
        ALT_DOSAGE == 0 ~ "HOM_REF",
        ALT_DOSAGE == 1 ~ "HET",
        ALT_DOSAGE == 2 ~ "HOM_ALT",
        is.na(ALT_DOSAGE) ~ "MISSING")
    ) %>%
    dplyr::group_by(MARKERS, POP_ID, ALT_DOSAGE) %>%
    dplyr::tally(.) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = MARKERS + POP_ID ~ ALT_DOSAGE,
      value.var = "n"
    ) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(dplyr::across(where(is.integer), .fns = replace_zero))

  if (!rlang::has_name(pop, "HET")) {
    pop %<>% dplyr::mutate(HET = 0)
  }
  if (!rlang::has_name(pop, "HOM_REF")) {
    pop %<>% dplyr::mutate(HOM_REF = 0)
  }
  if (!rlang::has_name(pop, "HOM_ALT")) {
    pop %<>% dplyr::mutate(HOM_ALT = 0)
  }
  pop %<>% dplyr::mutate(N = HOM_REF + HET + HOM_ALT)

  if (!rlang::has_name(pop, "MISSING")) {
    pop <- dplyr::mutate(pop, MISSING = as.integer("0"))
  }

  data <- dplyr::bind_rows(
    pop,
    dplyr::mutate(pop, POP_ID = "OVERALL") %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(dplyr::across(.cols = tidyselect::everything(), .fns = sum))
    ) %>%
    dplyr::arrange(MARKERS, POP_ID)
  pop <- NULL

  data  %<>%
    dplyr::mutate(
      FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
      FREQ_REF = 1 - FREQ_ALT,
      FREQ_HET = HET / (2 * N),
      FREQ_HOM_REF_O = HOM_REF / N,
      FREQ_HET_O = HET / N,
      FREQ_HOM_ALT_O = HOM_ALT / N,
      FREQ_HOM_REF_E = FREQ_REF^2,
      FREQ_HET_E = 2 * FREQ_REF * FREQ_ALT,
      FREQ_HOM_ALT_E = FREQ_ALT^2,
      N_HOM_REF_EXP = N * FREQ_HOM_REF_E,
      N_HET_EXP = N * FREQ_HET_E,
      N_HOM_ALT_EXP = N * FREQ_HOM_ALT_E,
      HOM_REF_Z_SCORE = (HOM_REF - N_HOM_REF_EXP) / sqrt(N * FREQ_HOM_REF_E * (1 - FREQ_HOM_REF_E)),
      HOM_HET_Z_SCORE = (HET - N_HET_EXP) / sqrt(N * FREQ_HET_E * (1 - FREQ_HET_E)),
      HOM_ALT_Z_SCORE = (HOM_ALT - N_HOM_ALT_EXP) / sqrt(N * FREQ_HOM_ALT_E * (1 - FREQ_HOM_ALT_E))
    ) %>%
    dplyr::bind_cols(rd) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE))

  rd <- NULL

  readr::write_tsv(x = data, file = file.path(path.folder, "genotypes.summary.tsv"))
  return(data)
}#End summarise_genotypes

# Generate genotypes frequencies boundaries


plot_het_outliers <- function(data, path.folder = NULL) {
  res <- list() # to store results
  n.pop <- dplyr::n_distinct(data$POP_ID)
  res$het.summary <- summarise_genotypes(data, path.folder = path.folder)

  # prepare data for figure
  freq.summary <- dplyr::bind_cols(
    res$het.summary %>%
      dplyr::select(MARKERS, POP_ID, HOM_REF = FREQ_HOM_REF_O, HET = FREQ_HET_O, HOM_ALT = FREQ_HOM_ALT_O) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -c("POP_ID", "MARKERS"),
        names_to = "GENOTYPES",
        values_to = "OBSERVED"
      ) %>%
      dplyr::arrange(MARKERS, POP_ID),
    res$het.summary %>%
      dplyr::select(MARKERS, POP_ID, HOM_REF = FREQ_HOM_REF_E, HET = FREQ_HET_E, HOM_ALT = FREQ_HOM_ALT_E) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -c("POP_ID", "MARKERS"),
        names_to = "GENOTYPES",
        values_to = "EXPECTED"
      ) %>%
      dplyr::arrange(MARKERS, POP_ID) %>%
      dplyr::select(EXPECTED)
  ) %>%
    dplyr::mutate(GENOTYPES = factor(
      GENOTYPES,
      levels = c("HOM_REF", "HET", "HOM_ALT"),
      labels = c("Homozygote REF allele", "Heterozygote", "Homozygote ALT allele"),
      ordered = TRUE))

  # generate boundaries
  boundaries <- generate_geno_freq_boundaries() %>%
    dplyr::mutate(
      GENOTYPES = factor(
        GENOTYPES,
        levels = c("Homozygote REF allele", "Heterozygote", "Homozygote ALT allele")))

  res$gt.boundaries.plot <- ggplot2::ggplot(freq.summary , ggplot2::aes(x = EXPECTED, y = OBSERVED, colour = GENOTYPES)) +
    ggplot2::geom_jitter(alpha = 0.1, width = 0.01, height = 0.01) +
    ggplot2::geom_polygon(data = boundaries, fill = NA, linetype = "dashed", colour = "black") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "solid") +
    ggplot2::labs(x = "Genotypes (expected frequency) ", y = "Genotypes (observed frequency)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(POP_ID ~ GENOTYPES)

  if (!is.null(path.folder)) {
    ggplot2::ggsave(
      filename = file.path(path.folder, "markers.genotypes.boundaries.png"),
      plot = res$gt.boundaries.plot,
      width = 20,
      height = (n.pop + 1) * 5,# + 1 for overall always present
      dpi = 200, units = "cm",
      # useDingbats = FALSE,
      limitsize = FALSE)
  }
  return(res)
}#End plot_het_outliers



generate_geno_freq_boundaries <- function() {
  # first do it for the homozygote category
  Poe <- seq(0,1, by = 0.005)
  phat <- 1 - sqrt(Poe)
  minPo <- pmax(0, 1 - 2 * phat)
  maxPo <- 1 - phat

  # these are the values for the two homozygote categories
  homo_tib <- tibble::tibble(
    EXPECTED = rep(Poe, 4),
    OBSERVED = rep(c(minPo, maxPo), 2),
    GENOTYPES = as.character(rep(c("Homozygote REF allele", "Homozygote ALT allele"), each = length(Poe) * 2)))

  # now, it should be easy to get the max/min values for heterozygotes.
  # They will occur where one of the homozygotes is min or max.
  P1e <- 2 * phat * (1 - phat)
  maxP1 <- 2 * (1 - phat - minPo)
  minP1 <- 2 * (1 - phat - maxPo)

  het_tib <- tibble::tibble(
    EXPECTED = rep(P1e, 2),
    OBSERVED = c(minP1, maxP1),
    GENOTYPES = "Heterozygote")

  dplyr::bind_rows(homo_tib, het_tib)
}#End generate_geno_freq_boundaries


estimate_m <- function(
  data,
  nreps = 200,
  m_init = stats::runif(1),
  a0 = 0.5,
  a1 = 0.5,
  sm = 0.005
) {
  D <- dplyr::select(data, INDIVIDUALS, MARKERS, ALT_DOSAGE) %>%
    dplyr::group_by(INDIVIDUALS) %>%
    tidyr::pivot_wider(data = ., names_from = "MARKERS", values_from = "ALT_DOSAGE") %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-INDIVIDUALS) %>%
    as.matrix(.)

  stopifnot(m_init > 0 & m_init < 1)

  D[is.na(D)] <- -1

  # get the N variables
  N0 <- colSums(D == 0)
  N1 <- colSums(D == 1)
  N2 <- colSums(D == 2)

  # initialize the Zs to the Ns
  Z0 <- N0
  Z1 <- N1
  Z2 <- N2

  # make some place to return the m values visited
  m <- rep(NA, nreps)
  m[1] <- m_init

  # then do the sweeps
  for (r in 2:nreps) {

    # new estimate of frequency of the "1" allele from Gibbs sampling
    p <- stats::rbeta(n = length(Z0),
                      shape1 = a1 + 2 * Z2 + Z1,
                      shape2 = a0 + 2 * Z0 + Z1)

    # propose then accept or reject a new value for m
    mprop <- m[r - 1] + stats::rnorm(1, 0, sm)
    reject <- TRUE  # reject it unless we don't
    if (mprop > 0 & mprop < 1) {
      numer <- sum(N0 * log((1 - p)^2 + mprop * p * (1 - p)) +
                     N1 * log((1 - mprop) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + mprop * p * (1 - p)))
      denom <- sum(N0 * log((1 - p)^2 + m[r - 1] * p * (1 - p)) +
                     N1 * log((1 - m[r - 1]) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + m[r - 1] * p * (1 - p)))
      if (log(stats::runif(1)) < numer - denom) {
        reject <- FALSE
      }
    }
    if (reject == FALSE) {
      m[r] <- mprop
    } else {
      m[r] <- m[r - 1]
    }

    # new values for Z from Gibbs sampling
    A0 <- stats::rbinom(n = length(N0), size = N0, prob = (m[r] * p) / (1 - p + m[r] * p))
    A2 <- stats::rbinom(n = length(N2), size = N2, prob = (m[r] * (1 - p)) / (p + m[r] * (1 - p)))

    Z0 <- N0 - A0
    Z1 <- N1 + A0 + A2
    Z2 <- N2 - A2

  }
  # return m, and eventually I need to also return the final Zs and the Ns
  # and I may as well return a new 012 file with "corrected" genotypes, which
  # I can make by broadcasting the Zs around, for example...

  # inferring/realizing/simulating genotypes. I can simulate these from their posterior
  # given the estimated allele freq and the observed genotype.  To do this I will cycle
  # over the columns (the snps) in D, and for each one, I will compute the posterior of the
  # the genotype given the observed genotype (only have to for 0's and 2's) and then I will
  # sample from those posteriors.  We have a separate function that does this
  ret <- list()
  ret$simmed_genos <- simulate_genos_from_posterior(D, p, m[nreps])

  # compute an overall genotyping error rate
  diff <- ret$simmed_genos != D
  diff[D == -1] <- NA
  ret$overall_geno_err_est <- mean(diff, na.rm = TRUE)

  # return the trace of m values
  ret$m <- m

  ret
}#End estimate_m


simulate_genos_from_posterior <- function(D, p, m) {
  stopifnot(length(m) == 1)

  glist <- lapply(1:ncol(D), function(i) {
    obs <- D[, i] # the observed genotypes
    pl <- p[i]  # the alle freq at locus i
    post0 <- c(
      (1 - pl) / (1 - pl + m * pl),  # posterior that observed 0 is truly a 0
      (m * pl) / (1 - pl + m * pl)   # posterior that observed 0 is truly a 1
    )
    post2 <- c(
      (m * (1 - pl)) / (pl + m * (1 - pl)),  # posterior that observed 2 is truly a 1
      pl / (pl + m * (1 - pl))               # posterior that observed 2 is truly a 2
    )
    obs[obs == 0] <- sample(x = c(0, 1), size = sum(obs == 0), replace = TRUE, prob = post0)
    obs[obs == 2] <- sample(x = c(1, 2), size = sum(obs == 2), replace = TRUE, prob = post2)
    obs
  })

  # then turn it into a matrix with the same dimensions and dimnames as D
  ret <- matrix(unlist(glist), nrow = nrow(D))
  dimnames(ret) <- dimnames(D)
  ret
}#End simulate_genos_from_posterior
