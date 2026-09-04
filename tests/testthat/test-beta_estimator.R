beta_gds_fixture <- function(genome) {
  skip_if_not_installed("SeqArray")
  prefix <- tempfile()
  vcf <- paste0(prefix, ".vcf")
  gds <- paste0(prefix, ".gds")
  samples <- unique(genome$INDIVIDUALS)
  markers <- unique(genome$MARKERS)
  lines <- c(
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=10000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", samples), collapse = "\t")
  )
  for (i in seq_along(markers)) {
    chunk <- genome[genome$MARKERS == markers[i], ]
    dose <- chunk$ALT_DOSAGE[match(samples, chunk$INDIVIDUALS)]
    gt <- c("0/0", "0/1", "1/1")[dose + 1L]
    gt[is.na(gt)] <- "./."
    lines <- c(lines, paste(c("1", i, markers[i], "A", "C", ".",
                            "PASS", ".", "GT", gt), collapse = "\t"))
  }
  writeLines(lines, vcf)
  SeqArray::seqVCF2GDS(vcf, gds, verbose = FALSE)
  withr::defer(unlink(c(vcf, gds)), envir = parent.frame())
  gds
}

beta_gds_estimate <- function(genome, ...) {
  beta_estimator(beta_gds_fixture(genome),
                 strata = dplyr::distinct(genome, INDIVIDUALS, STRATA), ...)
}

test_that("GDS filters and caller-owned handles survive success and failure", {
  genome <- tidyr::crossing(
    MARKERS = c("m1", "m2"),
    INDIVIDUALS = c("a1", "a2", "b1", "b2", "c1", "c2")
  ) |>
    dplyr::mutate(STRATA = substr(INDIVIDUALS, 1, 1),
                  ALT_DOSAGE = rep(c(0, 1, 1, 2, 0, 2), 2))
  path <- beta_gds_fixture(genome)
  strata <- dplyr::distinct(genome, INDIVIDUALS, STRATA)
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, sample.id = c("a1", "a2", "b1", "b2"),
                        variant.id = 1L, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  result <- beta_estimator(gds, strata = strata, verbose = FALSE)
  expect_equal(result$beta$N_MARKERS, c(1L, 1L))
  expect_equal(sort(result$beta$STRATA), c("a", "b"))
  expect_equal(SeqArray::seqGetFilter(gds), before)
  expect_error(beta_estimator(gds, strata = strata, populations = "a",
                              verbose = FALSE), "two populations")
  expect_equal(SeqArray::seqGetFilter(gds), before)
  expect_equal(SeqArray::seqGetData(gds, "variant.id"), 1L)
})

test_that("path input closes the handle it opens", {
  genome <- tibble::tibble(
    MARKERS = "m1", INDIVIDUALS = c("a1", "a2", "b1", "b2"),
    STRATA = c("A", "A", "B", "B"), ALT_DOSAGE = c(0, 1, 1, 2)
  )
  path <- beta_gds_fixture(genome)
  strata <- dplyr::distinct(genome, INDIVIDUALS, STRATA)
  beta_estimator(path, strata = strata, verbose = FALSE)
  gds <- SeqArray::seqOpen(path, readonly = FALSE)
  expect_s4_class(gds, "SeqVarGDSClass")
  SeqArray::seqClose(gds)
  expect_error(beta_estimator(path, strata = strata, populations = "A",
                              verbose = FALSE), "two populations")
  gds <- SeqArray::seqOpen(path, readonly = FALSE)
  expect_s4_class(gds, "SeqVarGDSClass")
  SeqArray::seqClose(gds)
})

test_that("beta_estimator uses a common marker set", {
  genome <- tibble::tribble(
    ~MARKERS, ~INDIVIDUALS, ~STRATA, ~ALT_DOSAGE,
    "m1", "a1", "A", 0,
    "m1", "a2", "A", 1,
    "m1", "b1", "B", 2,
    "m1", "b2", "B", 1,
    "m2", "a1", "A", 0,
    "m2", "a2", "A", 0,
    "m2", "b1", "B", NA,
    "m2", "b2", "B", NA
  )

  result <- beta_gds_estimate(genome, verbose = FALSE, chunk.size = 1L)

  expect_equal(unique(result$beta$N_MARKERS), 1L)
  expect_equal(nrow(result$between_populations), 1L)
  expect_equal(sort(unique(result$within_population$STRATA)), c("A", "B"))
})

test_that("population selection creates an explicit two-population analysis", {
  genome <- tidyr::crossing(
    MARKERS = c("m1", "m2"),
    INDIVIDUALS = c("a1", "a2", "b1", "b2", "c1", "c2")
  ) |>
    dplyr::mutate(
      STRATA = rep(c("A", "B", "C"), each = 2L)[
        match(INDIVIDUALS, c("a1", "a2", "b1", "b2", "c1", "c2"))
      ],
      ALT_DOSAGE = rep(c(0, 1, 1, 2, 0, 2), 2L)
    )

  result <- beta_estimator(
    beta_gds_fixture(genome),
    strata = dplyr::distinct(genome, INDIVIDUALS, STRATA),
    populations = c("A", "C"),
    verbose = FALSE
  )

  expect_equal(sort(result$beta$STRATA), c("A", "C"))
  expect_equal(nrow(result$between_populations), 1L)
  expect_equal(result$between_populations$STRATA_1, "A")
  expect_equal(result$between_populations$STRATA_2, "C")
})

test_that("beta_estimator rejects non-GDS input", {
  expect_error(beta_estimator(data.frame(), verbose = FALSE), "GDS")
  expect_error(beta_estimator("input.vcf", verbose = FALSE), "GDS")
  expect_error(beta_estimator("missing.gds", verbose = FALSE), "does not exist")
})

test_that("beta_estimator reproduces the hierfstat fs.dosage estimator", {
  # Fixed benchmark independently cross-checked against hierfstat::fs.dosage().
  genome <- tibble::tribble(
    ~MARKERS, ~INDIVIDUALS, ~STRATA, ~ALT_DOSAGE,
    "m1", "a1", "A", 0,  "m1", "a2", "A", 1,
    "m1", "b1", "B", 1,  "m1", "b2", "B", 2,
    "m1", "c1", "C", 0,  "m1", "c2", "C", 2,
    "m2", "a1", "A", 0,  "m2", "a2", "A", 0,
    "m2", "b1", "B", 1,  "m2", "b2", "B", 1,
    "m2", "c1", "C", 2,  "m2", "c2", "C", 2,
    "m3", "a1", "A", 1,  "m3", "a2", "A", 2,
    "m3", "b1", "B", 0,  "m3", "b2", "B", 1,
    "m3", "c1", "C", NA, "m3", "c2", "C", 1
  )

  result <- beta_gds_estimate(genome, verbose = FALSE)
  observed <- stats::setNames(result$beta$BETA, result$beta$STRATA)

  expect_equal(
    unname(observed[c("A", "B", "C")]),
    c(0.4418604651162791, 0.1627906976744186, 0.1627906976744186)
  )
})

test_that("random.mating reproduces the hierfstat betas estimator", {
  genome <- tibble::tribble(
    ~MARKERS, ~INDIVIDUALS, ~STRATA, ~ALT_DOSAGE,
    "m1", "a1", "A", 0,  "m1", "a2", "A", 1,
    "m1", "b1", "B", 1,  "m1", "b2", "B", 2,
    "m1", "c1", "C", 0,  "m1", "c2", "C", 2,
    "m2", "a1", "A", 0,  "m2", "a2", "A", 0,
    "m2", "b1", "B", 1,  "m2", "b2", "B", 1,
    "m2", "c1", "C", 2,  "m2", "c2", "C", 2,
    "m3", "a1", "A", 1,  "m3", "a2", "A", 2,
    "m3", "b1", "B", 0,  "m3", "b2", "B", 1,
    "m3", "c1", "C", NA, "m3", "c2", "C", 1
  )

  result <- beta_estimator(
    beta_gds_fixture(genome),
    strata = dplyr::distinct(genome, INDIVIDUALS, STRATA),
    random.mating = TRUE,
    chunk.size = 1L,
    verbose = FALSE
  )
  observed <- stats::setNames(result$beta$BETA, result$beta$STRATA)

  expect_equal(
    unname(observed[c("A", "B", "C")]),
    c(3 / 7, 1 / 21, 1 / 21)
  )
  expect_true(all(result$beta$ESTIMATOR == "ALLELE_FREQUENCY_RANDOM_MATING"))
})

test_that("beta_estimator requires STRATA", {
  genome <- tibble::tribble(
    ~MARKERS, ~INDIVIDUALS, ~GROUP, ~ALT_DOSAGE,
    "m1", "a1", "A", 0,
    "m1", "a2", "A", 1,
    "m1", "b1", "B", 1,
    "m1", "b2", "B", 2
  )

  expect_error(beta_estimator(beta_gds_fixture(genome), verbose = FALSE),
               "strata|STRATA|metadata")
})
