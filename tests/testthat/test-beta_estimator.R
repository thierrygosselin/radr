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

  result <- beta_estimator(genome, verbose = FALSE, chunk.size = 1L)

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
    genome,
    populations = c("A", "C"),
    verbose = FALSE
  )

  expect_equal(sort(result$beta$STRATA), c("A", "C"))
  expect_equal(nrow(result$between_populations), 1L)
  expect_equal(result$between_populations$STRATA_1, "A")
  expect_equal(result$between_populations$STRATA_2, "C")
})

test_that("beta_estimator rejects invalid dosage", {
  genome <- tibble::tribble(
    ~MARKERS, ~INDIVIDUALS, ~STRATA, ~ALT_DOSAGE,
    "m1", "a", "A", 0,
    "m1", "b", "B", 3
  )

  expect_error(beta_estimator(genome, verbose = FALSE), "0, 1, 2")
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

  result <- beta_estimator(genome, verbose = FALSE)
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
    genome,
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

  expect_error(beta_estimator(genome, verbose = FALSE), "STRATA")
})
