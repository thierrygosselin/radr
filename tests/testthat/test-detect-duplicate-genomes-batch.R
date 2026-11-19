test_that("batch duplicate filtering applies a supplied threshold without prompting", {
  skip_if_not_installed("amap")
  output <- tempfile("detect-duplicates-")
  dir.create(output)
  previous.dir <- setwd(output)
  on.exit({
    setwd(previous.dir)
    unlink(output, recursive = TRUE)
  }, add = TRUE)
  if (utils::packageVersion("testthat") >= "3.1.7") {
    testthat::local_mocked_bindings(
      question = function(...) stop("Duplicate filtering prompted in batch mode"),
      .package = "tgbase"
    )
  }

  ids <- sprintf("%03d", seq_len(6L))
  dosages <- cbind(
    c(0L, 1L, 2L, 1L, 0L, 2L),
    c(0L, 1L, NA, 1L, 0L, 2L),
    c(2L, 1L, 0L, 0L, 2L, 0L),
    c(2L, 0L, 1L, 1L, 2L, 1L),
    c(0L, 2L, 1L, 0L, 1L, 2L),
    c(1L, 2L, 0L, 2L, 1L, 0L)
  )
  genome <- tidyr::expand_grid(
    INDIVIDUALS = ids,
    POS = seq_len(nrow(dosages))
  ) |>
    dplyr::mutate(
      MARKERS = paste0("m", POS),
      CHROM = "1",
      LOCUS = MARKERS,
      STRATA = rep(c("A", "B"), each = 3L * nrow(dosages)),
      ALT_DOSAGE = as.vector(dosages)
    )

  diagnostics <- suppressMessages(detect_duplicate_genomes(
    genome, interactive.filter = FALSE, dup.threshold = 0,
    parallel.core = 1L, verbose = FALSE
  ))
  expect_setequal(diagnostics$INDIVIDUALS, ids)
  expect_length(list.files(
    pattern = "^individuals\\.pairwise\\.dist\\.tsv$", recursive = TRUE
  ), 1L)
  expect_length(list.files(
    pattern = "^blacklist\\.id\\.similar\\.tsv$", recursive = TRUE
  ), 0L)

  filtered <- suppressMessages(detect_duplicate_genomes(
    genome, interactive.filter = FALSE, dup.threshold = 0.01,
    parallel.core = 1L, verbose = FALSE
  ))
  expect_setequal(filtered$INDIVIDUALS, ids[-2L])
  blacklist <- list.files(
    pattern = "^blacklist\\.id\\.similar\\.tsv$", recursive = TRUE
  )
  expect_length(blacklist, 1L)
  expect_identical(readr::read_tsv(blacklist, col_types = "c")$INDIVIDUALS, "002")
})
