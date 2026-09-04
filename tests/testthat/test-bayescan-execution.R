testthat::test_that("BayeScan process errors stop before results are read", {
  testthat::skip_on_os("windows")
  executable <- Sys.which("false")
  testthat::skip_if(!nzchar(executable))
  folder <- tempfile(pattern = "bayescan path with spaces ")
  dir.create(folder)
  input <- file.path(folder, "input file.txt")
  writeLines(c("[loci]=1", "[populations]=2"), input)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  testthat::expect_error(bayescan_one(data = input, pr_odds = 10,
    path.folder = folder, file.date = "test", bayescan.path = executable),
    "BayeScan failed")
})

testthat::test_that("one individual per population is a supported subsample", {
  m <- data.frame(STRATA = rep(c("A", "B"), each=2),
    INDIVIDUALS = letters[1:4])
  x <- subsampling_data(1, m, subsample=1, random.seed=42)
  testthat::expect_equal(nrow(x), 2)
  testthat::expect_equal(sort(unique(x$STRATA)), c("A","B"))
})
