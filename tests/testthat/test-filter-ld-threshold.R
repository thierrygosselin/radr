test_that("squared-correlation thresholds are converted only for SNPRelate", {
  threshold <- 0.5
  expect_equal(.ld_pruning_parameters("r2", threshold),
               list(method = "r", threshold = sqrt(0.5)))
  expect_identical(threshold, 0.5)
  for (method in c("r", "corr", "composite", "dprime")) {
    expect_identical(.ld_pruning_parameters(method, threshold),
                     list(method = method, threshold = threshold))
  }
  expect_equal(.ld_pruning_parameters("r2", 0)$threshold, 0)
  expect_equal(.ld_pruning_parameters("r2", 1)$threshold, 1)
})

test_that("SNPRelate retains a pair below the requested squared-correlation cutoff", {
  skip_if_not_installed("SNPRelate")
  # Two linked biallelic markers with r near 0.6 and r2 near 0.36.
  withr::local_seed(20260921)
  n <- 2000L
  a <- matrix(rbinom(n * 2L, 1L, 0.5), n, 2L)
  b <- a
  flip <- matrix(runif(n * 2L) < 0.2, n, 2L)
  b[flip] <- 1L - b[flip]
  genotype <- rbind(rowSums(a), rowSums(b))
  path <- tempfile(fileext = ".gds")
  on.exit(unlink(path), add = TRUE)
  SNPRelate::snpgdsCreateGeno(path, genotype,
    sample.id = as.character(seq_len(n)), snp.id = 1:2,
    snp.chromosome = c(1L, 1L), snp.position = c(1L, 2L), snpfirstdim = TRUE)
  gds <- SNPRelate::snpgdsOpen(path)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE, after = FALSE)
  r <- SNPRelate::snpgdsLDMat(gds, slide = 1, method = "r", verbose = FALSE)$LD[1, 1]
  expect_gt(abs(r), 0.5)
  expect_lt(r^2, 0.5)
  prune <- function(method, threshold) {
    parameters <- .ld_pruning_parameters(method, threshold)
    unlist(SNPRelate::snpgdsLDpruning(gds, method = parameters$method,
      ld.threshold = parameters$threshold, start.pos = "first", verbose = FALSE),
      use.names = FALSE)
  }
  expect_length(prune("r2", 0.5), 2L)
  expect_length(prune("r", 0.5), 1L)
  expect_length(prune("r2", 0.3), 1L)
})
