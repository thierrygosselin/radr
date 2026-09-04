test_that("normalization separates library size from local depth and preserves NA", {
  x <- matrix(20, 12, 20)
  x[1:3, 1] <- 40
  x[4, 2] <- 0
  x[5, 2] <- NA
  x[7:12, ] <- x[7:12, ] * 3
  out <- .cnv_normalize(x, 6, 10)
  expect_equal(out$relative[1:3, 1], rep(2, 3))
  expect_equal(out$relative[4, 2], 0)
  expect_true(is.na(out$relative[5, 2]))
  expect_equal(out$relative[, 3], rep(1, 12))
  expect_error(.cnv_normalize(matrix(0, 12, 20), 6, 10), "baseline")
  expect_equal(.cnv_eta(c(1, 1, 2, 2), c("A", "A", "B", "B")), 1)
  expect_true(is.na(.cnv_eta(rep(1, 4), c("A", "A", "B", "B"))))
})

test_that("loci are unique tuples and not parsed marker strings", {
  m <- data.frame(VARIANT_ID = c(3, 1, 2), CHROM = c("a", "a", "b"),
                  LOCUS = c("x", "x", "x"))
  expect_equal(.cnv_representatives(m, c("CHROM", "LOCUS"))$VARIANT_ID, 1:2)
  expect_error(.cnv_representatives(m, "ABSENT"), "metadata")
})

test_that("GDS CNV analysis reads depth and restores active filters", {
  skip_if_not_installed("SeqArray")
  vcf <- tempfile(fileext = ".vcf"); path <- tempfile(fileext = ".gds")
  samples <- paste0("s", 1:12)
  records <- vapply(1:20, function(i) paste(c(
    "1", i, paste0("m", i), "A", "G", ".", "PASS", ".", "GT:DP:AD",
    ifelse(i == 1 & 1:12 <= 3, "0/1:40:20,20", "0/1:20:10,10")
  ), collapse = "\t"), character(1))
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", samples), collapse = "\t"), records), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  gds <- SeqArray::seqOpen(path)
  on.exit({SeqArray::seqClose(gds); unlink(c(vcf, path))}, add = TRUE)
  before <- SeqArray::seqGetFilter(gds)
  meta <- data.frame(INDIVIDUALS = samples, STRATA = rep(c("A", "B"), each = 6),
                     PLATE = rep(c("P1", "P2"), each = 6))
  out <- detect_cnv(gds, meta, locus.columns = "VARIANT_ID", batch.columns = "PLATE",
                    min.loci = 10, min.samples = 6, write.files = FALSE,
                    verbose = FALSE, chunk.size = 7)
  expect_equal(out$statistics$STATUS[1], "CANDIDATE_DEPTH_VARIATION")
  expect_equal(out$statistics$N_HIGH[1], 3)
  expect_equal(out$statistics$HET_ALLELE_IMBALANCE, rep(0, 20))
  expect_equal(nrow(out$relative.depth), 240)
  expect_true(out$active.selection.restored)
  expect_identical(SeqArray::seqGetFilter(gds), before)
  expect_equal(nrow(out$stratum.statistics), 40)
  expect_true(all(out$stratum.statistics$NUMBER_CALLED == 6))
  expect_true(all(out$technical.risks$PARALOG_FLAG == "NOT_ASSESSED"))
  expect_true(all(out$technical.risks$READ_MAPPING_AUDIT == "NOT_ASSESSED"))
  expect_true(all(out$evidence$ALLELE_BALANCE == "AVAILABLE_BY_STRATUM"))
  expect_error(detect_cnv(gds, meta, locus.columns = "ABSENT",
                          write.files = FALSE, verbose = FALSE), "metadata")
  expect_identical(SeqArray::seqGetFilter(gds), before)
  folder <- tempfile("cnv-output-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = 1:18, sample.id = samples[1:10],
                         verbose = FALSE)
  filtered <- SeqArray::seqGetFilter(gds)
  saved <- detect_cnv(gds, meta, locus.columns = "VARIANT_ID",
    min.loci = 10, min.samples = 6, path.folder = folder, internal = TRUE,
    write.files = TRUE, verbose = FALSE)
  expect_equal(nrow(saved$statistics), 18)
  expect_equal(nrow(saved$samples), 10)
  expect_true(file.exists(file.path(folder, "cnv_diagnostic.png")))
  expect_true(file.exists(file.path(folder, "cnv_parameters.rds")))
  expect_identical(SeqArray::seqGetFilter(gds), filtered)
})

test_that("risk flags preserve unknown checks and match by locus, not row", {
  reps <- data.frame(LOCUS_INDEX = 1:3, LOCUS = c("a", "b", "c"))
  annotation <- data.frame(LOCUS = c("b", "a"),
    PARALOG_FLAG = c(FALSE, TRUE), REPEAT_OVERLAP = c(NA, FALSE),
    INDEPENDENT_SUPPORT = c(TRUE, NA))
  result <- .cnv_risk_annotations(reps, "LOCUS", annotation)
  expect_equal(result$PARALOG_FLAG, c("FLAGGED", "NOT_FLAGGED", "NOT_ASSESSED"))
  expect_equal(result$REPEAT_OVERLAP, c("NOT_FLAGGED", "NOT_ASSESSED", "NOT_ASSESSED"))
  expect_equal(result$INDEPENDENT_SUPPORT, c("NOT_ASSESSED", "SUPPORTED", "NOT_ASSESSED"))
  expect_true(all(result$LOW_MAPPABILITY == "NOT_ASSESSED"))
  expect_error(.cnv_risk_annotations(reps, "LOCUS", annotation[c(1, 1), ]),
               "one row")
  annotation$PARALOG_FLAG <- c(0, 1)
  expect_error(.cnv_risk_annotations(reps, "LOCUS", annotation), "logical")
})
