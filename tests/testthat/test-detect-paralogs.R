test_that("HD statistics use called genotypes and classify evidence", {
  dosage <- cbind(
    candidate = c(rep(1, 8), NA, NA),
    balanced = c(rep(1, 4), rep(0, 6)),
    insufficient = c(1, rep(0, 9))
  )
  ref <- cbind(
    candidate = c(rep(18, 8), NA, NA),
    balanced = c(rep(10, 4), rep(0, 6)),
    insufficient = c(10, rep(0, 9))
  )
  alt <- cbind(
    candidate = c(rep(2, 8), NA, NA),
    balanced = c(rep(10, 4), rep(0, 6)),
    insufficient = c(10, rep(0, 9))
  )
  result <- .paralog_chunk_statistics(
    variant.id = 1:3,
    dosage = dosage,
    ref.depth = ref,
    alt.depth = alt,
    groups = rep(c("A", "B"), each = 5),
    by.strata = FALSE,
    min.call.rate = 0.7,
    min.samples.per.group = 5,
    min.heterozygotes = 3,
    min.heterozygote.depth = 20,
    heterozygosity.threshold = 0.6,
    deviation.threshold = 4
  )

  expect_equal(result$H_ALL[1], 0.8)
  expect_equal(result$H_CALLED[1], 1)
  expect_equal(result$STATUS[1], "PARALOG_CANDIDATE")
  expect_equal(result$STATUS[2], "SINGLETON_LIKE")
  expect_equal(result$STATUS[3], "INSUFFICIENT_DATA")
  expect_match(
    result$INSUFFICIENT_REASON[3],
    "too_few_heterozygotes_with_depth"
  )
})

test_that("detect_paralogs reads AD, stratifies, and restores filters", {
  skip_if_not_installed("SeqArray")

  sample.id <- paste0("sample", 1:20)
  candidate <- rep("0/1:18,2", 20)
  balanced <- c(rep("0/1:10,10", 8), rep("0/0:20,0", 12))
  low.call <- c(rep("0/1:10,10", 2), rep("./.:.,.", 18))
  records <- list(candidate, balanced, low.call)
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", sample.id
    ), collapse = "\t"),
    vapply(seq_along(records), function(i) paste(c(
      "1", i * 100, paste0("marker", i), "A", "G", ".", "PASS", ".",
      "GT:AD", records[[i]]
    ), collapse = "\t"), character(1))
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds.file)
  on.exit({
    try(SeqArray::seqClose(gds), silent = TRUE)
    unlink(c(vcf, gds.file))
  }, add = TRUE)

  SeqArray::seqSetFilter(gds, variant.id = 1:2, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  metadata <- data.frame(
    INDIVIDUALS = sample.id,
    STRATA = rep(c("A", "B"), each = 10)
  )
  parent <- tempfile("paralogs-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)

  result <- detect_paralogs(
    data = gds,
    strata = metadata,
    min.call.rate = 0.7,
    min.samples.per.group = 5,
    min.heterozygotes = 3,
    min.heterozygote.depth = 20,
    heterozygosity.threshold = 0.6,
    deviation.threshold = 4,
    chunk.size = 1,
    save.plots = FALSE,
    verbose = FALSE,
    path.folder = parent
  )
  after <- SeqArray::seqGetFilter(gds)

  expect_s3_class(result, "detect_paralogs")
  expect_true(result$active.selection.restored)
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(nrow(result$statistics), 2)
  expect_equal(nrow(result$strata.statistics), 6)
  expect_equal(result$statistics$STATUS[1], "PARALOG_CANDIDATE")
  expect_true(result$statistics$ANY_PARALOG_CANDIDATE[1])
  expect_equal(result$statistics$H_CALLED[1], 1)
  expect_equal(nrow(result$candidates), 1)
  expect_true(all(c(
    "paralog_marker_statistics.tsv",
    "paralog_stratum_statistics.tsv",
    "candidate_paralog_markers.tsv",
    "paralog_thresholds.tsv"
  ) %in% list.files(result$path.folder)))
})

test_that("detect_paralogs controls validate clearly", {
  expect_error(
    .paralog_check_probability(1.1, "threshold"),
    "between zero and one"
  )
  expect_error(
    .paralog_check_count(1.5, "chunk.size", 1),
    "whole number"
  )
  expect_error(
    .paralog_check_nonnegative(-1, "depth"),
    "non-negative"
  )
})
