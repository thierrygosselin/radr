test_that("rare support follows the observed minor allele", {
  dosage <- cbind(
    alt_supported = c(1, rep(0, 9)),
    ref_weak = c(1, rep(2, 9)),
    common = c(rep(1, 6), rep(0, 4))
  )
  ref <- cbind(
    alt_supported = c(10, rep(20, 9)),
    ref_weak = c(1, rep(0, 9)),
    common = c(rep(10, 6), rep(20, 4))
  )
  alt <- cbind(
    alt_supported = c(10, rep(0, 9)),
    ref_weak = c(19, rep(20, 9)),
    common = c(rep(10, 6), rep(0, 4))
  )
  result <- .rare_chunk(
    ids = 1:3, dosage = dosage, ref = ref, alt = alt,
    groups = rep(c("A", "B"), each = 5), by.strata = FALSE,
    max.mac = 3, min.call.rate = 0.8, min.samples = 5,
    min.completeness = 0.8, min.depth = 3, min.total.depth = 3,
    max.low.fraction = 0.5, min.balance = 0.2, min.hets = 1
  )
  expect_equal(result$MINOR_ALLELE, c("ALT", "REF", "ALT"))
  expect_equal(result$MINOR_ALLELE_COUNT, c(1, 1, 6))
  expect_equal(result$STATUS[1], "SUPPORTED_RARE")
  expect_equal(result$STATUS[2], "DEPTH_AND_BALANCE_WARNING")
  expect_equal(result$STATUS[3], "NOT_RARE")
})

test_that("rare support reads AD, stratifies, and restores filters", {
  skip_if_not_installed("SeqArray")
  sample.id <- paste0("sample", 1:20)
  records <- list(
    c("0/1:10,10", rep("0/0:20,0", 19)),
    c("0/1:19,1", rep("0/0:20,0", 19)),
    c(rep("0/1:10,10", 8), rep("0/0:20,0", 12))
  )
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", sample.id), collapse = "\t"),
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
    INDIVIDUALS = sample.id, STRATA = rep(c("A", "B"), each = 10)
  )
  parent <- tempfile("rare-support-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  result <- detect_rare_alleles(
    data = gds, strata = metadata, max.minor.allele.count = 3,
    min.call.rate = 0.8, min.samples.per.group = 5,
    min.depth.completeness = 0.8, min.minor.read.depth = 3,
    min.total.minor.depth = 3, chunk.size = 1, save.plots = FALSE,
    verbose = FALSE, path.folder = parent
  )
  after <- SeqArray::seqGetFilter(gds)
  expect_s3_class(result, "detect_rare_alleles")
  expect_true(result$active.selection.restored)
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(nrow(result$statistics), 2)
  expect_equal(nrow(result$strata.statistics), 6)
  expect_equal(result$statistics$STATUS,
               c("SUPPORTED_RARE", "DEPTH_AND_BALANCE_WARNING"))
  expect_equal(nrow(result$candidates), 1)
  expect_true(all(c(
    "rare_allele_marker_statistics.tsv",
    "rare_allele_stratum_statistics.tsv",
    "rare_allele_review_candidates.tsv", "rare_allele_thresholds.tsv"
  ) %in% list.files(result$path.folder)))
})

test_that("rare support controls reject invalid values", {
  expect_error(
    detect_rare_alleles("missing.gds", min.heterozygote.balance = 0.5),
    "smaller than 0.5"
  )
})
