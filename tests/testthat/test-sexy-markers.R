test_that("sexy_markers detects directional signals and restores filters", {
  skip_if_not_installed("SeqArray")

  sample.id <- c(paste0("F", 1:6), paste0("M", 1:6))
  female.missing <- c(rep("./.", 6), rep("0/1", 6))
  male.missing <- c(rep("0/1", 6), rep("./.", 6))
  x.het <- c(rep("0/1", 6), rep("0/0", 6))
  z.het <- c(rep("0/0", 6), rep("0/1", 6))
  background <- rep(c("0/0", "0/1"), 6)
  records <- list(female.missing, male.missing, x.het, z.het, background)

  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", sample.id), collapse = "\t"),
    vapply(seq_along(records), function(i) paste(c(
      "1", i * 100, paste0("marker", i), "A", "G", ".", "PASS", ".",
      "GT", records[[i]]
    ), collapse = "\t"), character(1))
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds.file)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  SeqArray::seqSetFilter(gds, variant.id = 1:4, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  strata <- data.frame(
    INDIVIDUALS = sample.id,
    SEX = rep(c("F", "M"), each = 6),
    BATCH = rep(c("A", "B"), 6)
  )
  result <- sexy_markers(
    data = gds,
    strata = strata,
    sex.column = "SEX",
    min.samples.per.sex = 5,
    save.plots = FALSE,
    verbose = FALSE,
    path.folder = tempdir()
  )
  after <- SeqArray::seqGetFilter(gds)

  expect_s3_class(result, "sexy_markers")
  expect_true(result$active_selection_restored)
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(nrow(result$statistics), 4)
  expect_true(result$statistics$candidate_y_like[1])
  expect_true(result$statistics$candidate_w_like[2])
  expect_true(result$statistics$candidate_x_heterozygosity[3])
  expect_true(result$statistics$candidate_z_heterozygosity[4])
  expect_false(result$depth_available)
  expect_true(all(c(
    "sex_marker_statistics.tsv", "candidate_sex_markers.tsv",
    "sex_sample_summary.tsv", "sex_metadata_audit.tsv"
  ) %in% list.files(result$path.folder)))
})

test_that("sex labels and numeric arguments are validated", {
  expect_equal(
    .sex_normalise(c("female", "F", "MALE", "m", "unknown", NA)),
    c("F", "F", "M", "M", "U", "U")
  )
  expect_error(
    .sex_check_probability(1.2, "threshold"),
    "zero to one"
  )
  expect_error(
    .sex_check_count(1.5, "chunk.size", 1L),
    "whole number"
  )
})
