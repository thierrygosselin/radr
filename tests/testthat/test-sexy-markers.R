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
  expect_equal(nrow(result$assignment_panel), 2)
  expect_equal(
    result$assignment_panel$ASSIGNMENT_DIRECTION,
    c("Y-like", "W-like")
  )
  expect_true(all(c(
    "sex_marker_statistics.tsv", "candidate_sex_markers.tsv",
    "sex_assignment_panel.tsv", "sex_sample_summary.tsv",
    "sex_metadata_audit.tsv"
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

test_that("SilicoDArT markers use presence and not heterozygosity", {
  dosage <- cbind(
    y_like = c(rep(0, 6), rep(2, 6)),
    w_like = c(rep(2, 6), rep(0, 6))
  )
  female <- rep(c(TRUE, FALSE), each = 6)
  male <- !female
  statistics <- .sex_chunk_statistics(
    dosage = dosage,
    depth = NULL,
    depth.scale = rep(1, 12),
    female = female,
    male = male,
    coverage.threshold = 1,
    dominant = c(TRUE, TRUE)
  )

  expect_equal(unname(statistics$female_presence), c(0, 1))
  expect_equal(unname(statistics$male_presence), c(1, 0))
  expect_equal(
    statistics$presence_source,
    rep("silicodart_presence", 2)
  )
  expect_true(all(is.na(statistics$female_heterozygosity)))
  expect_true(all(is.na(statistics$male_heterozygosity)))
  expect_true(all(is.na(statistics$heterozygosity_p)))

  with.missing <- dosage[, 2, drop = FALSE]
  with.missing[1, 1] <- NA_real_
  missing.statistics <- .sex_chunk_statistics(
    dosage = with.missing,
    depth = NULL,
    depth.scale = rep(1, 12),
    female = female,
    male = male,
    coverage.threshold = 1,
    dominant = TRUE
  )
  expect_equal(unname(missing.statistics$female_presence), 1)
})

test_that("metadata audit excludes identifier fields", {
  metadata <- data.frame(
    INDIVIDUALS = paste0("I", 1:12),
    TARGET_ID = paste0("T", 1:12),
    ID_SEQ = 1:12,
    STRATA_SEQ = rep(1:2, each = 6),
    FILTERS = "whitelist",
    SEX = rep(c("F", "M"), each = 6),
    SEX_RADR = rep(c("F", "M"), each = 6),
    PLATE = rep(c("A", "B"), each = 6)
  )
  audit <- .sex_metadata_audit(metadata, "SEX")
  expect_identical(audit$variable, "PLATE")
})
