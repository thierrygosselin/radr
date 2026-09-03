test_that("summary statistics are correct and filters are restored", {
  skip_if_not_installed("SeqArray")

  sample.id <- paste0("sample", 1:8)
  records <- list(
    c("0/0", "0/0", "0/1", "./.", "1/1", "1/1", "0/1", "0/1"),
    c("0/0", "0/0", "0/0", "0/1", "0/0", "0/0", "0/0", "0/1"),
    rep("0/0", 8)
  )
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", sample.id
    ), collapse = "\t"),
    vapply(seq_along(records), function(i) paste(c(
      "1", i * 100, paste0("marker", i), "A", "G", ".", "PASS", ".",
      "GT", records[[i]]
    ), collapse = "\t"), character(1))
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds.file)
  on.exit({
    try(SeqArray::seqClose(gds), silent = TRUE)
    unlink(c(vcf, gds.file))
  }, add = TRUE)

  SeqArray::seqSetFilter(gds, variant.id = 1:3, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  metadata <- data.frame(
    INDIVIDUALS = sample.id,
    STRATA = rep(c("A", "B"), each = 4)
  )
  output <- tempfile("genomic-summary-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)

  result <- summarise_genomic_data(
    data = gds, strata = metadata, chunk.size = 1L,
    verbose = FALSE, path.folder = output
  )
  after <- SeqArray::seqGetFilter(gds)

  expect_s3_class(result, "summarise_genomic_data")
  expect_true(result$active.selection.restored)
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(nrow(result$marker.statistics), 3)
  expect_equal(nrow(result$stratum.statistics), 9)

  pooled1 <- result$stratum.statistics |>
    dplyr::filter(.data$VARIANT_ID == 1, .data$GROUP == "OVERALL")
  expect_equal(pooled1$NUMBER_CALLED, 7)
  expect_equal(pooled1$CALL_RATE, 7 / 8, tolerance = 1e-6)
  expect_equal(pooled1$ALT_FREQUENCY, 0.5)
  expect_equal(pooled1$MINOR_ALLELE, "ALT")

  group.a <- result$stratum.statistics |>
    dplyr::filter(.data$VARIANT_ID == 1, .data$GROUP == "A")
  expect_equal(group.a$ALT_FREQUENCY, 1 / 6, tolerance = 1e-5)
  expect_equal(group.a$POOLED_MINOR_ALLELE_FREQUENCY, 1 / 6,
               tolerance = 1e-5)
  monomorphic <- result$marker.statistics |>
    dplyr::filter(.data$VARIANT_ID == 3)
  expect_true(is.na(monomorphic$FIS))
  expect_true(all(c(
    "genomic_marker_statistics.tsv",
    "genomic_stratum_statistics.tsv",
    "genomic_group_summary.tsv"
  ) %in% list.files(result$path.folder)))
})

test_that("summary helper distinguishes ALT frequency from MAF", {
  dosage <- matrix(c(2, 2, 2, 1), ncol = 1)
  result <- .summary_chunk_statistics(
    variant.id = 10, dosage = dosage,
    groups = rep("A", 4), by.strata = FALSE
  )
  expect_equal(result$ALT_FREQUENCY, 7 / 8)
  expect_equal(result$MINOR_ALLELE_FREQUENCY, 1 / 8)
  expect_equal(result$MINOR_ALLELE, "REF")
})
