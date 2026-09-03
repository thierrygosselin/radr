test_that("native scanner finds and classifies perfect tandem repeats", {
  markers <- data.frame(
    VARIANT_ID = 1:4,
    MARKERS = paste0("marker", 1:4),
    SEQUENCE = c(
      "GGATATATATATATCC",
      "ACGACGACGACGACGTT",
      "AAAAAAAAAA",
      "ACGTAGCTAGCATCG"
    )
  )
  hits <- .ssr_scan_sequences(
    markers, lengths = 1:3, repeats = c(10, 6, 5),
    canonicalize = TRUE
  )

  expect_equal(hits$MARKERS, c("marker1", "marker2", "marker3"))
  expect_equal(hits$MOTIF_LENGTH, c(2, 3, 1))
  expect_equal(hits$REPEAT_COUNT, c(6, 5, 10))
  expect_equal(hits$START, c(3, 1, 1))
  expect_equal(hits$CANONICAL_MOTIF, c("AT", "ACG", "A"))
})

test_that("compound motifs are assigned to their primitive repeat", {
  markers <- data.frame(
    VARIANT_ID = 1L, MARKERS = "marker1",
    SEQUENCE = "ATATATATATAT"
  )
  hits <- .ssr_scan_sequences(
    markers, lengths = c(2, 4), repeats = c(6, 3),
    canonicalize = TRUE
  )
  expect_equal(nrow(hits), 1)
  expect_equal(hits$MOTIF, "AT")
  expect_equal(hits$PRIMITIVE_MOTIF, "AT")
})

test_that("detect_microsatellites uses active GDS markers", {
  skip_if_not_installed("SeqArray")
  sample.id <- c("sample1", "sample2")
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", sample.id), collapse = "\t"),
    "1\t100\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
    "1\t200\tm2\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
    "1\t300\tm3\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(
    gds.file, readonly = FALSE, allow.duplicate = TRUE
  )
  on.exit({
    try(SeqArray::seqClose(gds), silent = TRUE)
    unlink(c(vcf, gds.file))
  }, add = TRUE)
  genometranslator::update_genome_gds(
    gds = gds, node.name = "markers.meta",
    value = data.frame(
      VARIANT_ID = 1:3, MARKERS = paste0("marker", 1:3),
      CHROM = "1", LOCUS = 1:3, POS = c(100, 200, 300),
      SEQUENCE = c("ATATATATATAT", "ACGTACGT", NA_character_)
    ),
    sync = FALSE
  )
  SeqArray::seqSetFilter(gds, variant.id = 1:2, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  parent <- tempfile("ssr-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)

  result <- detect_microsatellites(
    data = gds, motif.lengths = 2, min.repeats = 6,
    save.plots = FALSE, verbose = FALSE, path.folder = parent
  )
  after <- SeqArray::seqGetFilter(gds)

  expect_s3_class(result, "detect_microsatellites")
  expect_true(result$active.selection.restored)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(nrow(result$marker.status), 2)
  expect_equal(result$marker.status$SSR_STATUS,
               c("SSR_DETECTED", "NO_SSR_DETECTED"))
  expect_equal(nrow(result$candidates), 1)
  expect_equal(nrow(result$without.detected.repeats), 1)
  expect_equal(nrow(result$sequence.unavailable), 0)
  expect_true(all(c(
    "microsatellite_hits.tsv", "microsatellite_marker_status.tsv",
    "candidate_microsatellite_markers.tsv",
    "markers_without_detected_microsatellites.tsv",
    "markers_without_sequence.tsv"
  ) %in% list.files(result$path.folder)))
})

test_that("microsatellite thresholds validate clearly", {
  expect_error(.ssr_validate_lengths(c(2, 2)), "unique positive integers")
  expect_error(.ssr_validate_repeats(c(4, 5), 1:3), "one.*per motif")
  expect_equal(unname(.ssr_validate_repeats(c(`2` = 6, `3` = 5), 3:2)),
               c(5L, 6L))
})
