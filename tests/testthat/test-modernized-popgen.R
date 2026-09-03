make_popgen_gds <- function(readonly = TRUE) {
  sample.id <- paste0("sample", seq_len(12))
  records <- list(
    rep("0/1", 12),
    c(rep("0/0", 6), rep("1/1", 6)),
    c(rep("0/0", 4), rep("0/1", 4), rep("1/1", 4))
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
  list(
    gds = SeqArray::seqOpen(gds.file, readonly = readonly),
    files = c(vcf, gds.file), sample.id = sample.id,
    metadata = data.frame(
      INDIVIDUALS = sample.id,
      STRATA = rep(c("A", "B"), each = 6)
    )
  )
}

test_that("modern population diagnostics are GDS-native and restore filters", {
  skip_if_not_installed("SeqArray")
  fixture <- make_popgen_gds()
  on.exit({
    try(SeqArray::seqClose(fixture$gds), silent = TRUE)
    unlink(fixture$files)
  }, add = TRUE)
  SeqArray::seqSetFilter(fixture$gds, variant.id = 1:2, verbose = FALSE)
  before <- SeqArray::seqGetFilter(fixture$gds)

  biallelic <- detect_biallelic_problems(
    fixture$gds, write.files = FALSE, verbose = FALSE
  )
  private <- private_alleles(
    fixture$gds, fixture$metadata, write.files = FALSE, verbose = FALSE
  )
  diversity <- pi(
    fixture$gds, fixture$metadata, write.files = FALSE, verbose = FALSE
  )
  fh <- ibdg_fh(
    fixture$gds, fixture$metadata, write.files = FALSE, verbose = FALSE
  )
  het <- detect_het_outliers(
    fixture$gds, fixture$metadata, estimate.miscall = FALSE,
    write.files = FALSE, verbose = FALSE
  )
  after <- SeqArray::seqGetFilter(fixture$gds)

  expect_equal(nrow(biallelic$statistics), 2)
  expect_true(all(biallelic$statistics$RECORD_TYPE == "BIALLELIC"))
  expect_true(nrow(private$private.alleles) > 0)
  expect_true(all(c("A", "B", "OVERALL") %in%
                    diversity$population.statistics$GROUP))
  expect_equal(nrow(fh$individual.statistics), 12)
  expect_equal(nrow(het$statistics), 6)
  expect_true(all(vapply(
    list(biallelic, private, diversity, fh, het),
    function(x) isTRUE(x$active.selection.restored), logical(1)
  )))
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
})

test_that("FH uses count-scale observed and expected homozygotes", {
  fixture <- make_popgen_gds()
  on.exit({
    try(SeqArray::seqClose(fixture$gds), silent = TRUE)
    unlink(fixture$files)
  }, add = TRUE)
  result <- ibdg_fh(
    fixture$gds, fixture$metadata, write.files = FALSE, verbose = FALSE
  )$individual.statistics
  expect_true(all(result$OBSERVED_HOMOZYGOTES <= result$NUMBER_CALLED))
  expect_true(all(result$EXPECTED_HOMOZYGOTES_STRATUM <=
                    result$NUMBER_CALLED + 1e-8))
  expect_true(all(is.finite(result$FH_STRATUM)))
})

test_that("modern filters update the GDS audit metadata", {
  skip_if_not_installed("SeqArray")
  cases <- list(
    fis = function(gds, metadata, path) filter_fis(
      gds, metadata, fis.min.threshold = -0.2, fis.max.threshold = 0.2,
      min.call.rate = 0, verbose = FALSE, path.folder = path
    ),
    het = function(gds, metadata, path) filter_het(
      gds, metadata, marker.range = c(0, 0.8),
      min.marker.call.rate = 0, verbose = FALSE, path.folder = path
    ),
    hwe = function(gds, metadata, path) filter_hwe(
      gds, metadata, p.threshold = 0.1, adjustment = "none",
      min.samples = 2, min.call.rate = 0, verbose = FALSE,
      path.folder = path
    )
  )
  for (name in names(cases)) {
    fixture <- make_popgen_gds(readonly = FALSE)
    path <- tempfile(paste0(name, "-"))
    dir.create(path)
    suppressWarnings(cases[[name]](fixture$gds, fixture$metadata, path))
    marker <- genometranslator::extract_markers_metadata(
      fixture$gds, whitelist = FALSE
    )
    expect_true(any(marker$FILTERS == paste0("filter.", name)), info = name)
    expect_true(any(grepl("filter", list.files(path, recursive = TRUE))),
                info = name)
    SeqArray::seqClose(fixture$gds)
    unlink(c(fixture$files, path), recursive = TRUE)
  }
})

test_that("private haplotypes require and use phased multi-variant loci", {
  skip_if_not_installed("SeqArray")
  sample.id <- paste0("sample", 1:4)
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", sample.id
    ), collapse = "\t"),
    paste(c(
      "1", "100", "marker1", "A", "G", ".", "PASS", ".", "GT",
      "0|0", "0|0", "1|1", "1|1"
    ), collapse = "\t"),
    paste(c(
      "1", "110", "marker2", "C", "T", ".", "PASS", ".", "GT",
      "0|0", "0|0", "1|1", "1|1"
    ), collapse = "\t")
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds.file, readonly = FALSE)
  on.exit({
    try(SeqArray::seqClose(gds), silent = TRUE)
    unlink(c(vcf, gds.file))
  }, add = TRUE)
  marker <- genometranslator::extract_markers_metadata(
    gds, whitelist = FALSE
  )
  marker$LOCUS <- "locus1"
  invisible(genometranslator::update_genome_gds(
    gds = gds, node.name = "markers.meta", value = marker, sync = TRUE
  ))
  metadata <- data.frame(
    INDIVIDUALS = sample.id, STRATA = c("A", "A", "B", "B")
  )
  result <- private_haplotypes(
    gds, metadata, write.files = FALSE, verbose = FALSE
  )
  expect_equal(nrow(result$private.haplotypes), 2)
  expect_setequal(result$private.haplotypes$HAPLOTYPE, c("A|C", "G|T"))
  expect_true(result$active.selection.restored)
})

test_that("short-distance LD keeps one complete marker record per locus", {
  fixture <- make_popgen_gds(readonly = FALSE)
  on.exit({
    try(SeqArray::seqClose(fixture$gds), silent = TRUE)
    unlink(fixture$files)
  }, add = TRUE)
  marker <- genometranslator::extract_markers_metadata(
    fixture$gds, whitelist = FALSE
  )
  marker$LOCUS[1:2] <- "shared_locus"
  invisible(genometranslator::update_genome_gds(
    gds = fixture$gds, node.name = "markers.meta",
    value = marker, sync = TRUE
  ))
  path <- tempfile("ld-")
  dir.create(path)
  suppressWarnings(suppressMessages(filter_ld(
    fixture$gds, interactive.filter = FALSE,
    filter.short.ld = "first", filter.long.ld = NULL,
    verbose = FALSE, path.folder = path
  )))
  result <- genometranslator::extract_markers_metadata(
    fixture$gds, whitelist = FALSE
  )
  expect_equal(sum(result$FILTERS == "filter.short.ld"), 1)
  expect_equal(result$FILTERS[result$VARIANT_ID == 1], "whitelist")
})

test_that("filters open GDS filepaths with write access", {
  fixture <- make_popgen_gds()
  SeqArray::seqClose(fixture$gds)
  path <- tempfile("fis-path-")
  dir.create(path)
  result <- suppressWarnings(filter_fis(
    fixture$files[[2L]], fixture$metadata,
    fis.min.threshold = -0.2, fis.max.threshold = 0.2,
    min.call.rate = 0, verbose = FALSE, path.folder = path
  ))
  expect_s4_class(result, "SeqVarGDSClass")
  expect_true(any(
    genometranslator::extract_markers_metadata(
      result, whitelist = FALSE
    )$FILTERS == "filter.fis"
  ))
  SeqArray::seqClose(result)
  unlink(c(fixture$files, path), recursive = TRUE)
})
