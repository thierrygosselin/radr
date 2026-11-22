test_that("combined lists persist, audit both axes and reset to the QC baseline", {
  previous <- getwd()
  folder <- tempfile("list-filter-")
  dir.create(folder)
  setwd(folder)
  gds <- NULL
  on.exit({
    if (!is.null(gds)) SeqArray::seqClose(gds)
    setwd(previous)
    unlink(folder, recursive = TRUE)
  }, add = TRUE)

  file.copy(SeqArray::seqExampleFileName("gds"), "fixture.gds")
  gds <- SeqArray::seqOpen("fixture.gds", readonly = FALSE)
  samples <- SeqArray::seqGetData(gds, "sample.id")
  variants <- SeqArray::seqGetData(gds, "variant.id")
  markers <- sprintf("%08d", variants)
  genometranslator::update_genome_gds(gds, node.name = "individuals.meta", value =
    tibble::tibble(INDIVIDUALS = samples, STRATA = "A",
                   FILTERS = ifelse(samples == samples[1], "prior.qc", "whitelist")))
  genometranslator::update_genome_gds(gds, node.name = "markers.meta", value =
    tibble::tibble(VARIANT_ID = variants, MARKERS = markers, CHROM = "1", LOCUS = markers,
                   FILTERS = ifelse(variants == variants[1], "prior.qc", "whitelist")))
  SeqArray::seqClose(gds)
  gds <- NULL

  expect_identical(filter_lists("fixture.gds"), "fixture.gds")
  expect_error(filter_lists(), "Supply")
  expect_error(filter_lists("fixture.gds", blacklist.individuals = samples[2],
                            whitelist.markers = "missing_name"), "Unknown")
  expect_error(filter_lists("fixture.gds", whitelist.markers = character()), "No active")
  expect_length(list.dirs(recursive = FALSE), 0L)
  writeLines(c("MARKERS", markers[3], markers[3]), "inversion.tsv")
  gds <- filter_lists("fixture.gds", samples[1:5], samples[4], variants[1:5], "inversion.tsv")
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), samples[c(2, 3, 5)])
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[c(2, 4, 5)])
  manifest <- list.files(pattern = "filter_audit_manifest.tsv", recursive = TRUE)
  expect_length(manifest, 1L)
  audit <- readr::read_tsv(manifest, show_col_types = FALSE)
  expect_equal(audit$INDIVIDUALS_KEPT, 3)
  expect_equal(audit$MARKERS_KEPT, 3)
  SeqArray::seqClose(gds)
  gds <- genometranslator::read_genome("fixture.gds", verbose = FALSE)
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[c(2, 4, 5)])
  expect_equal(genometranslator::list_filters(gds, history = TRUE, verbose = FALSE)$history$FILTERS,
               "filter_lists")
  genometranslator::reset_filters(gds, filter.lists = TRUE)
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), samples[-1])
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[-1])

  # A temporary selection on the connection must survive applying and resetting a list.
  SeqArray::seqSetFilter(gds, sample.id = samples[2:4], variant.id = variants[2:4], verbose = FALSE)
  writeLines(c("variant_id", as.character(variants[2])), "loadings.tsv")
  expect_error(filter_lists(gds, blacklist.markers = "loadings.tsv"), "VARIANT_ID or MARKERS")
  writeLines(c("individual", samples[2:4]), "homokaryotypes.tsv")
  expect_error(filter_lists(gds, whitelist.individuals = "homokaryotypes.tsv"), "INDIVIDUALS")
  writeLines(c("VARIANT_ID", as.character(variants[2])), "loadings.tsv")
  writeLines(c("INDIVIDUALS", samples[2:4]), "homokaryotypes.tsv")
  filter_lists(gds, whitelist.individuals = "homokaryotypes.tsv", blacklist.markers = "loadings.tsv")
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[3:4])
  genometranslator::reset_filters(gds, filter.lists = TRUE)
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[2:4])
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), samples[2:4])
  # An audit failure occurs after both metadata updates and the history append.
  # The original selections and the cumulative history must still be recoverable.
  if (utils::packageVersion("testthat") >= "3.1.7") {
    before <- genometranslator::list_filters(gds, history = TRUE, verbose = FALSE)
    testthat::local_mocked_bindings(
      .write_filter_audit = function(...) stop("simulated audit failure"), .package = "radr"
    )
    expect_error(filter_lists(gds, blacklist.individuals = samples[2],
                              blacklist.markers = variants[2]), "simulated audit failure")
    expect_identical(SeqArray::seqGetData(gds, "variant.id"), variants[2:4])
    expect_identical(SeqArray::seqGetData(gds, "sample.id"), samples[2:4])
    expect_equal(genometranslator::list_filters(gds, history = TRUE, verbose = FALSE), before)
  }
})
