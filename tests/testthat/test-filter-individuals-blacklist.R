test_that("sample blacklists apply by exact ID and preserve existing GDS filters", {
  output <- tempfile("individual-blacklist-")
  dir.create(output)
  previous.dir <- setwd(output)
  gds <- NULL
  on.exit({
    if (!is.null(gds)) SeqArray::seqClose(gds)
    setwd(previous.dir)
    unlink(output, recursive = TRUE)
  }, add = TRUE)

  ids <- c("001", "002", "003", "004", "005", "sample_06")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", ids), collapse = "\t"),
    paste(c("1", "100", "m1", "A", "G", ".", "PASS", ".", "GT",
            "0/0", "0/1", "1/1", "0/0", "0/1", "1/1"), collapse = "\t"),
    paste(c("1", "200", "m2", "A", "G", ".", "PASS", ".", "GT",
            "0/1", "1/1", "0/0", "0/1", "1/1", "0/0"), collapse = "\t"),
    paste(c("1", "300", "m3", "A", "G", ".", "PASS", ".", "GT",
            "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"), collapse = "\t")
  ), "fixture.vcf")
  SeqArray::seqVCF2GDS("fixture.vcf", "fixture.gds", verbose = FALSE)
  gds <- SeqArray::seqOpen("fixture.gds", readonly = FALSE)

  genometranslator::update_genome_gds(
    gds = gds, node.name = "individuals.meta",
    value = data.frame(
      INDIVIDUALS = ids,
      STRATA = rep("A", length(ids)),
      FILTERS = c("filter.individuals.missing", rep("whitelist", length(ids) - 1L))
    )
  )
  genometranslator::update_genome_gds(
    gds = gds, node.name = "markers.meta",
    value = data.frame(
      VARIANT_ID = 1:3,
      CHROM = "1",
      LOCUS = paste0("m", 1:3),
      POS = 100L * 1:3,
      FILTERS = c("whitelist", "whitelist", "filter.genomic.range")
    )
  )
  SeqArray::seqSetFilter(gds, sample.id = ids[-c(1L, 5L)], variant.id = 1L,
                         verbose = FALSE)
  if (utils::packageVersion("testthat") >= "3.1.7") {
    testthat::local_mocked_bindings(
      question = function(...) stop("A sample blacklist requested interactive input."),
      .package = "tgbase"
    )
    testthat::local_mocked_bindings(
      generate_stats = function(...) stop("A sample blacklist ran sample QC."),
      .package = "genometranslator"
    )
  }

  before <- genometranslator::extract_individuals_metadata(gds)
  expect_error(filter_individuals(
    gds, blacklist.id = c("002", "unknown_id"), parallel.core = 1L, verbose = FALSE
  ), "unknown_id")
  writeLines(c("WRONG_COLUMN", "002"), "malformed.tsv")
  expect_error(filter_individuals(
    gds, blacklist.id = "malformed.tsv", parallel.core = 1L, verbose = FALSE
  ), "INDIVIDUALS")
  expect_identical(genometranslator::extract_individuals_metadata(gds), before)

  expect_identical(filter_individuals(
    gds, blacklist.id = c("001", "002", "sample_06", "002"), parallel.core = 1L, verbose = FALSE
  ), gds)
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), c("003", "004"))
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), 1L)

  expect_error(filter_individuals(
    gds, blacklist.id = c("003", "004"), parallel.core = 1L, verbose = FALSE
  ), "every active sample")
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), c("003", "004"))

  filter_individuals(gds, blacklist.id = "002", parallel.core = 1L, verbose = FALSE)
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), c("003", "004"))

  writeLines(c("INDIVIDUALS", "003"), "blacklist.tsv")
  filter_individuals(
    gds, interactive.filter = FALSE, blacklist.id = "blacklist.tsv", parallel.core = 1L, verbose = FALSE
  )
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), "004")
  expect_length(SeqArray::seqGetData(gds, "variant.id"), 0L)
  expect_identical(genometranslator::extract_markers_metadata(gds)$FILTERS,
                   c("filter.monomorphic", "whitelist", "filter.genomic.range"))

  reports <- list.files(
    pattern = "^blacklist\\.individuals\\.by_id\\.tsv$", recursive = TRUE
  )
  expect_length(reports, 2L)
  expect_setequal(
    unlist(lapply(reports, function(path) {
      readr::read_tsv(path, col_types = readr::cols(.default = "c"))$INDIVIDUALS
    })),
    c("002", "003", "sample_06")
  )

  SeqArray::seqClose(gds)
  gds <- genometranslator::read_genome("fixture.gds", verbose = FALSE)
  excluded <- genometranslator::extract_individuals_metadata(gds, blacklist = TRUE)
  expect_identical(excluded$INDIVIDUALS, c("001", "002", "003", "sample_06"))
  expect_identical(
    excluded$FILTERS,
    c("filter.individuals.missing", rep("filter.individuals.blacklist", 3L))
  )
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), c("004", "005"))
  expect_identical(SeqArray::seqGetData(gds, "variant.id"), 2L)
})
