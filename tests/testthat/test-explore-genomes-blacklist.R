test_that("explore_genomes applies the ID blacklist before its first filter", {
  if (utils::packageVersion("testthat") < "3.1.7") {
    skip("This regression requires local_mocked_bindings.")
  }
  output <- tempfile("explore-blacklist-")
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
            "0/0", "0/1", "./.", "0/0", "0/0", "0/0"), collapse = "\t"),
    paste(c("1", "200", "m2", "A", "G", ".", "PASS", ".", "GT",
            "0/1", "1/1", "./.", "0/1", "1/1", "0/0"), collapse = "\t"),
    paste(c("1", "300", "m3", "A", "G", ".", "PASS", ".", "GT",
            "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"), collapse = "\t")
  ), "fixture.vcf")
  SeqArray::seqVCF2GDS("fixture.vcf", "fixture.gds", verbose = FALSE)
  gds <- SeqArray::seqOpen("fixture.gds", readonly = FALSE)
  genometranslator::genome_gds_skeleton(gds)
  genometranslator::update_genome_gds(gds, node.name = "data.source", value = "dart")

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
      MARKERS = paste0("m", 1:3),
      POS = 100L * 1:3,
      FILTERS = c("whitelist", "whitelist", "filter.genomic.range")
    )
  )

  writeLines(c("INDIVIDUALS", "001", "002"), "blacklist.tsv")
  expected.samples <- ids[-(1:2)]
  SeqArray::seqClose(gds)
  gds <- NULL
  expect_true(file.copy("fixture.gds", "interactive.gds"))

  calls <- logical()
  testthat::local_mocked_bindings(
    filter_dart_reproducibility = function(data, interactive.filter, ...) {
      calls <<- c(calls, interactive.filter)
      expect_identical(SeqArray::seqGetData(data, "sample.id"), expected.samples)
      expect_identical(SeqArray::seqGetData(data, "variant.id"), 2L)
      excluded <- genometranslator::extract_individuals_metadata(data, blacklist = TRUE)
      expect_identical(excluded$INDIVIDUALS, ids[1:2])
      expect_identical(
        excluded$FILTERS, c("filter.individuals.missing", "filter.individuals.blacklist")
      )
      stop("Reached the DArT reproducibility filter")
    },
    .package = "radr"
  )

  gds <- genometranslator::read_genome("fixture.gds", verbose = FALSE)
  expect_error(
    explore_genomes(
      data = gds, blacklist.id = c("001", "002"),
      interactive.filter = FALSE, parallel.core = 1L, verbose = FALSE
    ),
    "Reached the DArT reproducibility filter"
  )
  expect_error(
    explore_genomes(
      data = "interactive.gds", blacklist.id = "blacklist.tsv",
      interactive.filter = TRUE, parallel.core = 1L, verbose = FALSE
    ),
    "Reached the DArT reproducibility filter"
  )
  expect_identical(calls, c(FALSE, TRUE))
})
