# Integration smoke test for the GDS-only inversion scan.

set.seed(701)
n.samples <- 30L
n.variants <- 1200L
sample.id <- paste0("sample", seq_len(n.samples))
chromosome <- rep(c("1", "2"), each = n.variants / 2L)
position <- ave(seq_len(n.variants), chromosome, FUN = seq_along) * 1000L

dosage <- matrix(
  stats::rbinom(n.samples * n.variants, size = 2L, prob = 0.35),
  nrow = n.samples,
  ncol = n.variants
)

# Two adjacent Chr2 windows carry three strongly differentiated haplotypes.
inversion.index <- 801:1000
arrangement <- rep(0:2, each = n.samples / 3L)
for (j in inversion.index) {
  error <- stats::rbinom(n.samples, size = 1L, prob = 0.03)
  dosage[, j] <- pmin(2L, pmax(0L, arrangement + ifelse(error == 1L, 1L, 0L)))
}
dosage[sample(
  length(dosage), size = floor(length(dosage) * 0.01)
)] <- NA_integer_

vcf <- tempfile(fileext = ".vcf")
gds.file <- tempfile(fileext = ".gds")
on.exit(unlink(c(vcf, gds.file)), add = TRUE)

header <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", sample.id), collapse = "\t")
)
gt <- c("0/0", "0/1", "1/1")
records <- vapply(seq_len(n.variants), function(j) {
  calls <- ifelse(is.na(dosage[, j]), "./.", gt[dosage[, j] + 1L])
  paste(c(
    chromosome[j], position[j], paste0("v", j), "A", "C", ".", "PASS", ".",
    "GT", calls
  ), collapse = "\t")
}, character(1))
writeLines(c(header, records), vcf)

SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  library(radr)
} else {
  load_all <- get("load_all", envir = asNamespace("pkgload"))
  load_all(".", quiet = TRUE)
}

result.parent <- tempfile("detect-inversions-")
dir.create(result.parent)
on.exit(unlink(result.parent, recursive = TRUE), add = TRUE)

sample.metadata <- data.frame(
  INDIVIDUALS = sample.id,
  STRATA = rep(c("batch_A", "batch_B", "batch_C"), each = n.samples / 3L)
)

result <- suppressWarnings(detect_inversions(
  data = gds.file,
  strata = sample.metadata,
  window.snps = 100L,
  sensitivity.window.snps = c(100L, 200L),
  outlier.quantile = 0.80,
  min.candidate.windows = 1L,
  known.regions = data.frame(
    chromosome = c("2", "1"),
    start = c(250000, 1000),
    end = c(350000, 50000),
    type = c("putative_centromere", "assembly_gap")
  ),
  return.ld = TRUE,
  save.plots = TRUE,
  plot.formats = "png",
  verbose = FALSE,
  path.folder = result.parent
))

stopifnot(inherits(result, "detect_inversions"))
stopifnot(nrow(result$windows) == 12L)
stopifnot(all(result$windows$n_input_snps == 100L))
stopifnot(sum(
  result$windows$candidate_window & result$windows$chromosome == "2" &
    result$windows$start >= 201000 & result$windows$end <= 400000
) == 2L)
stopifnot(nrow(result$candidates) >= 1L)
stopifnot(all(c(
  "size_bp", "cluster_separation", "cluster_compactness",
  "middle_heterozygosity_excess", "homokaryotype_mean_ld_r2",
  "flanking_mean_ld_r2", "boundary_contrast", "internal_transition_max",
  "known_region_overlap", "evidence_score", "evidence_strength",
  "candidate_class", "alternative_explanations"
) %in% names(result$candidates)))
candidate <- result$candidates[result$candidates$chromosome == "2", ][1, ]
stopifnot(candidate$three_cluster_evidence)
stopifnot(candidate$smallest_cluster_n == 10L)
stopifnot(candidate$cluster_separation > 1)
stopifnot(candidate$middle_heterozygosity_excess > 0)
stopifnot(candidate$regional_mean_ld_r2 > candidate$homokaryotype_mean_ld_r2)
stopifnot(candidate$regional_mean_ld_r2 > candidate$flanking_mean_ld_r2)
stopifnot(candidate$boundary_contrast > 0)
stopifnot(candidate$known_region_overlap == "putative_centromere")
stopifnot(candidate$n_known_region_overlaps == 1L)
stopifnot(candidate$evidence_strength == "strong")
stopifnot(nrow(result$sensitivity) > 0L)
stopifnot(all(c(
  "candidate_id", "individual", "arrangement", "arrangement_dosage",
  "arrangement_confidence"
) %in% names(result$arrangement.genotypes)))
stopifnot(all(as.character(result$arrangement.genotypes$arrangement) %in%
  c("AA", "AB", "BB")))
stopifnot(all(result$arrangement.genotypes$arrangement_confidence >= 0.5))
stopifnot(nrow(result$homokaryotype.whitelist) > 0L)
stopifnot(dir.exists(result$path.folder))
stopifnot(file.exists(file.path(result$path.folder, "inversion_windows.tsv")))
stopifnot(file.exists(file.path(
  result$path.folder, "inversion_arrangement_genotypes.tsv"
)))
stopifnot(file.exists(file.path(
  result$path.folder, "homokaryotypes_all_candidates_whitelist.tsv"
)))
stopifnot(file.exists(file.path(result$path.folder, "window_scores.png")))
stopifnot(any(grepl("_ld[.]png$", result$output.files$files)))
stopifnot(any(grepl(
  "_metadata_by_arrangement[.]tsv$", result$output.files$files
)))
stopifnot(any(grepl(
  "_pca_metadata_STRATA[.]png$", result$output.files$files
)))
stopifnot(nrow(result$diagnostics[[1]]$metadata_audit) == 1L)
stopifnot(nrow(result$diagnostics[[1]]$metadata_contingency) > 0L)
stopifnot(all(vapply(result$diagnostics, function(x) {
  all(c("scores", "cluster_summary", "mean_ld_r2", "ld_summary",
        "ld_matrices", "cluster_separation", "cluster_compactness",
        "middle_heterozygosity_excess", "homokaryotype_mean_ld_r2") %in%
      names(x))
}, logical(1))))
