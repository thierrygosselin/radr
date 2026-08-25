#' Detect candidate inversion-associated genomic regions
#'
#' Scan a diploid, biallelic GDS for genomic windows whose local population
#' structure differs from the genomic background. The implementation follows
#' the local-PCA principle: each window is represented by a low-rank covariance
#' matrix among individuals, distances are calculated between those matrices,
#' and classical multidimensional scaling (MDS) is used to identify unusual
#' windows. Candidate regions are then summarised using regional PCA,
#' heterozygosity, and linkage disequilibrium (LD).
#'
#' This is a screening method. A candidate region is not proof of a physical
#' inversion, and the returned coordinates describe an inversion-associated
#' haploblock rather than validated breakpoints. Long-read, read-pair, split-read,
#' cytogenetic, or genetic-map evidence is needed for breakpoint confirmation.
#'
#' @section Chromosome-specific windows:
#' Windows are constructed independently within each chromosome, linkage group,
#' or scaffold. A window can therefore never contain markers from two linkage
#' groups. Local PCA is not performed once on an entire linkage group: it is
#' performed separately for every SNP window within that linkage group. The
#' resulting window-level covariance summaries are compared across all valid
#' windows in the requested scan. Use `chromosome` to restrict that comparison
#' to one or more linkage groups. Standard genomic-position figures are faceted
#' by chromosome or linkage group in the results folder.
#'
#' @section Missing genotypes and RADseq data:
#' SNPs below `min.call.rate` are excluded separately within each window. For the
#' remaining SNPs, each missing genotype is replaced temporarily by that SNP's
#' mean dosage before covariance PCA is calculated. Mean imputation is
#' equivalent to giving a missing sample the observed allele-frequency
#' expectation at that SNP: after centring, it contributes zero rather than an
#' invented homozygous or heterozygous deviation. The GDS is not modified, and
#' observed calls are retained for heterozygosity calculations.
#'
#' This simple imputation keeps all samples in a common PCA space, but it does not
#' correct non-random missingness. In RADseq and other reduced-representation
#' data, allele dropout at restriction-site polymorphisms, uneven depth, library
#' quality, lanes, plates, populations, marker panels, or alignment quality may
#' correlate with biological groups. Mean imputation can then shrink affected
#' samples toward the window centre, alter covariance, weaken a real signal, or
#' create a batch-associated local signal. Before interpreting candidates, plot
#' call rate and depth against chromosome position, compare PCA clusters with
#' batch and plate metadata, repeat the scan at stricter `min.call.rate` values,
#' and confirm that candidate windows remain after removing problematic samples
#' or markers. Imputation makes the matrix computable; it does not make missing
#' data unbiased.
#'
#' LD uses observed genotypes with pairwise-complete correlations; missing LD
#' genotypes are not mean-imputed.
#'
#' @section Output and plotting:
#' Following other `radr` detection functions, each call creates a dated
#' `detect_inversions` results folder in the working directory (or below the
#' parent supplied with `path.folder`). It records the function arguments,
#' window and candidate tables, individual PCA scores, cluster summaries, LD
#' summaries, and standard diagnostic plots. PNG and PDF are written by default.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param chromosome Optional chromosome or scaffold names to scan. By default,
#'   all chromosomes represented by at least one complete window are scanned.
#' @param window.snps Number of SNPs per window.
#' @param step.snps Number of SNPs between consecutive window starts. Defaults
#'   to `window.snps` (non-overlapping windows).
#' @param window.bp Optional fixed physical window size in base pairs. When
#'   supplied, physical windows are used instead of fixed-SNP windows.
#' @param step.bp Distance in base pairs between physical window starts.
#'   Defaults to `window.bp`.
#' @param sensitivity.window.snps Additional fixed-SNP window sizes used for a
#'   sensitivity analysis. These runs summarise PC1 variance and LD along the
#'   genome without replacing the primary candidate scan. Use `NULL` to skip
#'   this additional work.
#' @param n.pcs Number of local covariance axes retained per window.
#' @param mds.axes Number of MDS axes used to score unusual windows.
#' @param outlier.quantile Quantile of the robust window score used as the
#'   candidate threshold.
#' @param min.window.snps Minimum number of usable polymorphic SNPs required in
#'   a window after missing-data and variance checks.
#' @param min.call.rate Minimum genotype call rate required for a SNP within a
#'   window.
#' @param min.candidate.windows Minimum number of consecutive candidate windows
#'   required to form a candidate region.
#' @param cluster.k Number of regional PCA clusters. The biological expectation
#'   for a common polymorphic inversion is often three, representing the two
#'   homokaryotypes and their heterokaryotype, but this is diagnostic rather
#'   than proof.
#' @param ld.max.snps Maximum number of evenly spaced SNPs used for each regional
#'   LD matrix. This bounds memory use without changing the GDS input.
#' @param return.ld Logical indicating whether sampled regional LD matrices are
#'   retained in the result.
#' @param save.plots Logical indicating whether standard PDF and PNG diagnostic
#'   figures are written to the results folder.
#' @param plot.formats One or more of `"png"` and `"pdf"`.
#' @param random.seed Integer seed used for reproducible k-means clustering.
#' @param verbose Logical indicating whether progress messages are printed.
#' @param ... Further arguments for the standard radr workflow. Use
#'   `path.folder` to choose the parent results directory.
#'
#' @return An object of class `detect_inversions` containing:
#'   \itemize{
#'   \item `windows`: coordinates, quality statistics, MDS coordinates, robust
#'     scores, and candidate flags for all analysed windows;
#'   \item `candidates`: one row per contiguous candidate region;
#'   \item `diagnostics`: regional PCA scores, cluster assignments,
#'     heterozygosity summaries, and optional LD matrices;
#'   \item `sensitivity`: optional summaries for additional fixed-SNP window
#'     sizes;
#'   \item `path.folder` and `output.files`: locations of written results;
#'   \item `settings`: the effective analysis settings.
#'   }
#'
#' @export
detect_inversions <- function(
    data,
    chromosome = NULL,
    window.snps = 100L,
    step.snps = window.snps,
    window.bp = NULL,
    step.bp = window.bp,
    sensitivity.window.snps = c(100L, 250L, 500L, 1000L),
    n.pcs = 2L,
    mds.axes = 2L,
    outlier.quantile = 0.99,
    min.window.snps = max(10L, n.pcs + 2L),
    min.call.rate = 0.8,
    min.candidate.windows = 1L,
    cluster.k = 3L,
    ld.max.snps = 500L,
    return.ld = FALSE,
    save.plots = TRUE,
    plot.formats = c("png", "pdf"),
    random.seed = 42L,
    verbose = TRUE,
    ...
) {
  .start <- tgbase::startup(
    package = "radr",
    f.name = "detect_inversions",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "internal"),
    deprecated = NULL,
    verbose = FALSE
  )

  path.folder <- tgbase::generate_folder(
    folder = "detect_inversions",
    path.folder = path.folder,
    internal = internal,
    file.date = file.date,
    prefix.int = TRUE,
    verbose = verbose
  )
  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_detect_inversions_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  if (missing(data)) rlang::abort("Argument `data` is required.")

  integer.args <- list(
    window.snps = window.snps,
    step.snps = step.snps,
    n.pcs = n.pcs,
    mds.axes = mds.axes,
    min.window.snps = min.window.snps,
    min.candidate.windows = min.candidate.windows,
    cluster.k = cluster.k,
    ld.max.snps = ld.max.snps,
    random.seed = random.seed
  )
  bad.integer <- vapply(
    integer.args,
    function(x) length(x) != 1L || !is.numeric(x) || is.na(x) ||
      !is.finite(x) || x < 1 || x != as.integer(x),
    logical(1)
  )
  if (any(bad.integer)) {
    rlang::abort(paste0(
      "These arguments must be positive whole numbers: ",
      paste(names(integer.args)[bad.integer], collapse = ", "), "."
    ))
  }
  integer.args <- lapply(integer.args, as.integer)
  window.snps <- integer.args$window.snps
  step.snps <- integer.args$step.snps
  n.pcs <- integer.args$n.pcs
  mds.axes <- integer.args$mds.axes
  min.window.snps <- integer.args$min.window.snps
  min.candidate.windows <- integer.args$min.candidate.windows
  cluster.k <- integer.args$cluster.k
  ld.max.snps <- integer.args$ld.max.snps
  random.seed <- integer.args$random.seed

  if (!is.null(window.bp)) {
    if (length(window.bp) != 1L || !is.numeric(window.bp) || is.na(window.bp) ||
        !is.finite(window.bp) || window.bp < 1) {
      rlang::abort("`window.bp` must be NULL or one positive number.")
    }
    window.bp <- as.numeric(window.bp)
    if (is.null(step.bp)) step.bp <- window.bp
    if (length(step.bp) != 1L || !is.numeric(step.bp) || is.na(step.bp) ||
        !is.finite(step.bp) || step.bp < 1) {
      rlang::abort("`step.bp` must be one positive number when `window.bp` is used.")
    }
    step.bp <- as.numeric(step.bp)
  }
  if (!is.null(sensitivity.window.snps)) {
    if (!is.numeric(sensitivity.window.snps) || anyNA(sensitivity.window.snps) ||
        any(!is.finite(sensitivity.window.snps)) ||
        any(sensitivity.window.snps < 1) ||
        any(sensitivity.window.snps != as.integer(sensitivity.window.snps))) {
      rlang::abort("`sensitivity.window.snps` must contain positive whole numbers.")
    }
    sensitivity.window.snps <- sort(unique(as.integer(sensitivity.window.snps)))
  }

  if (is.null(window.bp) && min.window.snps > window.snps) {
    rlang::abort("`min.window.snps` cannot exceed `window.snps`.")
  }
  .inversion_check_probability(outlier.quantile, "outlier.quantile", open = TRUE)
  .inversion_check_probability(min.call.rate, "min.call.rate", open = FALSE)
  if (!is.logical(return.ld) || length(return.ld) != 1L || is.na(return.ld)) {
    rlang::abort("`return.ld` must be TRUE or FALSE.")
  }
  if (!is.logical(save.plots) || length(save.plots) != 1L || is.na(save.plots)) {
    rlang::abort("`save.plots` must be TRUE or FALSE.")
  }
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }
  set.seed(random.seed)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    rlang::abort("`verbose` must be TRUE or FALSE.")
  }

  opened.here <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data)) {
    gds <- SeqArray::seqOpen(data)
    opened.here <- TRUE
  } else {
    rlang::abort("`data` must be a GDS filepath or an open SeqVarGDSClass object.")
  }
  on.exit({
    if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)

  SeqArray::seqFilterPush(gds)
  on.exit(try(SeqArray::seqFilterPop(gds), silent = TRUE), add = TRUE)

  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  chrom <- as.character(SeqArray::seqGetData(gds, "chromosome"))
  position <- suppressWarnings(as.numeric(SeqArray::seqGetData(gds, "position")))

  if (length(sample.id) < max(cluster.k, n.pcs + 1L)) {
    rlang::abort("Too few samples for the requested PCA and clustering settings.")
  }
  if (length(variant.id) != length(chrom) || length(chrom) != length(position)) {
    rlang::abort("GDS variant IDs, chromosomes, and positions have inconsistent lengths.")
  }
  if (anyNA(chrom) || any(!nzchar(chrom)) || anyNA(position)) {
    rlang::abort("All scanned GDS variants require chromosome and position values.")
  }

  if (!is.null(chromosome)) {
    chromosome <- unique(as.character(chromosome))
    unknown <- setdiff(chromosome, unique(chrom))
    if (length(unknown) > 0L) {
      rlang::abort(paste0(
        "Chromosome(s) not found in the active GDS variants: ",
        paste(unknown, collapse = ", "), "."
      ))
    }
    keep <- chrom %in% chromosome
    variant.id <- variant.id[keep]
    chrom <- chrom[keep]
    position <- position[keep]
  }

  marker.order <- order(.inversion_chromosome_order(chrom), position, variant.id)
  marker.table <- data.frame(
    variant_id = variant.id[marker.order],
    chromosome = chrom[marker.order],
    position = position[marker.order],
    stringsAsFactors = FALSE
  )
  windows <- .inversion_make_windows(
    marker.table = marker.table,
    window.snps = window.snps,
    step.snps = step.snps,
    window.bp = window.bp,
    step.bp = step.bp
  )
  if (length(windows) < 3L) {
    rlang::abort(
      "At least three complete chromosome-specific windows are required."
    )
  }

  if (verbose) {
    message(
      "Analysing ", length(windows), " windows across ",
      length(unique(vapply(windows, `[[`, character(1), "chromosome"))),
      " chromosome(s)..."
    )
  }

  window.results <- lapply(seq_along(windows), function(i) {
    w <- windows[[i]]
    dosage <- .inversion_get_dosage(gds, w$variant_id, sample.id)
    local <- .inversion_local_covariance(
      dosage = dosage,
      n.pcs = n.pcs,
      min.call.rate = min.call.rate,
      min.window.snps = min.window.snps
    )
    local$window_id <- i
    local$chromosome <- w$chromosome
    local$start <- w$start
    local$end <- w$end
    local$n_input_snps <- length(w$variant_id)
    local$variant_id <- w$variant_id
    local
  })

  valid <- vapply(window.results, function(x) isTRUE(x$valid), logical(1))
  if (sum(valid) < 3L) {
    rlang::abort(
      "Fewer than three windows passed call-rate and polymorphism checks."
    )
  }
  if (verbose && any(!valid)) {
    message(sum(!valid), " window(s) failed quality checks and were not scored.")
  }

  valid.results <- window.results[valid]
  distance.matrix <- .inversion_covariance_distances(valid.results)
  mds.k <- min(mds.axes, nrow(distance.matrix) - 1L)
  mds <- stats::cmdscale(stats::as.dist(distance.matrix), k = mds.k)
  if (is.null(dim(mds))) mds <- matrix(mds, ncol = 1L)
  colnames(mds) <- paste0("MDS", seq_len(ncol(mds)))
  robust.score <- .inversion_robust_score(mds)
  threshold <- as.numeric(stats::quantile(
    robust.score,
    probs = outlier.quantile,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  ))

  window.table <- do.call(rbind, lapply(window.results, function(x) {
    data.frame(
      window_id = x$window_id,
      chromosome = x$chromosome,
      start = x$start,
      end = x$end,
      n_input_snps = x$n_input_snps,
      n_used_snps = x$n_used_snps,
      mean_call_rate = x$mean_call_rate,
      mean_ld_r2 = x$mean_ld_r2,
      valid = x$valid,
      stringsAsFactors = FALSE
    )
  }))
  window.table$robust_score <- NA_real_
  window.table$candidate_window <- FALSE
  valid.index <- which(valid)
  window.table$robust_score[valid.index] <- robust.score
  window.table$candidate_window[valid.index] <- robust.score >= threshold
  for (j in seq_len(ncol(mds))) {
    window.table[[colnames(mds)[j]]] <- NA_real_
    window.table[[colnames(mds)[j]]][valid.index] <- mds[, j]
  }

  candidate.regions <- .inversion_candidate_regions(
    window.table = window.table,
    min.windows = min.candidate.windows
  )

  diagnostics <- vector("list", nrow(candidate.regions))
  if (nrow(candidate.regions) > 0L) {
    if (verbose) {
      message(nrow(candidate.regions), " candidate region(s) selected for diagnostics.")
    }
    for (i in seq_len(nrow(candidate.regions))) {
      region.windows <- window.table$window_id[
        window.table$chromosome == candidate.regions$chromosome[i] &
          window.table$start <= candidate.regions$end[i] &
          window.table$end >= candidate.regions$start[i]
      ]
      region.variant.id <- unique(unlist(
        lapply(window.results[region.windows], `[[`, "variant_id"),
        use.names = FALSE
      ))
      dosage <- .inversion_get_dosage(gds, region.variant.id, sample.id)
      diagnostics[[i]] <- .inversion_region_diagnostics(
        dosage = dosage,
        sample.id = sample.id,
        cluster.k = cluster.k,
        min.call.rate = min.call.rate,
        ld.max.snps = ld.max.snps,
        return.ld = return.ld
      )
      diagnostics[[i]]$candidate_id <- candidate.regions$candidate_id[i]
      diagnostics[[i]]$chromosome <- candidate.regions$chromosome[i]
      diagnostics[[i]]$start <- candidate.regions$start[i]
      diagnostics[[i]]$end <- candidate.regions$end[i]
      diagnostics[[i]]$ld_position <- marker.table$position[
        match(diagnostics[[i]]$ld_variant_id, marker.table$variant_id)
      ]
      candidate.regions$n_samples[i] <- nrow(diagnostics[[i]]$scores)
      candidate.regions$n_snps[i] <- diagnostics[[i]]$n_snps
      candidate.regions$regional_mean_ld_r2[i] <- diagnostics[[i]]$mean_ld_r2
      candidate.regions$three_cluster_pattern[i] <-
        cluster.k == 3L && diagnostics[[i]]$all_clusters_present
      candidate.regions$heterozygote_like_middle_cluster[i] <-
        diagnostics[[i]]$heterozygote_like_middle_cluster
    }
  }

  sensitivity <- .inversion_window_sensitivity(
    gds = gds,
    marker.table = marker.table,
    sample.id = sample.id,
    window.sizes = sensitivity.window.snps,
    min.call.rate = min.call.rate,
    min.window.snps = min.window.snps
  )

  output.files <- .inversion_write_outputs(
    path.folder = path.folder,
    window.table = window.table,
    candidate.regions = candidate.regions,
    diagnostics = diagnostics,
    sensitivity = sensitivity,
    threshold = threshold,
    save.plots = save.plots,
    plot.formats = plot.formats,
    verbose = verbose
  )

  if (!return.ld) {
    diagnostics <- lapply(diagnostics, function(x) {
      x$ld_matrices <- NULL
      x
    })
  }

  out <- list(
    windows = tibble::as_tibble(window.table),
    candidates = tibble::as_tibble(candidate.regions),
    diagnostics = diagnostics,
    sensitivity = sensitivity,
    path.folder = path.folder,
    output.files = output.files,
    settings = list(
      chromosome = chromosome,
      window.snps = window.snps,
      step.snps = step.snps,
      window.bp = window.bp,
      step.bp = step.bp,
      sensitivity.window.snps = sensitivity.window.snps,
      n.pcs = n.pcs,
      mds.axes = mds.k,
      outlier.quantile = outlier.quantile,
      score.threshold = threshold,
      min.window.snps = min.window.snps,
      min.call.rate = min.call.rate,
      min.candidate.windows = min.candidate.windows,
      cluster.k = cluster.k,
      ld.max.snps = ld.max.snps,
      random.seed = random.seed
    )
  )
  class(out) <- c("detect_inversions", class(out))
  out
}

.inversion_check_probability <- function(x, name, open) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (open) valid <- valid && x > 0 && x < 1
  if (!open) valid <- valid && x >= 0 && x <= 1
  if (!valid) {
    interval <- if (open) "strictly between 0 and 1" else "between 0 and 1"
    rlang::abort(paste0("`", name, "` must be ", interval, "."))
  }
  invisible(x)
}

.inversion_chromosome_order <- function(x) {
  stripped <- sub("^(chr|chromosome)", "", x, ignore.case = TRUE)
  numeric.part <- suppressWarnings(as.numeric(stripped))
  ifelse(is.na(numeric.part), Inf, numeric.part) * 1e6 +
    match(x, unique(x))
}

.inversion_make_windows <- function(
    marker.table, window.snps, step.snps, window.bp = NULL, step.bp = window.bp
) {
  split.markers <- split(marker.table, marker.table$chromosome, drop = TRUE)
  windows <- unlist(lapply(split.markers, function(x) {
    x <- x[order(x$position, x$variant_id), , drop = FALSE]
    if (!is.null(window.bp)) {
      first <- floor(min(x$position) / step.bp) * step.bp
      starts <- seq(first, max(x$position), by = step.bp)
      out <- lapply(starts, function(start) {
        idx <- which(x$position >= start & x$position < start + window.bp)
        if (!length(idx)) return(NULL)
        list(
          chromosome = x$chromosome[idx[1L]],
          start = start,
          end = start + window.bp - 1,
          variant_id = x$variant_id[idx]
        )
      })
      return(Filter(Negate(is.null), out))
    }
    if (nrow(x) < window.snps) return(list())
    starts <- seq.int(1L, nrow(x) - window.snps + 1L, by = step.snps)
    lapply(starts, function(first) {
      idx <- seq.int(first, first + window.snps - 1L)
      list(
        chromosome = x$chromosome[first],
        start = min(x$position[idx]),
        end = max(x$position[idx]),
        variant_id = x$variant_id[idx]
      )
    })
  }), recursive = FALSE)
  windows
}

.inversion_get_dosage <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  dosage <- SeqArray::seqGetData(gds, "$dosage_alt")
  dosage <- as.matrix(dosage)

  if (nrow(dosage) == length(sample.id)) {
    out <- dosage
  } else if (ncol(dosage) == length(sample.id)) {
    out <- t(dosage)
  } else {
    rlang::abort("The GDS dosage matrix dimensions do not match its sample IDs.")
  }
  storage.mode(out) <- "double"
  rownames(out) <- sample.id
  colnames(out) <- as.character(SeqArray::seqGetData(gds, "variant.id"))
  out
}

.inversion_prepare_dosage <- function(dosage, min.call.rate) {
  call.rate <- colMeans(!is.na(dosage))
  keep <- is.finite(call.rate) & call.rate >= min.call.rate
  dosage <- dosage[, keep, drop = FALSE]
  call.rate <- call.rate[keep]
  if (ncol(dosage) == 0L) {
    return(list(dosage = dosage, observed = dosage, call.rate = call.rate))
  }

  observed <- dosage
  means <- colMeans(dosage, na.rm = TRUE)
  for (j in seq_len(ncol(dosage))) {
    missing <- is.na(dosage[, j])
    if (any(missing)) dosage[missing, j] <- means[j]
  }
  variance <- apply(dosage, 2L, stats::var)
  keep <- is.finite(variance) & variance > sqrt(.Machine$double.eps)
  list(
    dosage = dosage[, keep, drop = FALSE],
    observed = observed[, keep, drop = FALSE],
    call.rate = call.rate[keep]
  )
}

.inversion_mean_ld <- function(dosage) {
  if (ncol(dosage) < 2L) return(NA_real_)
  r <- suppressWarnings(stats::cor(dosage, use = "pairwise.complete.obs"))
  values <- r[upper.tri(r)]^2
  if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
}

.inversion_local_covariance <- function(
    dosage, n.pcs, min.call.rate, min.window.snps
) {
  prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
  dosage <- prepared$dosage
  base <- list(
    valid = FALSE,
    n_used_snps = ncol(dosage),
    mean_call_rate = if (length(prepared$call.rate)) {
      mean(prepared$call.rate)
    } else {
      NA_real_
    },
    mean_ld_r2 = .inversion_mean_ld(prepared$observed),
    values = NULL,
    vectors = NULL
  )
  if (ncol(dosage) < min.window.snps) return(base)

  centred <- sweep(dosage, 2L, colMeans(dosage), "-")
  covariance <- tcrossprod(centred) / max(1, ncol(centred) - 1L)
  eig <- eigen(covariance, symmetric = TRUE)
  keep <- which(eig$values > sqrt(.Machine$double.eps))
  keep <- utils::head(keep, n.pcs)
  if (length(keep) < n.pcs) return(base)
  values <- eig$values[keep]
  values <- values / sum(abs(values))

  base$valid <- TRUE
  base$values <- values
  base$vectors <- eig$vectors[, keep, drop = FALSE]
  base
}

.inversion_covariance_distances <- function(results) {
  n <- length(results)
  out <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      cross <- crossprod(results[[i]]$vectors, results[[j]]$vectors)
      cross.term <- sum(
        outer(results[[i]]$values, results[[j]]$values) * cross^2
      )
      squared <- sum(results[[i]]$values^2) +
        sum(results[[j]]$values^2) - 2 * cross.term
      out[i, j] <- out[j, i] <- sqrt(max(0, squared))
    }
  }
  out
}

.inversion_robust_score <- function(x) {
  centre <- apply(x, 2L, stats::median, na.rm = TRUE)
  scale <- apply(x, 2L, stats::mad, constant = 1.4826, na.rm = TRUE)
  fallback <- apply(x, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <-
    fallback[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)]
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
  z <- sweep(sweep(x, 2L, centre, "-"), 2L, scale, "/")
  sqrt(rowSums(z^2))
}

.inversion_candidate_regions <- function(window.table, min.windows) {
  candidate <- window.table[window.table$candidate_window, , drop = FALSE]
  empty <- data.frame(
    candidate_id = character(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    n_windows = integer(),
    max_robust_score = numeric(),
    n_samples = integer(),
    n_snps = integer(),
    regional_mean_ld_r2 = numeric(),
    three_cluster_pattern = logical(),
    heterozygote_like_middle_cluster = logical(),
    stringsAsFactors = FALSE
  )
  if (nrow(candidate) == 0L) return(empty)

  candidate <- candidate[order(
    .inversion_chromosome_order(candidate$chromosome),
    candidate$window_id
  ), , drop = FALSE]
  previous <- c(FALSE,
    candidate$chromosome[-1L] == candidate$chromosome[-nrow(candidate)] &
      candidate$window_id[-1L] == candidate$window_id[-nrow(candidate)] + 1L
  )
  group <- cumsum(!previous)
  groups <- split(candidate, group)
  groups <- groups[vapply(groups, nrow, integer(1)) >= min.windows]
  if (length(groups) == 0L) return(empty)

  rows <- lapply(seq_along(groups), function(i) {
    x <- groups[[i]]
    data.frame(
      candidate_id = paste0("INV-CAND-", i),
      chromosome = x$chromosome[1],
      start = min(x$start),
      end = max(x$end),
      n_windows = nrow(x),
      max_robust_score = max(x$robust_score),
      n_samples = NA_integer_,
      n_snps = NA_integer_,
      regional_mean_ld_r2 = NA_real_,
      three_cluster_pattern = FALSE,
      heterozygote_like_middle_cluster = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.inversion_region_diagnostics <- function(
    dosage, sample.id, cluster.k, min.call.rate, ld.max.snps, return.ld
) {
  prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
  dosage <- prepared$dosage
  if (ncol(dosage) < 2L) {
    rlang::abort("A candidate region contains fewer than two usable SNPs.")
  }

  centred <- sweep(dosage, 2L, colMeans(dosage), "-")
  pca <- stats::prcomp(centred, center = FALSE, scale. = FALSE, rank. = 2L)
  scores <- as.data.frame(pca$x[, seq_len(min(2L, ncol(pca$x))), drop = FALSE])
  names(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$individual <- sample.id

  pc1 <- pca$x[, 1L]
  distinct.pc1 <- length(unique(signif(pc1, digits = 12L)))
  if (distinct.pc1 >= cluster.k) {
    cluster <- stats::kmeans(
      pc1, centers = cluster.k, nstart = 50L
    )$cluster
  } else {
    cluster <- rep(1L, length(pc1))
  }
  cluster.means <- tapply(pc1, cluster, mean)
  ordered.cluster <- order(cluster.means)
  cluster.labels <- match(cluster, ordered.cluster)
  scores$cluster <- factor(cluster.labels, levels = seq_len(cluster.k))

  observed <- prepared$observed
  heterozygosity <- rowMeans(observed == 1, na.rm = TRUE)
  scores$heterozygosity <- heterozygosity
  heterozygosity.by.cluster <- tapply(
    heterozygosity, scores$cluster, mean, na.rm = TRUE
  )
  middle <- ceiling(cluster.k / 2)
  heterozygote.like <- cluster.k == 3L &&
    length(heterozygosity.by.cluster) == 3L &&
    is.finite(heterozygosity.by.cluster[middle]) &&
    heterozygosity.by.cluster[middle] == max(heterozygosity.by.cluster)

  ld.index <- unique(round(seq(1, ncol(dosage), length.out = min(
    ncol(dosage), ld.max.snps
  ))))
  observed.ld <- observed[, ld.index, drop = FALSE]
  ld.all <- .inversion_ld_matrix(observed.ld)
  mean.ld <- .inversion_ld_summary(ld.all)

  cluster.levels <- levels(scores$cluster)
  ld.by.cluster <- lapply(cluster.levels, function(level) {
    .inversion_ld_matrix(observed.ld[scores$cluster == level, , drop = FALSE])
  })
  names(ld.by.cluster) <- paste0("cluster_", cluster.levels)

  outer <- c(1L, cluster.k)
  outer.counts <- table(scores$cluster)[outer]
  common.outer <- outer[which.max(outer.counts)]
  common.name <- paste0("cluster_", common.outer)
  ld.common <- ld.by.cluster[[common.name]]
  ld.summary <- tibble::tibble(
    group = c("all", names(ld.by.cluster), "common_homokaryotype"),
    n_samples = c(
      nrow(observed.ld),
      as.integer(table(scores$cluster)),
      as.integer(table(scores$cluster)[common.outer])
    ),
    mean_r2 = c(
      mean.ld,
      vapply(ld.by.cluster, .inversion_ld_summary, numeric(1)),
      .inversion_ld_summary(ld.common)
    )
  )

  list(
    n_snps = ncol(dosage),
    scores = tibble::as_tibble(scores),
    cluster_summary = tibble::tibble(
      cluster = factor(seq_len(cluster.k)),
      n = as.integer(table(scores$cluster)),
      mean_heterozygosity = as.numeric(heterozygosity.by.cluster)
    ),
    all_clusters_present = all(table(scores$cluster) > 0L),
    heterozygote_like_middle_cluster = heterozygote.like,
    mean_ld_r2 = mean.ld,
    ld_variant_index = ld.index,
    ld_variant_id = colnames(observed)[ld.index],
    common_homokaryotype_cluster = common.outer,
    ld_summary = ld.summary,
    ld_matrices = c(
      list(all = ld.all),
      ld.by.cluster,
      list(common_homokaryotype = ld.common)
    )
  )
}

.inversion_ld_matrix <- function(dosage) {
  if (nrow(dosage) < 3L || ncol(dosage) < 2L) {
    return(matrix(NA_real_, nrow = ncol(dosage), ncol = ncol(dosage)))
  }
  suppressWarnings(
    stats::cor(dosage, use = "pairwise.complete.obs")^2
  )
}

.inversion_ld_summary <- function(ld) {
  if (ncol(ld) < 2L) return(NA_real_)
  values <- ld[upper.tri(ld)]
  if (!length(values) || all(!is.finite(values))) return(NA_real_)
  mean(values[is.finite(values)])
}

.inversion_window_sensitivity <- function(
    gds, marker.table, sample.id, window.sizes, min.call.rate,
    min.window.snps
) {
  empty <- tibble::tibble(
    window_snps = integer(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    n_used_snps = integer(),
    pc1_variance = numeric(),
    mean_ld_r2 = numeric()
  )
  if (is.null(window.sizes) || !length(window.sizes)) return(empty)

  rows <- lapply(window.sizes, function(size) {
    windows <- .inversion_make_windows(
      marker.table = marker.table,
      window.snps = size,
      step.snps = size
    )
    if (!length(windows)) return(NULL)
    do.call(rbind, lapply(windows, function(w) {
      dosage <- .inversion_get_dosage(gds, w$variant_id, sample.id)
      prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
      if (ncol(prepared$dosage) < min(min.window.snps, size)) return(NULL)
      pca <- stats::prcomp(
        prepared$dosage,
        center = TRUE,
        scale. = FALSE,
        rank. = 1L
      )
      variance <- pca$sdev^2
      data.frame(
        window_snps = size,
        chromosome = w$chromosome,
        start = w$start,
        end = w$end,
        n_used_snps = ncol(prepared$dosage),
        pc1_variance = variance[1L] / sum(variance),
        mean_ld_r2 = .inversion_mean_ld(prepared$observed),
        stringsAsFactors = FALSE
      )
    }))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)
  tibble::as_tibble(do.call(rbind, rows))
}

.inversion_write_outputs <- function(
    path.folder, window.table, candidate.regions, diagnostics, sensitivity,
    threshold, save.plots, plot.formats, verbose
) {
  files <- character()
  write_table <- function(x, filename) {
    path <- file.path(path.folder, filename)
    readr::write_tsv(x, path, na = "NA")
    files <<- c(files, path)
  }
  write_table(window.table, "inversion_windows.tsv")
  write_table(candidate.regions, "candidate_inversion_regions.tsv")
  if (nrow(sensitivity) > 0L) {
    write_table(sensitivity, "window_size_sensitivity.tsv")
  }
  for (i in seq_along(diagnostics)) {
    prefix <- diagnostics[[i]]$candidate_id
    write_table(
      diagnostics[[i]]$scores,
      paste0(prefix, "_individual_pca.tsv")
    )
    write_table(
      diagnostics[[i]]$cluster_summary,
      paste0(prefix, "_cluster_summary.tsv")
    )
    write_table(
      diagnostics[[i]]$ld_summary,
      paste0(prefix, "_ld_summary.tsv")
    )
  }

  plots <- list()
  if (save.plots) {
    midpoint <- (window.table$start + window.table$end) / 2
    plot.data <- transform(window.table, midpoint_mb = midpoint / 1e6)
    plots$window_scores <- ggplot2::ggplot(
      plot.data,
      ggplot2::aes(
        x = midpoint_mb, y = robust_score, colour = candidate_window
      )
    ) +
      ggplot2::geom_hline(yintercept = threshold, linetype = 2L) +
      ggplot2::geom_line(ggplot2::aes(group = chromosome), colour = "grey70") +
      ggplot2::geom_point(size = 1.7) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)",
        y = "Robust local-structure score",
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    if (all(c("MDS1", "MDS2") %in% names(window.table))) {
      plots$window_mds <- ggplot2::ggplot(
        window.table,
        ggplot2::aes(
          x = MDS1, y = MDS2, colour = chromosome,
          shape = candidate_window
        )
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::labs(
          colour = "Chromosome / LG",
          shape = "Candidate window"
        ) +
        ggplot2::theme_bw()
    }

    plots$window_ld <- ggplot2::ggplot(
      plot.data,
      ggplot2::aes(x = midpoint_mb, y = mean_ld_r2)
    ) +
      ggplot2::geom_line(colour = "grey50") +
      ggplot2::geom_point(
        ggplot2::aes(colour = candidate_window), size = 1.5
      ) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)", y = expression("Mean pairwise " * r^2),
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    plots$window_call_rate <- ggplot2::ggplot(
      plot.data,
      ggplot2::aes(x = midpoint_mb, y = mean_call_rate)
    ) +
      ggplot2::geom_line(colour = "grey50") +
      ggplot2::geom_point(
        ggplot2::aes(colour = candidate_window), size = 1.5
      ) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)", y = "Mean SNP call rate",
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    if (nrow(sensitivity) > 0L) {
      sensitivity.plot <- transform(
        sensitivity,
        midpoint_mb = (start + end) / 2e6,
        window_snps = factor(window_snps)
      )
      plots$window_sensitivity <- ggplot2::ggplot(
        sensitivity.plot,
        ggplot2::aes(
          x = midpoint_mb, y = pc1_variance, colour = window_snps,
          group = window_snps
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
        ggplot2::labs(
          x = "Genomic position (Mb)", y = "PC1 variance proportion",
          colour = "SNPs per window"
        ) +
        ggplot2::theme_bw()
    }

    for (i in seq_along(diagnostics)) {
      d <- diagnostics[[i]]
      plots[[paste0(d$candidate_id, "_pca")]] <- ggplot2::ggplot(
        d$scores,
        ggplot2::aes(x = PC1, y = PC2, colour = cluster)
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::labs(
          title = paste0(
            d$candidate_id, " | ", d$chromosome, ":",
            format(d$start, big.mark = ","), "-",
            format(d$end, big.mark = ",")
          ),
          colour = "Putative arrangement genotype"
        ) +
        ggplot2::theme_bw()

      plots[[paste0(d$candidate_id, "_heterozygosity")]] <-
        ggplot2::ggplot(
          d$scores,
          ggplot2::aes(x = cluster, y = heterozygosity, colour = cluster)
        ) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
        ggplot2::labs(
          x = "Putative arrangement genotype",
          y = "Observed heterozygosity"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

      if (length(d$ld_position) > 1L) {
        ld.all <- d$ld_matrices$all
        ld.common <- d$ld_matrices$common_homokaryotype
        grid <- expand.grid(
          row = seq_along(d$ld_position),
          column = seq_along(d$ld_position)
        )
        use.all <- grid$row <= grid$column
        grid$r2 <- ifelse(
          use.all,
          ld.all[cbind(grid$row, grid$column)],
          ld.common[cbind(grid$row, grid$column)]
        )
        grid$position1 <- d$ld_position[grid$row] / 1e6
        grid$position2 <- d$ld_position[grid$column] / 1e6
        grid$group <- ifelse(use.all, "All individuals", "Common homokaryotype")
        plots[[paste0(d$candidate_id, "_ld")]] <- ggplot2::ggplot(
          grid,
          ggplot2::aes(x = position1, y = position2, fill = r2)
        ) +
          ggplot2::geom_raster() +
          ggplot2::coord_equal() +
          ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
          ggplot2::labs(
            x = "Genomic position (Mb)", y = "Genomic position (Mb)",
            fill = expression(r^2),
            subtitle = paste0(
              "Upper triangle: all individuals; lower triangle: cluster ",
              d$common_homokaryotype_cluster
            )
          ) +
          ggplot2::theme_bw()
      }
    }

    for (name in names(plots)) {
      for (format in plot.formats) {
        path <- file.path(path.folder, paste0(name, ".", format))
        ggplot2::ggsave(
          filename = path,
          plot = plots[[name]],
          width = if (grepl("_pca$|_heterozygosity$", name)) 6 else 10,
          height = if (grepl("_pca$|_heterozygosity$", name)) 5 else 7,
          dpi = 300
        )
        files <- c(files, path)
      }
    }
  }
  if (verbose) message("Inversion results written to: ", path.folder)
  list(files = normalizePath(files, mustWork = FALSE), plots = plots)
}

#' @export
print.detect_inversions <- function(x, ...) {
  cat("Candidate inversion-associated region scan\n")
  cat("  Windows analysed:", sum(x$windows$valid), "of", nrow(x$windows), "\n")
  cat("  Candidate regions:", nrow(x$candidates), "\n")
  cat("  Results folder:", x$path.folder, "\n")
  if (nrow(x$candidates) > 0L) print(x$candidates)
  invisible(x)
}
