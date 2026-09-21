#' Inspect genomic principal components before filtering
#'
#' Compute principal component analysis (PCA) of active diploid, biallelic
#' genotype calls and inspect sample scores and genome-wide marker loadings.
#' This read-only diagnostic never blacklists samples or markers.
#'
#' @param data A Genomic Data Structure (GDS) path or open SeqVarGDSClass object.
#' @param strata Optional metadata data frame or tab-delimited path containing
#'   unique INDIVIDUALS. All active samples must be covered. Extra rows are
#'   ignored; missing metadata values are allowed and shown in grey. With NULL,
#'   use genometranslator metadata if available, otherwise sample identifiers.
#' @param metadata.columns Metadata columns to colour scores by. NULL selects
#'   all atomic metadata columns except identifiers and FILTERS. Computed
#'   missingness, heterozygosity and available depth are always plotted.
#' @param n.components Maximum number of components (default four).
#' @param chunk.size Maximum variants per read block.
#' @param max.cells Maximum cells in the sample covariance matrix or a genotype
#'   read block. This is a guard, not a total memory limit; eigendecomposition
#'   and temporary arrays require additional memory.
#' @param reference.genome Logical, or NULL to read the stored reference flag.
#'   Only TRUE enables physical-position loading plots. Otherwise active marker
#'   order is shown, explicitly not genomic distance. Loading plots are paginated
#'   at six sequences per page. Missing physical positions are reported in the
#'   marker table and omitted from physical-position plots.
#'   Without a reference, chromosome labels are ignored only in the loading
#'   figure and all markers share an unmapped panel. Original labels remain in
#'   the marker table. Adjacent plotted markers need not be physically linked.
#' @param write.files Save tables, settings and PNG figures in a standard folder.
#' @param verbose Display progress messages.
#' @param sensitivity Optional named list of explicit subset scenarios. Each
#'   scenario is a list containing sample.id and/or variant.id to retain, using
#'   active GDS identifiers (not row numbers). Omitted dimensions retain the
#'   baseline selection. Unknown, duplicated or inactive identifiers are errors.
#'   No subsets are chosen automatically.
#' @param ... Common output settings: path.folder and internal.
#'
#' @details Allele frequencies are estimated from called active genotypes.
#' Calls are centred at twice allele frequency and scaled by the square root
#' of twice frequency times one minus frequency. Missing calls are temporarily
#' assigned the marker mean (zero after centring). Partial calls are missing.
#' Markers without calls or without observed variation are omitted from the
#' numerical decomposition, but retained in the marker-status table. A sample
#' with no calls among usable markers causes an explicit error, not an invented
#' position. Genotypes and persistent filters are never changed.
#'
#' This implementation uses two blockwise passes and a dense sample covariance
#' matrix. It is not a genotype-likelihood method or a projection method.
#' Mean completion can generate misleading structure when missingness is
#' substantial or structured. Compare coordinates against biological strata,
#' missingness, depth, heterozygosity, plate, library, lane and project.
#' Metadata colourings are descriptive, not batch correction or causal tests.
#' Perfect population/batch confounding cannot be resolved from these plots.
#'
#' No linkage disequilibrium (LD) pruning, iterative region exclusion, outlier
#' deletion or automatic choice of population number is performed. Inspect
#' loadings along chromosomes: a regional peak may represent an inversion or
#' other linked region, not genome-wide population structure. Keep such regions
#' for separate [detect_inversions()] analysis. Compare explicit selections in
#' separate runs; do not silently discard rare populations or relatives.
#' Sampling-site balance does not establish balance of true genetic groups.
#'
#' Sensitivity scenarios refit allele frequencies and PCA on each subset.
#' Shared samples are matched by identifier, centred, and each score matrix is
#' divided by its Frobenius norm. Orthogonal Procrustes alignment permits axis
#' flips and rotations. Similarity is the sum of singular values of the cross
#' product of the normalized score matrices (one means identical geometry).
#' Per-sample displacement is measured after alignment in these normalized
#' coordinates. Scale changes are removed deliberately. Comparisons use the
#' smaller returned component count and require more shared samples than axes
#' and full column rank on shared samples; otherwise metrics are unavailable.
#' These are descriptive, in-sample comparisons, not independent validation,
#' tests of population number, or evidence that an exclusion is appropriate.
#'
#' Sample summaries describe all active markers, before numerical omissions.
#' Depth is mean finite nonnegative FORMAT/DP, including depth at missing calls;
#' unavailable depth is NA, never zero. Depth availability is reported separately.
#' Component variance proportions use total standardized variance, not only the
#' displayed components. Axis signs are arbitrary, and tied eigenspaces may rotate.
#'
#' @return List with samples, markers (including loadings and omission reasons),
#'   components, plots, settings, depth.status, path.folder and
#'   active.selection.restored. Scores use PC1, PC2, etc.; computed sample
#'   diagnostics use the PCA_ prefix. Metadata names that collide are rejected.
#'   When requested, sensitivity contains summary metrics, per-sample aligned
#'   displacements, and exact scenario selections; otherwise it is NULL.
#'   Sensitivity tables and selections are also saved when write.files is TRUE.
#' @references Prive F, Luu K, Blum MGB, McGrath JJ, Vilhjalmsson BJ (2020).
#' Efficient toolkit implementing best practices for principal component
#' analysis of population genetic data. Bioinformatics 36:4449-4457.
#' doi:10.1093/bioinformatics/btaa520
#' @export
#' @examples
#' \dontrun{
#' diagnostics <- diagnostic_pca("study.gds",
#'   metadata.columns = c("STRATA", "PLATE", "LIBRARY"))
#' diagnostics$plots$loadings
#' diagnostics$plots$scores_1_PCA_MISSING
#' # Use an open GDS object imported with genometranslator.
#' individuals <- genometranslator::extract_individuals_metadata(
#'   gds = your.gds, whitelist = TRUE)
#' markers <- genometranslator::extract_markers_metadata(
#'   gds = your.gds, whitelist = TRUE)
#' 
#' # Restrict metadata to the baseline's active identifiers.
#' # A metadata whitelist can be broader than a temporary GDS selection.
#' diagnostics <- diagnostic_pca(your.gds)
#' individuals <- dplyr::filter(
#'   individuals, INDIVIDUALS %in% diagnostics$samples$INDIVIDUALS)
#' markers <- dplyr::filter(
#'   markers, VARIANT_ID %in% diagnostics$markers$VARIANT_ID)
#' 
#' # After reviewing the figures, subset these tables explicitly as appropriate.
#' # Leaving them unchanged is a same-selection control, not a sensitivity test.
#' reviewed_individuals <- individuals
#' reviewed_markers <- markers
#' sensitivity <- diagnostic_pca(your.gds, sensitivity = list(
#'   reviewed_samples = list(sample.id = reviewed_individuals$INDIVIDUALS),
#'   reviewed_markers = list(variant.id = reviewed_markers$VARIANT_ID)))
#' sensitivity$sensitivity$summary
#' }
diagnostic_pca <- function(data, strata = NULL, metadata.columns = NULL,
    n.components = 4L, chunk.size = 1000L, max.cells = 25000000,
    reference.genome = NULL, write.files = TRUE, verbose = TRUE,
    sensitivity = NULL, ...) {
  .paralog_check_count(n.components, "n.components", 1)
  .paralog_check_count(chunk.size, "chunk.size", 1)
  .paralog_check_count(max.cells, "max.cells", 1)
  if (!is.null(reference.genome)) .paralog_check_flag(reference.genome, "reference.genome")
  context <- .analysis_gds_start(data, "diagnostic_pca", write.files, verbose, ...)
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  ids <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variants <- SeqArray::seqGetData(gds, "variant.id")
  n <- length(ids); m <- length(variants)
  scenarios <- .diagnostic_pca_scenarios(sensitivity, ids, variants)
  if (n < 3L || m < 2L) rlang::abort("PCA requires at least three samples and two active markers.")
  if (as.double(n)^2 > max.cells) rlang::abort("Sample covariance exceeds `max.cells`; select a smaller fitting set explicitly.")
  alleles <- SeqArray::seqGetData(gds, "allele")
  if (any(lengths(strsplit(alleles, ",", fixed = TRUE)) != 2L))
    rlang::abort("PCA currently requires biallelic genotype calls.")
  meta <- .diagnostic_pca_metadata(gds, strata, ids)
  if (any(grepl("^PC[0-9]+$|^PCA_", names(meta))))
    rlang::abort("Metadata names beginning PCA_ or matching PC followed by digits are reserved.")
  if (is.null(metadata.columns)) metadata.columns <- setdiff(names(meta), c("INDIVIDUALS", "FILTERS"))
  if (!is.character(metadata.columns) || anyNA(metadata.columns) ||
      !all(metadata.columns %in% names(meta))) rlang::abort("Unknown `metadata.columns`.")
  if (any(!vapply(meta[metadata.columns], is.atomic, logical(1))))
    rlang::abort("Plot metadata must contain atomic columns.")
  if (is.null(reference.genome)) {
    node <- gdsfmt::index.gdsn(gds, "genometranslator/reference.genome", silent = TRUE)
    reference.genome <- !is.null(node) && isTRUE(as.logical(gdsfmt::read.gdsn(node)))
  }
  markers <- tibble::tibble(VARIANT_ID = variants,
    CHROM = as.character(SeqArray::seqGetData(gds, "chromosome")),
    POSITION = suppressWarnings(as.numeric(SeqArray::seqGetData(gds, "position"))),
    ACTIVE_ORDER = seq_len(m))
  markers$POSITION_AVAILABLE <- is.finite(markers$POSITION) & markers$POSITION > 0
  if (reference.genome && !any(markers$POSITION_AVAILABLE)) {
    warning("No valid physical positions; plotting active marker order instead.")
    reference.genome <- FALSE
  }
  depth.available <- !is.null(gdsfmt::index.gdsn(gds, "annotation/format/DP", silent = TRUE))
  blocks <- split(seq_len(m), ceiling(seq_len(m) / min(chunk.size, max(1, floor(max.cells / (2*n))))))
  covariance <- matrix(0, n, n)
  called <- heterozygous <- usable.calls <- depth.sum <- depth.n <- numeric(n)
  frequency <- rep(NA_real_, m); count <- numeric(m); status <- rep("USED", m)
  if (verbose) message("Computing PCA covariance from ", n, " samples and ", m, " markers...")
  for (index in blocks) {
    x <- .diagnostic_pca_calls(gds, variants[index], n)
    called <- called + rowSums(!is.na(x))
    heterozygous <- heterozygous + rowSums(x == 1, na.rm = TRUE)
    count[index] <- colSums(!is.na(x))
    frequency[index] <- colSums(x, na.rm = TRUE) / (2 * count[index])
    valid <- count[index] > 0 & is.finite(frequency[index]) & frequency[index] > 0 & frequency[index] < 1
    # Constant heterozygotes have p=0.5 but no variation after centring.
    varying <- vapply(seq_len(ncol(x)), function(j) length(unique(x[!is.na(x[,j]),j])) > 1L, logical(1))
    valid <- valid & varying
    status[index[!valid]] <- ifelse(count[index[!valid]] == 0, "ALL_MISSING", "NO_VARIATION")
    if (any(valid)) {
      usable.calls <- usable.calls + rowSums(!is.na(x[,valid,drop=FALSE]))
      z <- .diagnostic_pca_standardize(x[,valid,drop=FALSE], frequency[index[valid]])
      covariance <- covariance + tcrossprod(z)
    }
    if (depth.available) {
      dp <- SeqArray::seqGetData(gds, "annotation/format/DP")
      # FORMAT fields may be declared in the header but absent at some sites.
      if (is.list(dp)) {
        lengths <- dp$length
        values <- dp$data
        if (!is.null(lengths) && length(lengths) == length(index) &&
            all(lengths %in% 0:1) && length(values) == n*sum(lengths)) {
          dp <- matrix(NA_real_,n,length(index))
          if(any(lengths==1)) dp[,lengths==1] <- values
        } else dp <- values
      }
      if (!length(dp)) next
      if (length(dp) != n*length(index)) rlang::abort("FORMAT/DP must contain one value per sample and marker.")
      dp <- matrix(dp, n, length(index)); ok <- is.finite(dp) & dp >= 0
      dp[!ok] <- 0
      depth.sum <- depth.sum + rowSums(dp); depth.n <- depth.n + rowSums(ok)
    }
  }
  if (any(usable.calls == 0)) rlang::abort(paste0("No usable genotype calls for: ", paste(ids[usable.calls==0],collapse=", ")))
  eig <- eigen(covariance, symmetric = TRUE)
  positive <- which(eig$values > max(eig$values) * 1e-10)
  k <- min(n.components, length(positive), n-1L)
  if (k < 1L) rlang::abort("No variable PCA dimensions remain.")
  axes <- paste0("PC",seq_len(k)); u <- eig$vectors[,seq_len(k),drop=FALSE]
  singular <- sqrt(eig$values[seq_len(k)])
  scores <- sweep(u,2,singular,"*"); colnames(scores) <- axes
  loadings <- matrix(NA_real_,m,k,dimnames=list(NULL,axes))
  for (index in blocks) {
    selected <- index[status[index]=="USED"]
    if (!length(selected)) next
    z <- .diagnostic_pca_standardize(.diagnostic_pca_calls(gds,variants[selected],n),frequency[selected])
    loadings[selected,] <- sweep(crossprod(z,u),2,singular,"/")
  }
  samples <- dplyr::bind_cols(meta, tibble::as_tibble(scores)) |>
    dplyr::mutate(PCA_MISSING=1-called/m, PCA_HETEROZYGOSITY=heterozygous/called,
                  PCA_DEPTH=ifelse(depth.n>0,depth.sum/depth.n,NA_real_))
  markers <- dplyr::bind_cols(markers,tibble::tibble(STATUS=status,
    CALLED=count,ALT_FREQUENCY=frequency),tibble::as_tibble(loadings))
  components <- tibble::tibble(COMPONENT=axes,
    VARIANCE_PROPORTION=eig$values[seq_len(k)]/sum(pmax(eig$values,0)))
  plots <- .diagnostic_pca_plots(samples,markers,components,
    unique(c(metadata.columns,"PCA_MISSING","PCA_HETEROZYGOSITY",if(any(depth.n>0)) "PCA_DEPTH")),
    reference.genome)
  settings <- list(n.components=n.components,returned.components=k,
    chunk.size=chunk.size,max.cells=max.cells,reference.genome=reference.genome,
    metadata.columns=metadata.columns,sample.ids=ids,variant.ids=variants,
    imputation="Temporary pooled mean; original GDS unchanged",session=sessionInfo())
  depth.status <- if(!depth.available) "UNAVAILABLE" else if(!any(depth.n>0)) "NO_VALID_VALUES" else "AVAILABLE"
  sensitivity.result <- NULL
  if(length(scenarios)) {
    comparisons <- displacements <- vector("list",length(scenarios))
    for(j in seq_along(scenarios)) {
      scenario <- scenarios[[j]]; label <- names(scenarios)[j]
      if(verbose) message("PCA sensitivity: ",label)
      SeqArray::seqSetFilter(gds,sample.id=scenario$sample.id,
        variant.id=scenario$variant.id,verbose=FALSE)
      fit <- diagnostic_pca(gds,strata=strata,metadata.columns=metadata.columns,
        n.components=n.components,chunk.size=chunk.size,max.cells=max.cells,
        reference.genome=reference.genome,write.files=FALSE,verbose=FALSE)
      comparison <- .diagnostic_pca_compare(samples,fit$samples,axes,
        fit$components$COMPONENT,label)
      comparisons[[j]] <- dplyr::mutate(comparison$summary,
        N_SAMPLES=nrow(fit$samples),N_MARKERS=nrow(fit$markers),
        N_USABLE_MARKERS=sum(fit$markers$STATUS=="USED"))
      displacements[[j]] <- comparison$samples
    }
    sensitivity.result <- list(summary=dplyr::bind_rows(comparisons),
      samples=dplyr::bind_rows(displacements),selections=scenarios)
    available <- dplyr::filter(sensitivity.result$samples,is.finite(.data$DISPLACEMENT))
    if(nrow(available)) plots$sensitivity <- ggplot2::ggplot(available,
      ggplot2::aes(x=.data$SCENARIO,y=.data$DISPLACEMENT))+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.15,height=0,seed=1),alpha=.6)+
      ggplot2::theme_bw()+ggplot2::labs(x="Explicit subset scenario",
        y="Aligned sample displacement (normalized coordinates)",
        subtitle="Shared samples only; descriptive sensitivity, not independent validation")
  }
  .analysis_gds_restore(context)
  result <- list(samples=samples,markers=markers,components=components,plots=plots,
    settings=settings,depth.status=depth.status,path.folder=context$path.folder,
    active.selection.restored=context$restored,sensitivity=sensitivity.result)
  if (write.files) {
    for (name in c("samples","markers","components"))
      readr::write_tsv(result[[name]],file.path(context$path.folder,paste0("pca_",name,".tsv")))
    saveRDS(list(settings=settings,depth.status=depth.status),file.path(context$path.folder,"pca_settings.rds"))
    if(!is.null(sensitivity.result)) {
      readr::write_tsv(sensitivity.result$summary,file.path(context$path.folder,"pca_sensitivity_summary.tsv"))
      readr::write_tsv(sensitivity.result$samples,file.path(context$path.folder,"pca_sensitivity_samples.tsv"))
      saveRDS(scenarios,file.path(context$path.folder,"pca_sensitivity_selections.rds"))
    }
    for(name in names(plots)) ggplot2::ggsave(file.path(context$path.folder,paste0(name,".png")),
      plots[[name]],width=10,height=7,dpi=150)
    if(verbose) message("PCA diagnostics written to: ",basename(context$path.folder))
  }
  result
}

.diagnostic_pca_metadata <- function(gds,strata,ids) {
  if(is.null(strata)) {
    node <- gdsfmt::index.gdsn(gds,"genometranslator/individuals.meta",silent=TRUE)
    if(is.null(node)) return(tibble::tibble(INDIVIDUALS=ids))
    strata <- genometranslator::extract_individuals_metadata(gds,whitelist=FALSE)
  } else if(is.character(strata) && length(strata)==1L) {
    strata <- readr::read_tsv(strata,col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
  }
  if(!is.data.frame(strata) || anyDuplicated(names(strata)) || !"INDIVIDUALS" %in% names(strata) ||
     anyNA(strata$INDIVIDUALS) || anyDuplicated(strata$INDIVIDUALS) ||
     !all(ids %in% strata$INDIVIDUALS)) rlang::abort("Metadata must uniquely cover every active INDIVIDUALS identifier.")
  tibble::as_tibble(strata[match(ids,strata$INDIVIDUALS),,drop=FALSE])
}

.diagnostic_pca_calls <- function(gds,variants,n) {
  SeqArray::seqSetFilter(gds,variant.id=variants,verbose=FALSE)
  g <- SeqArray::seqGetData(gds,"genotype")
  if(!identical(dim(g),c(2L,as.integer(n),as.integer(length(variants)))) ||
     any(!is.na(g) & !g %in% 0:1)) rlang::abort("PCA requires diploid biallelic calls.")
  matrix(g[1,,]+g[2,,],n,length(variants))
}

.diagnostic_pca_standardize <- function(x,p) {
  z <- sweep(sweep(x,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/")
  z[is.na(z)] <- 0
  z
}

.diagnostic_pca_plots <- function(samples,markers,components,variables,reference) {
  axes <- components$COMPONENT; out <- list()
  out$scree <- ggplot2::ggplot(components,ggplot2::aes(x=.data$COMPONENT,y=.data$VARIANCE_PROPORTION))+
    ggplot2::geom_col()+ggplot2::labs(x="Principal component",y="Proportion of total variance")+ggplot2::theme_bw()
  long <- tidyr::pivot_longer(markers,tidyselect::all_of(axes),names_to="COMPONENT",values_to="LOADING") |>
    dplyr::filter(.data$STATUS=="USED")
  # Placeholder chromosome labels are not physical mapping information.
  # Collapse only the plotting copy; preserve original labels in the output table.
  if(!reference) long$CHROM <- "Unmapped markers"
  long$CHROM <- factor(long$CHROM,levels=unique(long$CHROM)[stringi::stri_order(unique(long$CHROM),numeric=TRUE)])
  long$LOCATION <- if(reference) ifelse(long$POSITION_AVAILABLE,long$POSITION,NA_real_) else long$ACTIVE_ORDER
  long <- dplyr::filter(long,is.finite(.data$LOCATION))
  sequences <- levels(droplevels(long$CHROM))
  if(length(sequences)) for(page in seq_len(ceiling(length(sequences)/6))) {
    selected <- sequences[seq.int((page-1)*6+1,min(page*6,length(sequences)))]
    panel <- long[long$CHROM %in% selected,,drop=FALSE]
    out[[if(page==1) "loadings" else paste0("loadings_",page)]] <- ggplot2::ggplot(panel,ggplot2::aes(x=.data$LOCATION,y=.data$LOADING))+
    ggplot2::geom_point(size=.5,alpha=.5)+
    ggplot2::facet_grid(COMPONENT~CHROM,scales="free_x")+
    ggplot2::labs(x=if(reference) "Physical position (bases)" else "Active marker order (not genomic distance)",
      y="Marker loading",title=if(reference) "Regional contributions to principal components" else
        "Marker contributions to principal components",
      subtitle=if(!reference) "Marker order does not imply physical proximity" else NULL)+ggplot2::theme_bw()
  }
  for(i in seq(1,length(axes),by=2)) {
    x <- if(length(axes)==1) ".sample_index" else axes[min(i,length(axes)-1)]
    y <- axes[min(i+1,length(axes))]
    for(j in seq_along(variables)) {
      variable <- variables[j]; frame <- samples
      frame$.sample_index <- seq_len(nrow(frame))
      frame$.colour <- frame[[variable]]
      if(is.numeric(frame$.colour)) frame$.colour[!is.finite(frame$.colour)] <- NA_real_
      if(!is.numeric(frame$.colour)) frame$.colour <- as.factor(frame$.colour)
      plot <- ggplot2::ggplot(frame,ggplot2::aes(x=.data[[x]],y=.data[[y]],colour=.data$.colour))+
        ggplot2::geom_point(size=1.8,alpha=.8)+ggplot2::theme_bw()+
        ggplot2::labs(x=if(x==".sample_index") "Active sample order" else x,
          y=y,colour=variable,title=paste("Sample coordinates by",variable),
          subtitle="Descriptive diagnostic: no sample removal or batch correction")
      if(is.numeric(frame$.colour)) plot <- plot+ggplot2::scale_colour_viridis_c(na.value="grey70")
      else plot <- plot+ggplot2::scale_colour_discrete(na.value="grey70",na.translate=TRUE)
      key <- paste0("scores_",i,"_",gsub("[^A-Za-z0-9_]","_",variable))
      if(key %in% names(out)) key <- paste0(key,"_",j)
      out[[key]] <- plot
    }
  }
  out
}
