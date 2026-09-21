# Resolve explicit identifiers against the baseline, never the full unfiltered GDS.
.diagnostic_pca_scenarios <- function(sensitivity, ids, variants) {
  if(is.null(sensitivity)) return(list())
  labels <- names(sensitivity)
  if(!is.list(sensitivity) || !length(sensitivity) || is.null(labels) ||
     anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels))
    rlang::abort("`sensitivity` must be a nonempty, uniquely named list.")
  lapply(sensitivity,function(scenario) {
    if(!is.list(scenario) || !length(scenario) || is.null(names(scenario)) ||
       anyDuplicated(names(scenario)) ||
       any(!names(scenario) %in% c("sample.id","variant.id")))
      rlang::abort("Each sensitivity scenario must specify sample.id and/or variant.id.")
    resolve <- function(selected,baseline,minimum,label) {
      if(is.null(selected)) return(baseline)
      if(!is.atomic(selected) || anyNA(selected) || anyDuplicated(selected) ||
         !all(selected %in% baseline) || length(selected)<minimum)
        rlang::abort(paste0("Invalid sensitivity ",label,": use unique active identifiers (minimum ",minimum,")."))
      baseline[baseline %in% selected]
    }
    list(sample.id=resolve(scenario$sample.id,ids,3,"sample.id"),
         variant.id=resolve(scenario$variant.id,variants,2,"variant.id"))
  })
}

# Normalized orthogonal Procrustes, matched on identifiers before any arithmetic.
.diagnostic_pca_compare <- function(baseline,subset,baseline.axes,subset.axes,label) {
  ids <- baseline$INDIVIDUALS[baseline$INDIVIDUALS %in% subset$INDIVIDUALS]
  k <- min(length(baseline.axes),length(subset.axes))
  summary <- tibble::tibble(SCENARIO=label,N_SHARED=length(ids),N_COMPONENTS=k,
    STATUS="INSUFFICIENT_SHARED_SAMPLES",SIMILARITY=NA_real_,
    RMS_DISPLACEMENT=NA_real_,MAX_DISPLACEMENT=NA_real_)
  samples <- tibble::tibble(SCENARIO=label,INDIVIDUALS=ids,DISPLACEMENT=NA_real_)
  if(k<1 || length(ids)<=k) return(list(summary=summary,samples=samples))
  a <- as.matrix(baseline[match(ids,baseline$INDIVIDUALS),baseline.axes[seq_len(k)],drop=FALSE])
  b <- as.matrix(subset[match(ids,subset$INDIVIDUALS),subset.axes[seq_len(k)],drop=FALSE])
  a <- scale(a,center=TRUE,scale=FALSE); b <- scale(b,center=TRUE,scale=FALSE)
  full.rank <- function(x) {
    values <- svd(x,nu=0,nv=0)$d
    max(values)>0 && sum(values>max(values)*1e-8)==ncol(x)
  }
  if(!all(is.finite(a)) || !all(is.finite(b)) || !full.rank(a) || !full.rank(b)) {
    summary$STATUS <- "RANK_DEFICIENT"
    return(list(summary=summary,samples=samples))
  }
  a <- a/sqrt(sum(a*a)); b <- b/sqrt(sum(b*b))
  decomposition <- svd(crossprod(b,a))
  aligned <- b %*% (decomposition$u %*% t(decomposition$v))
  samples$DISPLACEMENT <- unname(sqrt(rowSums((a-aligned)^2)))
  summary$STATUS <- "OK"
  summary$SIMILARITY <- min(1,max(0,sum(decomposition$d)))
  summary$RMS_DISPLACEMENT <- sqrt(mean(samples$DISPLACEMENT^2))
  summary$MAX_DISPLACEMENT <- max(samples$DISPLACEMENT)
  list(summary=summary,samples=samples)
}
