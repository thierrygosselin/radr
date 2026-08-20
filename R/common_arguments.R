#' Common arguments used by radr
#'
#' This documentation-only helper centralizes parameters shared by genomic
#' screening, filtering, and diagnostic functions.
#'
#' @name radr_common_arguments
#' @rdname radr_common_arguments
#' @keywords internal
#' @export
#'
#' @param interactive.filter Logical indicating whether an interactive filtering
#'   session may display diagnostics and ask for thresholds.
#'   Default: \code{interactive.filter = TRUE}.
#' @param gds A genome GDS file path or object supported by `genometranslator`.
#' @param data A tidy genomic data frame or another genomic object supported by
#'   the calling function.
#' @param parallel.core Number of workers available for parallel operations.
#'   Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param verbose Logical indicating whether progress messages are emitted.
#'   Default: \code{verbose = TRUE}.
#' @param random.seed Optional integer seed for reproducible operations.
#'   Default: \code{random.seed = NULL}.
#' @param ... Additional arguments passed to lower-level screening or filtering
#'   functions.
#'
#' @return `NULL`, invisibly. This function exists to share documentation.
radr_common_arguments <- function(
    interactive.filter = TRUE,
    gds = NULL,
    data = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE,
    random.seed = NULL,
    ...
) {
  invisible(NULL)
}
