#' Check radr dependencies
#'
#' Reports the required and optional R packages used by radr and checks whether
#' the optional \command{bcftools} and BayeScan executables are available. The
#' function is diagnostic: it does not install packages or modify a Conda
#' environment.
#'
#' @param verbose Logical. Print guidance for unavailable components.
#' Default: \code{verbose = TRUE}.
#'
#' @return A tibble with the component, source, requirement level, principal
#' use, and availability.
#' @export
radr_dependencies <- function(verbose = TRUE) {
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    rlang::abort("`verbose` must be TRUE or FALSE.")
  }

  required.cran <- c(
    "arrow", "carrier", "cli", "ComplexUpset", "data.table", "dplyr",
    "ggplot2", "HardyWeinberg", "magrittr", "matrixStats", "purrr",
    "readr", "rlang", "scales", "stringi", "tibble", "tidyr",
    "tidyselect", "vctrs", "vroom"
  )
  required.bioc <- c("gdsfmt", "SeqArray")
  required.github <- c("genometranslator", "tgbase")

  optional <- tibble::tribble(
    ~component, ~source, ~used_by,
    "adegenet", "CRAN", "selected population-genetic objects and methods",
    "amap", "CRAN", "distance calculations on tidy data",
    "fst", "CRAN", "FST-backed data workflows",
    "ragg", "CRAN", "fast PNG output from detect_ibm()",
    "SNPRelate", "Bioconductor", "LD pruning, LD statistics, and GDS IBS",
    "stringdist", "CRAN", "sequence and identifier distance utilities"
  )

  packages <- tibble::tibble(
    component = c(required.cran, required.bioc, required.github),
    source = c(
      rep("CRAN", length(required.cran)),
      rep("Bioconductor", length(required.bioc)),
      "GitHub: thierrygosselin/genometranslator",
      "GitHub: thierrygosselin/tgbase"
    ),
    required = TRUE,
    used_by = "core radr workflows"
  ) |>
    dplyr::bind_rows(
      dplyr::mutate(optional, required = FALSE) |>
        dplyr::select(component, source, required, used_by)
    )

  packages$available <- vapply(
    packages$component,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  executables <- tibble::tibble(
    component = c("bcftools", "BayeScan"),
    source = c("Conda/Bioconda or system PATH", "Conda/Bioconda or explicit path"),
    required = FALSE,
    used_by = c("VCF-level filter_*_vcf() functions", "run_bayescan()"),
    available = c(
      nzchar(Sys.which("bcftools")),
      any(nzchar(Sys.which(c("bayescan", "bayescan_2.1"))))
    )
  )

  result <- dplyr::bind_rows(packages, executables)

  if (verbose) {
    missing.required <- dplyr::filter(result, required, !available)
    missing.optional <- dplyr::filter(result, !required, !available)

    if (nrow(missing.required) == 0L) {
      message("All required radr components are available.")
    } else {
      message("Missing required radr components:")
      apply(missing.required, 1L, function(x) {
        message("  ", x[["component"]], " [", x[["source"]], "]")
      })
    }

    if (nrow(missing.optional) > 0L) {
      message("Optional components not currently available:")
      apply(missing.optional, 1L, function(x) {
        message(
          "  ", x[["component"]], " [", x[["source"]], "] — ",
          x[["used_by"]]
        )
      })
    }
    message("See the Installation section in the radr README for installation guidance.")
  }

  result
}
