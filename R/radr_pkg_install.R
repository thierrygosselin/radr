#' Legacy radr dependency helper
#'
#' `radr_pkg_install()` no longer installs or updates packages. Automatic
#' installation made library changes difficult to anticipate and obscured the
#' distinction between required and optional components. Use
#' [radr_dependencies()] for a read-only diagnostic and follow the installation
#' instructions in the radr README.
#'
#' @param check Retained for compatibility and ignored.
#' Default: \code{check = TRUE}.
#' @param minimal.install Retained for compatibility and ignored.
#' Default: \code{minimal.install = FALSE}.
#'
#' @return The dependency table returned by [radr_dependencies()].
#' @export
radr_pkg_install <- function(check = TRUE, minimal.install = FALSE) {
  .Deprecated("radr_dependencies", package = "radr")
  radr_dependencies(verbose = isTRUE(check))
}
