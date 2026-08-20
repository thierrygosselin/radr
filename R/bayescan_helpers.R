# BayeScan executable helpers --------------------------------------------------

find_conda_executable <- function(conda = NULL) {
  candidates <- c(
    conda,
    Sys.which("micromamba"),
    Sys.which("mamba"),
    Sys.which("conda"),
    file.path(path.expand("~"), "mambaforge", "bin", "mamba"),
    file.path(path.expand("~"), "mambaforge", "bin", "conda"),
    file.path(path.expand("~"), "miniforge3", "bin", "mamba"),
    file.path(path.expand("~"), "miniforge3", "bin", "conda")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1]], mustWork = TRUE)
}

locate_bayescan_in_conda <- function(conda.env, conda = NULL) {
  if (is.null(conda.env) || !nzchar(conda.env)) return(NULL)

  # An environment prefix can be supplied instead of its name.
  if (dir.exists(conda.env)) {
    candidates <- file.path(
      normalizePath(conda.env, mustWork = TRUE),
      "bin",
      c("bayescan", "bayescan_2.1")
    )
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates)) return(candidates[[1]])
  }

  conda <- find_conda_executable(conda)
  if (is.null(conda)) return(NULL)

  # Named environments normally live beside the Conda/Mamba installation.
  # Prefer direct lookup because `conda run` can require a writable lockfile.
  conda.base <- dirname(dirname(conda))
  candidates <- file.path(
    conda.base,
    "envs",
    conda.env,
    "bin",
    c("bayescan", "bayescan_2.1")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    return(normalizePath(candidates[[1]], mustWork = TRUE))
  }

  for (executable in c("bayescan", "bayescan_2.1")) {
    found <- suppressWarnings(system2(
      command = conda,
      args = c("run", "-n", conda.env, "which", executable),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(found, "status")
    if (is.null(status)) status <- 0L
    found <- found[file.exists(found)]
    if (identical(status, 0L) && length(found)) {
      return(normalizePath(found[[1]], mustWork = TRUE))
    }
  }
  NULL
}

#' Locate and validate BayeScan
#'
#' Searches an explicit path, the system \code{PATH}, and optionally a Conda
#' environment for a BayeScan 2.1 executable.
#'
#' @param bayescan.path Explicit path to a BayeScan executable. When
#' \code{NULL}, the executable is discovered automatically.
#' Default: \code{bayescan.path = NULL}.
#' @param conda.env Conda environment name or prefix containing BayeScan.
#' Default: \code{conda.env = "genomics"}.
#' @param conda Optional path to a Conda, Mamba, or Micromamba executable.
#' Default: \code{conda = NULL}.
#' @param verbose Logical. Display the detected executable.
#' Default: \code{verbose = TRUE}.
#'
#' @return The normalized BayeScan executable path.
#' @export
check_bayescan <- function(
  bayescan.path = NULL,
  conda.env = "genomics",
  conda = NULL,
  verbose = TRUE
) {
  candidates <- c(
    bayescan.path,
    unname(Sys.which("bayescan")),
    unname(Sys.which("bayescan_2.1"))
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates)) {
    executable <- normalizePath(candidates[[1]], mustWork = TRUE)
  } else {
    executable <- locate_bayescan_in_conda(
      conda.env = conda.env,
      conda = conda
    )
  }

  if (is.null(executable) || file.access(executable, mode = 1L) != 0L) {
    rlang::abort(paste0(
      "BayeScan was not found or is not executable. ",
      "Install it with radr::install_bayescan(), place it on PATH, ",
      "or supply bayescan.path."
    ))
  }

  help.output <- suppressWarnings(system2(
    command = executable,
    args = "-help",
    stdout = TRUE,
    stderr = TRUE
  ))
  if (!any(grepl("BayeScan 2.1", help.output, fixed = TRUE))) {
    rlang::abort(
      paste0("The detected executable is not a working BayeScan 2.1 program: ",
             executable)
    )
  }

  if (verbose) message("BayeScan executable: ", executable)
  executable
}

#' Install BayeScan with Conda or Mamba
#'
#' Creates a dedicated environment and installs the Bioconda BayeScan 2.1
#' package. Installation is performed only when this function is called
#' explicitly.
#'
#' @param conda.env Name of the environment to create.
#' Default: \code{conda.env = "genomics"}.
#' @param conda Optional path to a Conda, Mamba, or Micromamba executable.
#' Default: \code{conda = NULL}.
#' @param force Logical. Re-run the installation even when BayeScan is already
#' available in the environment. Default: \code{force = FALSE}.
#' @param verbose Logical. Display installation progress.
#' Default: \code{verbose = TRUE}.
#'
#' @return Invisibly returns the installed BayeScan executable path.
#' @export
install_bayescan <- function(
  conda.env = "genomics",
  conda = NULL,
  force = FALSE,
  verbose = TRUE
) {
  conda <- find_conda_executable(conda)
  if (is.null(conda)) {
    rlang::abort(
      "Conda, Mamba, or Micromamba was not found. Install Miniforge first."
    )
  }

  existing <- locate_bayescan_in_conda(conda.env = conda.env, conda = conda)
  if (!force && !is.null(existing)) {
    if (verbose) message("BayeScan is already installed: ", existing)
    return(invisible(existing))
  }

  conda.base <- dirname(dirname(conda))
  env.by.name <- file.path(conda.base, "envs", conda.env)
  env.exists <- dir.exists(conda.env) || dir.exists(env.by.name)
  env.option <- if (dir.exists(conda.env)) {
    c("--prefix", normalizePath(conda.env, mustWork = TRUE))
  } else {
    c("--name", conda.env)
  }
  action <- if (env.exists) "install" else "create"

  if (verbose) {
    action.message <- if (env.exists) "Installing into" else "Creating"
    message(action.message, " Conda environment: ", conda.env)
  }
  status <- system2(
    command = conda,
    args = c(
      action, "--yes", env.option,
      "--channel", "conda-forge", "--channel", "bioconda",
      "bayescan=2.1"
    )
  )
  if (!identical(status, 0L)) {
    rlang::abort("The Conda installation of BayeScan failed.")
  }

  executable <- check_bayescan(
    conda.env = conda.env,
    conda = conda,
    verbose = verbose
  )
  invisible(executable)
}
