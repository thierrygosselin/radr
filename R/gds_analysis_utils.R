.analysis_gds_start <- function(
    data, function.name, write.files = TRUE, verbose = TRUE, ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = function.name, verbose = verbose
  )
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), c("path.folder", "internal"))
  if (length(unknown)) {
    rlang::abort(paste0(
      "Unknown argument(s): ", paste(unknown, collapse = ", "), "."
    ))
  }
  .paralog_check_flag(write.files, "write.files")
  .paralog_check_flag(verbose, "verbose")
  parent <- dots$path.folder %||% getwd()
  internal <- isTRUE(dots$internal)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(parent, mustWork = TRUE)
  path.folder <- if (internal || !write.files) parent else radr_folder(
    rad.folder = paste0(function.name, "_", .start$file.date),
    path.folder = parent, prefix.int = TRUE
  )
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  opened <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    gds <- SeqArray::seqOpen(data)
    opened <- TRUE
  } else {
    rlang::abort(
      "`data` must be a GDS filepath or an open SeqVarGDSClass object."
    )
  }
  selection <- SeqArray::seqGetFilter(gds)
  SeqArray::seqFilterPush(gds)
  context <- new.env(parent = emptyenv())
  context$gds <- gds
  context$opened <- opened
  context$pushed <- TRUE
  context$restored <- FALSE
  context$selection <- selection
  context$startup <- .start
  context$path.folder <- path.folder
  context$verbose <- verbose
  context
}

.analysis_gds_restore <- function(context) {
  if (!isTRUE(context$pushed)) return(isTRUE(context$restored))
  SeqArray::seqFilterPop(context$gds)
  context$pushed <- FALSE
  after <- SeqArray::seqGetFilter(context$gds)
  context$restored <- identical(
    context$selection$sample.sel, after$sample.sel
  ) && identical(context$selection$variant.sel, after$variant.sel)
  if (!context$restored) {
    rlang::abort("The active GDS selection was not restored.")
  }
  TRUE
}

.analysis_gds_finish <- function(context) {
  if (isTRUE(context$pushed)) {
    try(.analysis_gds_restore(context), silent = TRUE)
  }
  if (isTRUE(context$opened)) {
    try(SeqArray::seqClose(context$gds), silent = TRUE)
  }
  tgbase::teardown(context$startup)
}

.analysis_metadata <- function(
    gds, strata = NULL, group.column = "STRATA", by.strata = TRUE
) {
  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  metadata <- .paralog_read_metadata(
    strata, gds, sample.id, group.column, by.strata
  )
  sample.id <- sample.id[sample.id %in% metadata$INDIVIDUALS]
  metadata <- metadata[
    match(sample.id, metadata$INDIVIDUALS), , drop = FALSE
  ]
  groups <- as.character(metadata[[group.column]])
  if (anyNA(groups) || any(!nzchar(trimws(groups)))) {
    rlang::abort(paste0("`", group.column, "` contains missing values."))
  }
  SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
  list(sample.id = sample.id, metadata = metadata, groups = groups)
}

.filter_gds_open <- function(data) {
  if (inherits(data, "SeqVarGDSClass")) return(list(gds = data, opened = FALSE))
  if (is.character(data) && length(data) == 1L && !is.na(data) &&
      file.exists(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    return(list(
      gds = SeqArray::seqOpen(data, readonly = FALSE), opened = TRUE
    ))
  }
  rlang::abort("`data` must be a GDS filepath or open SeqVarGDSClass object.")
}

.filter_gds_apply_markers <- function(
    gds, remove.variant.id, filter.label, path.folder, file.date,
    parameter.name, parameter.value, verbose = TRUE
) {
  parameters <- filter_parameters(
    generate = TRUE, initiate = TRUE, update = FALSE,
    parameter.obj = NULL, data = gds, path.folder = path.folder,
    file.date = file.date, internal = FALSE, verbose = verbose
  )
  metadata <- genometranslator::extract_markers_metadata(
    gds = gds, whitelist = FALSE
  )
  remove.variant.id <- unique(remove.variant.id)
  metadata <- metadata |>
    dplyr::mutate(
      FILTERS = dplyr::if_else(
        .data$VARIANT_ID %in% remove.variant.id &
          (is.na(.data$FILTERS) | .data$FILTERS == "whitelist"),
        filter.label, .data$FILTERS
      )
    )
  genometranslator::update_genome_gds(
    gds = gds, node.name = "markers.meta", value = metadata, sync = TRUE
  )
  filter_parameters(
    generate = FALSE, initiate = FALSE, update = TRUE,
    parameter.obj = parameters, data = gds,
    filter.name = filter.label, param.name = parameter.name,
    values = parameter.value, path.folder = path.folder,
    file.date = file.date, internal = FALSE, verbose = verbose
  )
  invisible(metadata)
}

.filter_gds_apply <- function(
    gds, remove.variant.id = integer(), remove.sample.id = character(),
    filter.label, path.folder, file.date, parameter.name,
    parameter.value, verbose = TRUE
) {
  parameters <- filter_parameters(
    generate = TRUE, initiate = TRUE, update = FALSE,
    parameter.obj = NULL, data = gds, path.folder = path.folder,
    file.date = file.date, internal = FALSE, verbose = verbose
  )
  markers <- genometranslator::extract_markers_metadata(gds, whitelist = FALSE)
  individuals <- genometranslator::extract_individuals_metadata(
    gds, whitelist = FALSE
  )
  markers <- markers |>
    dplyr::mutate(FILTERS = dplyr::if_else(
      .data$VARIANT_ID %in% remove.variant.id &
        (is.na(.data$FILTERS) | .data$FILTERS == "whitelist"),
      filter.label, .data$FILTERS
    ))
  individuals <- individuals |>
    dplyr::mutate(FILTERS = dplyr::if_else(
      as.character(.data$INDIVIDUALS) %in% remove.sample.id &
        (is.na(.data$FILTERS) | .data$FILTERS == "whitelist"),
      filter.label, .data$FILTERS
    ))
  genometranslator::update_genome_gds(
    gds = gds, node.name = "markers.meta", value = markers, sync = FALSE
  )
  genometranslator::update_genome_gds(
    gds = gds, node.name = "individuals.meta", value = individuals,
    sync = TRUE
  )
  filter_parameters(
    generate = FALSE, initiate = FALSE, update = TRUE,
    parameter.obj = parameters, data = gds,
    filter.name = filter.label, param.name = parameter.name,
    values = parameter.value, path.folder = path.folder,
    file.date = file.date, internal = FALSE, verbose = verbose
  )
  invisible(list(markers = markers, individuals = individuals))
}
