#' Exclude a chromosome or genomic interval from a GDS
#'
#' Blacklist active, whitelisted markers on a chromosome or linkage group.
#' The filter records `filter.genomic.range` in the genometranslator
#' `markers.meta` node and updates the active variant selection. Genotypes
#' remain in the file. Use a separate working copy to retain both panels.
#'
#' @param data A writable, open genometranslator GDS (`SeqVarGDSClass`).
#' @param chromosome One chromosome name matching the GDS exactly.
#' @param start,end Optional first and last base-pair positions (1-based,
#'   both inclusive). Omit both to exclude the entire chromosome.
#' @param verbose Print counts of newly excluded and remaining active markers.
#' @return The same open GDS object, modified in place.
#' @examples
#' \dontrun{
#' # Use the candidate coordinates from detect_inversions().
#' region <- scan$candidates[1, ]
#' gds <- filter_genomic_range(gds, region$chromosome, region$start, region$end)
#' # Whole linkage group, using its exact identifier in the GDS:
#' gds <- filter_genomic_range(gds, chromosome = "LG14")
#' }
#' @export
filter_genomic_range <- function(data, chromosome, start = NULL, end = NULL,
                                 verbose = TRUE) {
  if (!inherits(data, "SeqVarGDSClass")) {
    stop("Supply an open genometranslator GDS.", call. = FALSE)
  }
  if (length(chromosome) != 1L || is.na(chromosome) ||
      !is.atomic(chromosome) || !nzchar(as.character(chromosome))) {
    stop("Supply one chromosome identifier.", call. = FALSE)
  }
  chromosome <- as.character(chromosome)
  if (xor(is.null(start), is.null(end))) {
    stop("Supply both start and end, or omit both.", call. = FALSE)
  }
  valid.bp <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!is.null(start) && (!valid.bp(start) || !valid.bp(end) || start > end)) {
    stop("Expected whole base-pair coordinates with 1 <= start <= end.", call. = FALSE)
  }

  # Read unfiltered coordinate nodes so a repeated exclusion can recognise
  # the chromosome even when no markers on it are still selected.
  read.node <- function(name) gdsfmt::read.gdsn(gdsfmt::index.gdsn(data, name))
  chromosomes <- as.character(read.node("chromosome"))
  on.chromosome <- which(chromosomes == chromosome)
  if (!length(on.chromosome)) stop("Unknown chromosome: ", chromosome, call. = FALSE)
  if (!is.null(start)) {
    position <- suppressWarnings(as.numeric(read.node("position")[on.chromosome]))
    if (any(!is.finite(position) | position < 1 | position != floor(position))) {
      stop("Invalid marker positions on chromosome ", chromosome, ".", call. = FALSE)
    }
    on.chromosome <- on.chromosome[position >= start & position <= end]
  }
  if (is.null(gdsfmt::index.gdsn(data, "genometranslator/markers.meta", silent = TRUE))) {
    stop("Missing genometranslator/markers.meta.", call. = FALSE)
  }
  markers <- genometranslator::extract_markers_metadata(data, whitelist = FALSE)
  ids <- read.node("variant.id")
  if (anyDuplicated(markers$VARIANT_ID) || !setequal(ids, markers$VARIANT_ID)) {
    stop("Marker metadata and GDS variant IDs disagree.", call. = FALSE)
  }
  active <- SeqArray::seqGetData(data, "variant.id")
  exclude <- markers$FILTERS %in% "whitelist" &
    markers$VARIANT_ID %in% intersect(active, ids[on.chromosome])
  removed <- sum(exclude)
  retained <- active[active %in% markers$VARIANT_ID[markers$FILTERS %in% "whitelist" & !exclude]]
  if (removed && !length(retained)) {
    stop("This filter would exclude every active marker.", call. = FALSE)
  }
  if (removed) {
    markers$FILTERS <- as.character(markers$FILTERS)
    markers$FILTERS[exclude] <- "filter.genomic.range"
    genometranslator::update_genome_gds(
      gds = data, node.name = "markers.meta", value = markers, sync = FALSE
    )
  }
  if (!identical(active, retained)) {
    SeqArray::seqSetFilter(data, variant.id = retained, verbose = FALSE)
  }
  if (verbose) message("Excluded ", removed, " marker(s); ", length(retained),
                       " active marker(s) remain.")
  data
}
