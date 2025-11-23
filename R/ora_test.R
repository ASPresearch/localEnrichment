#' Over-Representation Analysis (ORA)
#'
#' Perform an over-representation test to identify sets
#' (pathways, chemical classes, ontological terms, etc.)
#' that are significantly enriched in a list of selected features.
#'
#' @param eset An \code{EnrichmentSet} object containing the mapping
#'   between sets and features.
#' @param selected A character vector with the identifiers of selected features
#'   (e.g. significant metabolites or genes).
#' @param background Optional character vector of all possible feature identifiers.
#'   If \code{NULL}, all unique features present in \code{eset} are used.
#' @param test Character string specifying the test type:
#'   \code{"fisher"} (default) or \code{"hypergeom"}.
#' @param p_adjust Method for multiple testing correction
#'   (see \code{p.adjust.methods}); default is \code{"BH"}.
#' @param min_set_size Minimum number of features required in a set
#'   (default = 3).
#' @param max_set_size Maximum number of features allowed in a set
#'   (default = 500).
#' @param quiet Logical; if \code{FALSE}, progress messages are printed.
#'
#' @return A tibble with one row per set, containing:
#'   \itemize{
#'     \item \code{set_name}: set name.
#'     \item \code{n_selected_in_set}: number of selected features in the set.
#'     \item \code{n_in_set}: total number of features in the set.
#'     \item \code{p_value}: raw p-value.
#'     \item \code{p_adj}: adjusted p-value.
#'     \item \code{fold_enrichment}: enrichment ratio.
#'     \item \code{overlap_features}: semicolon-separated list of overlapping IDs.
#'   }
#' @export
#' @examples
#' meta <- list(
#'   mapping_name = "toy_example",
#'   feature_id_type = "HMDB_ID",
#'   feature_species = "Homo sapiens",
#'   set_source = "ExampleDB",
#'   version = "v1",
#'   description = "Toy example for ORA"
#' )
#'
#' df <- data.frame(
#'   set_id = c("A", "B"),
#'   set_name = c("Pathway A", "Pathway B"),
#'   feature_ids = c("x1;x2;x3;x4", "x3;x4;x5;x6")
#' )
#'
#' eset <- EnrichmentSet(df, meta)
#' selected <- c("x1", "x3", "x6")
#'
#' ora_test(eset, selected)
ora_test <- function(eset, selected,
                     background = NULL,
                     test = "fisher",
                     p_adjust = "BH",
                     min_set_size = 3,
                     max_set_size = 500,
                     quiet = TRUE) {

  #--- 1. Validation ---------------------------------------------------------
  if (!inherits(eset, "EnrichmentSet")) {
    stop("Input 'eset' must be an EnrichmentSet object.")
  }
  if (!is.character(selected)) {
    stop("'selected' must be a character vector of feature identifiers.")
  }
  mapping <- as(eset, "list")
  if (is.null(background)) {
    background <- unique(unlist(mapping))
  }

  # remove NAs
  selected <- selected[!is.na(selected)]
  background <- background[!is.na(background)]

  N <- length(background)
  K <- length(unique(selected))

  if (!quiet) message("Background features: ", N, " | Selected: ", K)

  #--- 2. Iterate over sets --------------------------------------------------
  results <- purrr::map_dfr(names(mapping), function(sid) {

    feats <- unique(mapping[[sid]])
    feats <- feats[feats %in% background]
    m <- length(feats)
    if (m < min_set_size || m > max_set_size) return(NULL)

    overlap <- intersect(selected, feats)
    k <- length(overlap)
    if (k == 0) return(NULL)

    if (test == "fisher") {
      p <- fisher.test(matrix(c(k, m - k, K - k, N - m - K + k), nrow = 2))$p.value
    } else if (test == "hypergeom") {
      p <- phyper(k - 1, m, N - m, K, lower.tail = FALSE)
    } else {
      stop("Unsupported test type: use 'fisher' or 'hypergeom'.")
    }

    tibble::tibble(
      set_id = sid,
      n_selected_in_set = k,
      n_in_set = m,
      p_value = p,
      overlap_features = paste(overlap, collapse = ";")
    )
  })

  if (nrow(results) == 0) {
    warning("No sets passed the size thresholds or had overlaps.")
    return(results)
  }

  # Afegim el set_name "en cristià" des de l'EnrichmentSet
  set_info <- eset@data %>%
    dplyr::distinct(.data$set_id, .data$set_name)

  results <- results %>%
    dplyr::left_join(set_info, by = "set_id") %>%
    dplyr::relocate(set_id, set_name, .before = n_selected_in_set) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = p_adjust),
      fold_enrichment = (n_selected_in_set / length(selected)) /
        (n_in_set / length(background))
    ) %>%
    dplyr::arrange(p_value)

  attr(results, "metadata") <- eset@metadata
  class(results) <- c("EnrichmentResult", class(results))
  results
}
