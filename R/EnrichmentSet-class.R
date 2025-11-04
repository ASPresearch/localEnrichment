#' @importFrom utils globalVariables
utils::globalVariables(c("set_id", "set_name", "feature_ids"))

#' @importFrom methods as setAs
NULL

#' EnrichmentSet Class
#'
#' An S4 class representing a collection of sets (pathways, classes, ontologies, etc.)
#' together with their metadata and associated features.
#'
#' @slot metadata A named list with mapping information:
#'   \itemize{
#'     \item mapping_name: internal name for the mapping
#'     \item feature_id_type: type of feature IDs (e.g., HMDB_ID, ENSEMBL)
#'     \item feature_species: species if applicable
#'     \item set_source: database or source of sets
#'     \item version: version or date of database
#'     \item description: description of the mapping
#'   }
#' @slot data A data.frame containing at least:
#'   \itemize{
#'     \item set_id: unique identifier for each set
#'     \item set_name: human-readable name
#'     \item feature_ids: semicolon-separated string of feature IDs
#'   }
#'
#' @export
setClass(
  "EnrichmentSet",
  slots = c(
    metadata = "list",
    data = "data.frame"
  ),
  prototype = list(
    metadata = list(),
    data = data.frame()
  ),
  validity = function(object) {
    required <- c("set_id", "set_name", "feature_ids")
    missing <- setdiff(required, colnames(object@data))
    if (length(missing) > 0) {
      return(paste("Missing required columns:", paste(missing, collapse = ", ")))
    }
    TRUE
  }
)

#' Constructor for EnrichmentSet
#'
#' @param data A data.frame containing at least columns set_id, set_name, feature_ids.
#' @param metadata A named list with mapping metadata.
#' @param sep Character separator used in feature_ids (default = ";").
#'
#' @return An object of class EnrichmentSet.
#' @export
EnrichmentSet <- function(data, metadata, sep = ";") {
  required_cols <- c("set_id", "set_name", "feature_ids")
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  data$feature_ids <- gsub("\\s+", "", data$feature_ids)
  data$feature_list <- strsplit(data$feature_ids, sep)
  new("EnrichmentSet", metadata = metadata, data = data)
}

#' @title Build an EnrichmentSet from a mapping table
#' @description High-level constructor that takes a two-column mapping (ID ↔ Category)
#'   and metadata, and returns a valid `EnrichmentSet` object.
#'
#' @param data A data.frame containing at least `id_col` and `category_col`
#' @param id_col Column name with metabolite or feature IDs (e.g., "HMDB")
#' @param category_col Column name with categories or set names (e.g., "Pathway")
#' @param set_name (optional) A name for the mapping, used as metadata$mapping_name
#' @param source (optional) Source database (used as metadata$set_source)
#' @param species (optional) Species name (used as metadata$feature_species)
#' @param version (optional) Version string or date
#' @param description (optional) Description of the mapping
#' @param sep (optional) Separator for feature_ids; defaults to ";"
#'
#' @importFrom dplyr select all_of filter distinct group_by summarise mutate
#' @importFrom tibble tibble
#' @importFrom stats setNames
#'
#' @return An S4 object of class `EnrichmentSet`
#' @export
buildEnrichmentSet <- function(
    data,
    id_col,
    category_col,
    set_name = NULL,
    source = NULL,
    species = "Homo sapiens",
    version = NULL,
    description = NULL,
    sep = ";"
) {
  if (!id_col %in% colnames(data))
    stop(paste("Column", id_col, "not found in data."))
  if (!category_col %in% colnames(data))
    stop(paste("Column", category_col, "not found in data."))

  # --- Clean data and remove missing values ----------------------------
  df <- data %>%
    select(all_of(c(id_col, category_col))) %>%
    filter(!is.na(.data[[id_col]]), !is.na(.data[[category_col]])) %>%
    distinct()

  # --- Build canonical data.frame: one row per category ----------------
  mapping <- df %>%
    group_by(.data[[category_col]]) %>%
    summarise(
      feature_ids = paste(unique(.data[[id_col]]), collapse = sep),
      .groups = "drop"
    ) %>%
    mutate(
      set_id = make.unique(as.character(.data[[category_col]])),
      set_name = .data[[category_col]]
    ) %>%
    select(set_id, set_name, feature_ids)

  # --- Build metadata list --------------------------------------------
  metadata <- list(
    mapping_name = ifelse(is.null(set_name), category_col, set_name),
    feature_id_type = id_col,
    feature_species = species,
    set_source = source,
    version = version,
    description = description
  )

  # --- Call the low-level constructor ---------------------------------
  EnrichmentSet(mapping, metadata, sep = sep)
}

#' Validate EnrichmentSet structure
#'
#' Checks that an object is a valid `EnrichmentSet` and contains the required columns.
#' @param Eset An object of class `EnrichmentSet`.
#' @return Invisibly returns TRUE if valid, otherwise throws an error.
#' @export
validateEnrichmentSet <- function(Eset) {
  stopifnot(inherits(Eset, "EnrichmentSet"))
  required_cols <- c("set_id", "set_name", "feature_ids")
  if (!all(required_cols %in% colnames(Eset@data))) {
    stop("Invalid EnrichmentSet: missing columns ",
         paste(setdiff(required_cols, colnames(Eset@data)), collapse = ", "))
  }
  invisible(TRUE)
}
