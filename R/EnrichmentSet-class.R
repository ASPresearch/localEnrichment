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

#' Validate enrichment mapping (low-level)
#'
#' Ensures that the mapping table contains required columns and
#' normalizes the feature_ids column by removing whitespace.
#'
#' @param data data.frame with columns set_id, set_name, feature_ids
#' @param metadata list with mapping metadata
#'
#' @return validated and normalized mapping table
#' @keywords internal
validateEnrichment <- function(data, metadata) {

  required_cols <- c("set_id", "set_name", "feature_ids")
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Invalid enrichment mapping: missing columns ",
         paste(missing, collapse = ", "))
  }

  # normalize feature_ids
  data$feature_ids <- gsub("\\s+", "", data$feature_ids)

  # return data frame
  data
}


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

  # Normalize feature_ids (remove extra spaces, fix separators)
  data$feature_ids <- gsub("\\s+", "", data$feature_ids)
  data$feature_ids <- gsub(paste0(sep, "$"), "", data$feature_ids)  # remove trailing sep

  # Build feature_list (primary internal structure)
  data$feature_list <- lapply(strsplit(data$feature_ids, sep, fixed = TRUE), function(v) {
    v <- trimws(v)          # remove whitespace
    v <- v[v != ""]         # remove empty entries
    v <- unique(v)          # ensure uniqueness
    v
  })

  # Ensure set_id is unique
  if (anyDuplicated(data$set_id)) {
    warning("Duplicated set_id detected applying make.unique().")
    data$set_id <- make.unique(data$set_id)
  }

  new("EnrichmentSet", metadata = metadata, data = data)
}


#' Build an EnrichmentSet object from a mapping table
#'
#' High-level constructor for generating an EnrichmentSet from a mapping table.
#'
#' @param data A data.frame containing mapping information.
#' @param id_col Column name containing feature IDs (e.g., HMDB, PubChem).
#' @param category_col Column name with human-readable set labels (e.g., pathways).
#' @param set_id_col Optional column with stable set identifiers (e.g., KEGG ID, SMPDB ID).
#' @param set_name Optional name of the mapping (stored in metadata).
#' @param source Optional source database name.
#' @param species Species name, e.g. "Homo sapiens".
#' @param version Optional version string or date.
#' @param description Optional description of the enrichment mapping.
#' @param sep Separator used to collapse multiple features in `feature_ids`.
#'
#' @return An object of class \code{EnrichmentSet}.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   HMDB = c("H1","H2","H3"),
#'   Pathway = c("P1","P1","P2"),
#'   PathwayID = c("PW1","PW1","PW2")
#' )
#'
#' eset <- buildEnrichmentSet(
#'   data = df,
#'   id_col = "HMDB",
#'   category_col = "Pathway",
#'   set_id_col = "PathwayID"
#' )
buildEnrichmentSet <- function(
    data,
    id_col,
    category_col,
    set_id_col = NULL,
    set_name = NULL,
    source = NULL,
    species = "Homo sapiens",
    version = NULL,
    description = NULL,
    sep = ";"
) {
  stopifnot(id_col %in% colnames(data))
  stopifnot(category_col %in% colnames(data))
  if (!is.null(set_id_col)) stopifnot(set_id_col %in% colnames(data))

  # 1) Netejar i deduplicar
  df <- data %>%
    dplyr::filter(
      !is.na(.data[[id_col]]),
      !is.na(.data[[category_col]])
    )

  if (!is.null(set_id_col)) {
    df <- df %>%
      dplyr::filter(!is.na(.data[[set_id_col]])) %>%
      dplyr::distinct(.data[[id_col]], .data[[category_col]], .data[[set_id_col]])
  } else {
    df <- df %>%
      dplyr::distinct(.data[[id_col]], .data[[category_col]])
  }

  # 2) Construir mapping: set_id, set_name, feature_ids
  if (is.null(set_id_col)) {
    mapping <- df %>%
      dplyr::group_by(.data[[category_col]]) %>%
      dplyr::summarise(
        feature_ids = paste(unique(.data[[id_col]]), collapse = sep),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        set_id   = make.unique(as.character(.data[[category_col]])),
        set_name = as.character(.data[[category_col]])
      ) %>%
      dplyr::select(set_id, set_name, feature_ids)
  } else {
    mapping <- df %>%
      dplyr::group_by(.data[[set_id_col]], .data[[category_col]]) %>%
      dplyr::summarise(
        feature_ids = paste(unique(.data[[id_col]]), collapse = sep),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        set_id   = as.character(.data[[set_id_col]]),
        set_name = as.character(.data[[category_col]])
      ) %>%
      dplyr::select(set_id, set_name, feature_ids)
  }

  # 3) Metadata
  metadata <- list(
    mapping_name = ifelse(is.null(set_name), category_col, set_name),
    feature_id_type = id_col,
    feature_species = species,
    set_source = source,
    version = version,
    description = description
  )

  # 4) Validacio del mapping (funcio antiga del paquet)
  mapping <- validateEnrichment(mapping, metadata)

  # 5) Crear objecte S4
  Eset <- EnrichmentSet(mapping, metadata, sep = sep)

  # 6) Validacio de l’objecte S4 a nivell alt
  validateEnrichmentSet(Eset)

  return(Eset)
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

