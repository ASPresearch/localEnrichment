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
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.")
  library(dplyr)

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



#' @rdname EnrichmentSet
#' @aliases show,EnrichmentSet-method
#' @export
setMethod("show", "EnrichmentSet", function(object) {
  md <- object@metadata
  cat("EnrichmentSet:", md$mapping_name, "\n")
  cat("  Source:", md$set_source, "\n")
  cat("  Feature IDs:", md$feature_id_type, "\n")
  cat("  Number of sets:", nrow(object@data), "\n")
  cat("  Example set:", object@data$set_name[1], "\n")
})

#' @rdname EnrichmentSet
#' @aliases summary,EnrichmentSet-method
#' @export
setMethod("summary", "EnrichmentSet", function(object) {
  cat("EnrichmentSet summary:\n")
  cat("  Mapping name:", object@metadata$mapping_name, "\n")
  cat("  Source:", object@metadata$set_source, "\n")
  cat("  Feature IDs:", object@metadata$feature_id_type, "\n")
  cat("  Number of sets:", nrow(object@data), "\n")
  sizes <- lengths(object@data$feature_list)
  cat("  Mean set size:", mean(sizes), "features\n")
  cat("  Median set size:", median(sizes), "features\n")
})

#' @title Coerce EnrichmentSet to list
#' @description Converts an EnrichmentSet object to a named list of feature vectors,
#' suitable for enrichment analysis.
#' @name EnrichmentSet-coercion
#' @rdname EnrichmentSet
setAs("EnrichmentSet", "list", function(from) {
  data <- from@data
  sets_list <- setNames(
    lapply(strsplit(data$feature_ids, ";"), trimws),
    data$set_name
  )
  sets_list
})

#' @title Coerce EnrichmentSet to data.frame
#' @description Converts an EnrichmentSet object to a long-format data.frame
#'   with one row per (set_name, feature_id) pair.
#' @name EnrichmentSet-coercion
#' @rdname EnrichmentSet
#' @export
setAs("EnrichmentSet", "data.frame", function(from) {
  df <- from@data
  # Expand feature_ids into one row per element
  expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    tibble::tibble(
      set_id   = df$set_id[i],
      set_name = df$set_name[i],
      feature_id = unlist(strsplit(df$feature_ids[i], ";", fixed = TRUE)),
      stringsAsFactors = FALSE
    )
  }))

  # Add metadata columns if desired
  if (!is.null(from@metadata$set_source)) expanded$source <- from@metadata$set_source
  if (!is.null(from@metadata$feature_species)) expanded$species <- from@metadata$feature_species

  expanded
})

#' @export
as.data.frame.EnrichmentSet <- function(x, ...) {
  as(x, "data.frame")
}

# Quick helper to build an EnrichmentSet from a legacy list object
# (e.g. MetaboAnalyst KEGG/SMPDB data)
as.EnrichmentSet.list <- function(x, id_type = "HMDB_ID", species = "Homo sapiens") {
  if (!is.list(x) || !"sets" %in% names(x))
    stop("Input must be a list with a 'sets' element.")
  df <- tibble::tibble(
    set_id = make.unique(names(x$sets)),
    set_name = names(x$sets),
    feature_ids = vapply(x$sets, function(v) paste(unique(v), collapse = ";"), character(1))
  )
  meta <- list(
    mapping_name = x$name %||% "UnknownMapping",
    feature_id_type = id_type,
    feature_species = species,
    set_source = x$name %||% "UnknownSource",
    version = x$version %||% NA,
    description = x$description %||% NA
  )
  EnrichmentSet(df, meta)
}
`%||%` <- function(a,b) if(!is.null(a)) a else b


#' @title Filter EnrichmentSet by feature IDs
#' @description Returns a new EnrichmentSet containing only sets that include
#'   at least one of the provided feature IDs.
#' @param Eset An object of class `EnrichmentSet`
#' @param ids Character vector of feature IDs to retain
#' @return A filtered `EnrichmentSet` object
#' @export
filterEnrichmentSet <- function(Eset, ids) {
  stopifnot(inherits(Eset, "EnrichmentSet"))
  if (missing(ids) || length(ids) == 0)
    stop("You must provide a non-empty vector of IDs to filter by.")

  sets_list <- as(Eset, "list")
  sets_filtered <- Filter(function(v) any(v %in% ids), sets_list)

  message(length(sets_filtered), " sets retained out of ", length(sets_list))

  df <- tibble::tibble(
    set_id = make.unique(names(sets_filtered)),
    set_name = names(sets_filtered),
    feature_ids = vapply(
      sets_filtered,
      function(v) paste(unique(v), collapse = ";"),
      character(1)
    )
  )

  meta <- Eset@metadata
  EnrichmentSet(df, meta)
}

