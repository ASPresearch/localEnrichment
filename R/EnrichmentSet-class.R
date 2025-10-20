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
