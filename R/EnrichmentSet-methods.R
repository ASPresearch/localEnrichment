#' @rdname EnrichmentSet
#' @aliases show,EnrichmentSet-method
#' @param object An `EnrichmentSet` object.
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
#' @param object An `EnrichmentSet` object.
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
#' @name EnrichmentSet_as_list
#' @description Converts an EnrichmentSet object to a named list of feature vectors,
#' suitable for enrichment analysis.
setAs("EnrichmentSet", "list", function(from) {
  data <- from@data
  setNames(
    lapply(strsplit(data$feature_ids, ";"), trimws),
    data$set_name
  )
})

#' @title Coerce EnrichmentSet to data.frame
#' @name EnrichmentSet_as_data_frame
#' @description Converts an EnrichmentSet object to a long-format data.frame
#'   with one row per (set_name, feature_id) pair.
setAs("EnrichmentSet", "data.frame", function(from) {
  df <- from@data
  expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    tibble::tibble(
      set_id   = df$set_id[i],
      set_name = df$set_name[i],
      feature_id = unlist(strsplit(df$feature_ids[i], ";", fixed = TRUE))
    )
  }))
  if (!is.null(from@metadata$set_source)) expanded$source <- from@metadata$set_source
  if (!is.null(from@metadata$feature_species)) expanded$species <- from@metadata$feature_species
  expanded
})

# Ensure S4 and S3 coercion are consistent

setMethod("as.data.frame", "EnrichmentSet", function(x, ...) as(x, "data.frame"))

#' @title Coerce EnrichmentSet to data.frame (S3 wrapper)
#' @description Wrapper that allows `as.data.frame()` to be called directly
#' on `EnrichmentSet` objects.
#' @param x An `EnrichmentSet` object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with one row per (set_name, feature_id) pair.
#' @method as.data.frame EnrichmentSet
#' @export
as.data.frame.EnrichmentSet <- function(x, ...) {
  as(x, "data.frame")
}


#' @title Coerce EnrichmentSet to a condensed metabolite-set table
#' @name EnrichmentSet_as_MetaboliteSetDataFrame
#' @description Converts an EnrichmentSet object to a two-column data.frame
#'   with one row per set and a comma-separated string of metabolite IDs.
#' @param x An `EnrichmentSet` object.
#' @param id_sep Separator used to concatenate feature IDs (default = ",").
#' @return A data.frame with columns `set_name` and `Metabolites`.
#' @export
#' @export
as.MetaboliteSetDataFrame <- function(x, id_sep = ",") {
  stopifnot(inherits(x, "EnrichmentSet"))

  df <- x@data
  tibble::tibble(
    set_name = df$set_name,
    Metabolites = vapply(
      strsplit(df$feature_ids, ";", fixed = TRUE),
      function(v) paste(unique(v), collapse = id_sep),
      character(1)
    )
  )
}

