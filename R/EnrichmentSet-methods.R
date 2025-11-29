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
  sets <- lapply(strsplit(data$feature_ids, ";", fixed = TRUE), trimws)
  # clau interna = ID estable
  names(sets) <- data$set_id
  sets
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

#' @title Coerce EnrichmentSet to a condensed metabolite-set table
#' @name EnrichmentSet_as_MetaboliteSetDataFrame
#' @description Converts an EnrichmentSet object to a two- or three-column
#' data.frame with one row per set and a concatenated string of metabolites.
#' If `id_sep` is NULL, the output separator is automatically set to match
#' the one used in `feature_ids`.
#'
#' @param x An `EnrichmentSet` object.
#' @param id_type Which identifier to include?
#'   One of: "name" (default), "id", or "both".
#' @param id_sep Separator used to concatenate metabolite IDs.
#'   If NULL (default), the function uses the separator already present
#'   in `feature_ids` (either ";" or ","). If not NULL, it is used only
#'   as output separator, preserving the detected input separator.
#' @return A tibble with columns set_name / set_id and Metabolites.
#' @export
as.MetaboliteSetDataFrame <- function(x,
                                      id_type = c("name", "id", "both"),
                                      id_sep = NULL) {

  stopifnot(inherits(x, "EnrichmentSet"))
  id_type <- match.arg(id_type)

  df <- x@data

  ## --- Detectar separador d'ENTRADA (sep_in) -------------------------
  # Agafem el primer feature_ids no buit
  sample_str <- df$feature_ids[which(nzchar(df$feature_ids))[1]]

  if (grepl(";", sample_str, fixed = TRUE)) {
    sep_in <- ";"
  } else if (grepl(",", sample_str, fixed = TRUE)) {
    sep_in <- ","
  } else {
    sep_in <- ";"  # fallback raonable
  }

  ## --- Decidir separador de SORTIDA (sep_out) ------------------------
  if (is.null(id_sep)) {
    sep_out <- sep_in          # conservar el que ja porta l'EnrichmentSet
  } else {
    sep_out <- id_sep          # forçar-ne un altre ("," o ";")
  }

  ## --- Construir la columna de metabòlits ----------------------------
  metabolites <- vapply(
    strsplit(df$feature_ids, split = sep_in, fixed = TRUE),
    function(v) paste(unique(trimws(v)), collapse = sep_out),
    character(1)
  )

  ## --- Construir sortida segons id_type ------------------------------
  out <- dplyr::tibble(Metabolites = metabolites)

  if (id_type %in% c("name", "both")) {
    out$set_name <- df$set_name
  }
  if (id_type %in% c("id", "both")) {
    out$set_id <- df$set_id
  }

  out <- dplyr::relocate(out, dplyr::starts_with("set_"), Metabolites)

  out
}
