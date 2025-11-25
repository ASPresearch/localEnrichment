`%||%` <- function(a,b) if(!is.null(a)) a else b

#' @title Convert a list-style pathway object to EnrichmentSet
#' @description Converts a list object (like MetaboAnalystR SMPDB/KEGG data)
#'   into a formal `EnrichmentSet` with standardized slots.
#' @param x A list with elements `name`, `description`, `version`, and `sets`,
#'   where `sets` is a named list of character vectors (IDs).
#' @param id_type Character string describing the type of feature IDs (default = "HMDB_ID").
#' @param species Species name (default = "Homo sapiens").
#' @return An object of class \code{EnrichmentSet}.
#' @export
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

  df <- Eset@data

  # Ens quedem els sets on hi ha coincidències
  keep <- vapply(df$feature_list, function(v) any(v %in% ids), logical(1))
  df2 <- df[keep, ]

  message(sum(keep), " sets retained out of ", nrow(df))

  # Capar els feature_list al conjunt d'interès
  df2$feature_list <- lapply(df2$feature_list, function(v) intersect(v, ids))

  # Regenerar feature_ids
  df2$feature_ids <- vapply(
    df2$feature_list,
    function(v) paste(unique(v), collapse = ";"),
    character(1)
  )

  EnrichmentSet(df2, Eset@metadata)
}



