`%||%` <- function(a,b) if(!is.null(a)) a else b

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
  sets_filtered <- lapply(sets_filtered, function(v) intersect(v, ids))

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

  EnrichmentSet(df, Eset@metadata)
}

