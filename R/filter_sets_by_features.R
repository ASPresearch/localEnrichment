#' Filter an EnrichmentSet by feature overlap
#'
#' Keep only sets with at least `min_overlap` features present in `features`.
#'
#' @param eset An \code{EnrichmentSet} object.
#' @param features Character vector of feature IDs to check.
#' @param min_overlap Minimum number of overlapping features (default = 1).
#' @return A reduced \code{EnrichmentSet}.
#' @export
filter_sets_by_features <- function(eset, features, min_overlap = 1) {
  if (!inherits(eset, "EnrichmentSet"))
    stop("'eset' must be an EnrichmentSet.")
  if (!is.character(features))
    stop("'features' must be a character vector.")

  df <- eset@data
  ov <- vapply(df$feature_list, function(x) sum(x %in% features), numeric(1))
  keep <- ov >= min_overlap
  kept <- df[keep, , drop = FALSE]

  if (nrow(kept) == 0)
    warning("No sets retained after filtering (min_overlap = ", min_overlap, ").")

  new("EnrichmentSet", metadata = eset@metadata, data = kept)
}
