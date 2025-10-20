#' Validate a mapping table for enrichment analysis
#'
#' @param mapping A data frame with at least three columns:
#'   `set_id`, `set_name`, and `feature_ids`.
#' @param sep Character used to separate feature IDs within each set (default ";").
#' @return A cleaned and validated mapping data frame.
#' @export
validate_mapping <- function(mapping, sep = ";") {
  required <- c("set_id", "set_name", "feature_ids")
  missing <- setdiff(required, colnames(mapping))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  mapping <- mapping |> dplyr::mutate(
    feature_ids = gsub("\\s+", "", feature_ids),
    feature_list = strsplit(feature_ids, sep)
  )
  return(mapping)
}
