#' @importFrom magrittr %>%
NULL

#' Prepare a features list for enrichment or pathway analysis
#'
#' @description
#' Filters and extracts a list of features for enrichment analyses
#' (ORA, MSEA, ChemRich, etc.) based on correlation and adjusted p-value thresholds.
#'
#' @param features_table A data frame containing feature-level results for one dimension.
#'        Must include at least `correlation` and `p.adj`.
#' @param fdr_cutoff Numeric. Maximum adjusted p-value allowed (default = 1, i.e., no filter).
#' @param cor_cutoff Numeric. Minimum absolute correlation required (default = 0).
#' @param fields Character vector of column names to return (default = "Feature_ID").
#' @param save_file Logical. Whether to save the filtered list as a TSV file (default = FALSE).
#' @param prefix Character. Prefix for the output file name (default = inferred from `Dimension` column if available).
#' @param outdir Character. Output directory for the saved file (default = "results/features_lists").
#' @param colnames Logical. Whether to include column names in the saved TSV file (default = FALSE).
#'
#' @return A tibble containing the filtered features with the requested fields.
#' @export
#'
#' @examples
#' # Example feature table for one dimension
#' features_example <- data.frame(
#'   Feature_ID = paste0("HMDB0000", 161:170),
#'   correlation = c(0.35, 0.42, 0.12, 0.29, -0.25, -0.32, -0.18, 0.10, 0.05, -0.02),
#'   p.adj = c(0.01, 0.04, 0.30, 0.10, 0.25, 0.60, 0.80, 0.05, 0.90, 0.40),
#'   chemical_families = rep(c("Amino Acids", "Biogenic Amines"), each = 5),
#'   stringsAsFactors = FALSE
#' )
#'
#' # 1. Filtered list for ORA-style enrichment
#' ora_list <- prepare_features_list(
#'   features_table = features_example,
#'   fdr_cutoff = 0.25,
#'   cor_cutoff = 0.1,
#'   fields = c("Feature_ID")
#' )
#' ora_list
#'
#' # 2. Full list for MSEA-style enrichment
#' msea_list <- prepare_features_list(
#'   features_table = features_example,
#'   fdr_cutoff = 1,
#'   cor_cutoff = 0,
#'   fields = c("Feature_ID", "correlation", "chemical_families")
#' )
#' msea_list
prepare_features_list <- function(features_table,
                                  fdr_cutoff = 1,
                                  cor_cutoff = 0,
                                  fields = c("Feature_ID"),
                                  save_file = FALSE,
                                  prefix = NULL,
                                  outdir = "results/features_lists",
                                  colnames = FALSE) {

  # --- 1. Check that required columns exist ---------------------------------
  required <- c("correlation", "p.adj")
  missing_required <- setdiff(required, names(features_table))
  if (length(missing_required) > 0) {
    stop("Input table must contain columns: ", paste(required, collapse = ", "))
  }

  missing_fields <- setdiff(fields, names(features_table))
  if (length(missing_fields) > 0) {
    stop("Missing requested columns: ", paste(missing_fields, collapse = ", "))
  }

  # --- 2. Filter rows -------------------------------------------------------
  filtered <- features_table %>%
    dplyr::filter(p.adj <= fdr_cutoff, abs(correlation) >= cor_cutoff) %>%
    dplyr::arrange(dplyr::desc(abs(correlation)))

  if (nrow(filtered) == 0) {
    warning("No features passed the specified thresholds (fdr_cutoff = ",
            fdr_cutoff, ", cor_cutoff = ", cor_cutoff, ").")
  }

  # --- 3. Select desired fields --------------------------------------------
  result <- dplyr::select(filtered, dplyr::all_of(fields))

  # --- 4. Optionally save --------------------------------------------------
  if (save_file) {
    if (is.null(prefix)) {
      prefix <- if ("Dimension" %in% names(features_table)) {
        unique(features_table$Dimension)[1]
      } else {
        "features"
      }
    }
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    outname <- file.path(outdir, paste0(prefix, "_filtered.tsv"))
    readr::write_tsv(result, outname, col_names = colnames)
    message("Saved filtered list to: ", outname)
  }

  # --- 5. Return -----------------------------------------------------------
  tibble::as_tibble(result)
}

