# zzz.R — global variables and general imports
# ---------------------------------------------
# This file declares global variables and imports to avoid check() notes.

#' @keywords internal
#' @importFrom methods new
#' @importFrom stats fisher.test median p.adjust phyper setNames
#' @importFrom readr write_tsv
#' @importFrom magrittr %>%
NULL

# Declare global variables used in dplyr pipelines to silence R CMD check notes
utils::globalVariables(c(
  # common columns in enrichment outputs
  "p_value", "p_adj", "n_selected_in_set", "n_in_set", "fold_enrichment",
  "set_name", "feature_ids", "hit",
  # typical columns in features tables
  "correlation", "p.adj", "chemical_families",
  # internal variables used in dplyr
  ".data"
))
