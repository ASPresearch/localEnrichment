#' Quantitative Enrichment Analysis (QEA / GSEA-like)
#'
#' Performs a ranked enrichment analysis for all sets in an EnrichmentSet.
#'
#' @param eset An \code{EnrichmentSet} object.
#' @param scores Named numeric vector of feature scores (names = feature IDs).
#' @param nperm Number of permutations for p-value estimation (default = 1000).
#' @param p_adjust Method for multiple testing correction (default = "BH").
#' @param min_set_size Minimum number of features per set to test (default = 3).
#' @param max_set_size Maximum number of features per set (default = 500).
#' @param return_profiles Logical; if TRUE, returns running-score profiles per set.
#' @param quiet Logical; suppress progress messages.
#'
#' @return A list of class \code{"QEAResult"} with:
#'   \itemize{
#'     \item table: tibble with set_name, n_in_set, ES, p_value, p_adj, direction
#'     \item profiles (optional): list of data.frames (set_name, position, running_score, hit)
#'   }
#' @export
qea_test <- function(eset, scores,
                     nperm = 1000,
                     p_adjust = "BH",
                     min_set_size = 3,
                     max_set_size = 500,
                     return_profiles = FALSE,
                     quiet = TRUE) {

  if (!inherits(eset, "EnrichmentSet"))
    stop("'eset' must be an EnrichmentSet object.")
  if (is.null(names(scores)))
    stop("'scores' must be a named numeric vector with feature IDs as names.")

  # ranked universe
  scores_sorted <- sort(scores, decreasing = TRUE)
  all_features <- names(scores_sorted)
  N <- length(all_features)

  sets_list <- as(eset, "list")

  # helper to compute running-score + ES
  es_profile <- function(hits_idx, N, Nh) {
    hits <- integer(N); hits[hits_idx] <- 1L
    rs <- cumsum(hits / Nh - (!hits) / (N - Nh))
    ES <- if (max(rs) >= abs(min(rs))) max(rs) else min(rs)
    list(ES = ES, running_score = rs, hits = hits)
  }

  results <- list()
  profiles <- list()
  tested <- 0L
  skipped <- 0L

  for (set_name in names(sets_list)) {
    set_feats <- intersect(sets_list[[set_name]], all_features)
    Nh <- length(set_feats)
    if (Nh == 0L) { skipped <- skipped + 1L; next }
    if (Nh < min_set_size || Nh > max_set_size) { skipped <- skipped + 1L; next }

    tested <- tested + 1L
    hits_idx <- match(set_feats, all_features)

    # Observed ES
    obs <- es_profile(hits_idx, N, Nh)
    ES <- obs$ES
    direction <- ifelse(ES >= 0, "positive", "negative")

    # Permutations for empirical p-value
    if (nperm > 0) {
      perm_ES <- replicate(nperm, {
        # permutar les posicions dels hits dins l’univers ranquejat
        idx_perm <- sample.int(N, Nh, replace = FALSE)
        es_profile(idx_perm, N, Nh)$ES
      })
      p_value <- if (ES >= 0) mean(perm_ES >= ES) else mean(perm_ES <= ES)
      # assegura p > 0
      p_value <- max(p_value, 1 / (nperm + 1))
    } else {
      p_value <- NA_real_
    }

    results[[set_name]] <- data.frame(
      set_name = set_name,
      n_in_set = Nh,
      ES = ES,
      direction = direction,
      p_value = p_value
    )

    if (return_profiles) {
      profiles[[set_name]] <- data.frame(
        set_name = set_name,
        position = seq_len(N),
        running_score = obs$running_score,
        hit = obs$hits
      )
    }
  }

  if (length(results) == 0)
    stop("No sets had overlapping features or met size constraints.")

  res_table <- dplyr::bind_rows(results) |>
    dplyr::mutate(p_adj = p.adjust(p_value, method = p_adjust)) |>
    dplyr::arrange(p_value)

  if (!quiet) {
    message("QEA tested sets: ", tested,
            " | skipped (no overlap/size): ", skipped,
            " | returned: ", nrow(res_table))
  }

  out <- list(table = res_table)
  if (return_profiles) out$profiles <- profiles
  class(out) <- "QEAResult"
  out
}
