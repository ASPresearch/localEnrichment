#' Quantitative Enrichment Analysis (QEA / MSEA)
#'
#' Perform a quantitative enrichment analysis on a ranked list of features
#' (e.g., metabolites) using a running-score statistic similar to GSEA.
#'
#' Given an \code{EnrichmentSet} object and a named numeric vector of scores
#' (typically log2 fold changes, t-statistics, or other ranking measures),
#' this function computes an enrichment score (ES) for each set and an
#' empirical p-value based on permutations.
#'
#' @param eset An \code{EnrichmentSet} object containing the mapping from set
#'   identifiers to feature IDs.
#' @param scores A named numeric vector of scores for all features. The
#'   \strong{names} must be feature identifiers (e.g., PubChem, HMDB) and
#'   must be comparable to the feature IDs stored in \code{eset}. The vector
#'   will be internally sorted in decreasing order.
#' @param nperm Integer. Number of permutations used to estimate empirical
#'   p-values (default \code{1000}). Set to \code{0} to skip permutation-based
#'   p-values (p-values will be \code{NA}).
#' @param p_adjust Character string specifying the multiple testing correction
#'   method to apply to the permutation p-values. Passed to
#'   \code{\link[stats]{p.adjust}} (default \code{"BH"}).
#' @param min_set_size Minimum number of overlapping features required for a
#'   set to be tested (default \code{3}).
#' @param max_set_size Maximum number of overlapping features allowed for a
#'   set to be tested (default \code{500}).
#' @param return_profiles Logical; if \code{TRUE}, the function also returns
#'   the running-score profile for each tested set (useful for plotting
#'   enrichment profiles). Default \code{FALSE}.
#' @param quiet Logical; if \code{FALSE}, progress information is printed
#'   (number of tested and skipped sets, etc.). Default \code{TRUE}.
#'
#' @details
#' For each set, the function:
#' \enumerate{
#'   \item Intersects the set features with the ranked universe (names of
#'         \code{scores}) and filters by \code{min_set_size} and
#'         \code{max_set_size}.
#'   \item Computes a running-score profile where hits increase the score
#'         by \code{1/Nh} and misses decrease it by \code{1/(N - Nh)},
#'         with \code{Nh} the number of hits and \code{N} the total number
#'         of ranked features.
#'   \item Defines the observed Enrichment Score (ES) as the maximum (or
#'         minimum, for negative enrichment) of the running-score profile.
#'   \item Estimates a permutation-based p-value by comparing the observed
#'         ES to the ES distribution obtained by random permutations of
#'         feature positions.
#' }
#'
#' The resulting table includes the set identifier, set name, set size,
#' ES, enrichment direction (positive/negative), raw p-value and adjusted
#' p-value, and is sorted by \code{p_value}.
#'
#' @return An object of class \code{"QEAResult"}, which is a list with:
#' \itemize{
#'   \item \code{table}: a \code{data.frame} with one row per tested set, with
#'         columns \code{set_id}, \code{set_name}, \code{n_in_set},
#'         \code{ES}, \code{direction}, \code{p_value}, \code{p_adj}.
#'   \item \code{profiles} (optional): if \code{return_profiles = TRUE}, a
#'         named list of \code{data.frame}s with the running-score profile
#'         (columns \code{set_id}, \code{position}, \code{running_score},
#'         \code{hit}, \code{set_name}) for each tested set.
#' }
#'
#' @examples
#' \dontrun{
#' # Toy example
#' meta <- list(
#'   mapping_name    = "toy_QEA",
#'   feature_id_type = "HMDB",
#'   feature_species = "Homo sapiens",
#'   set_source      = "ExampleDB",
#'   version         = "v1",
#'   description     = "Toy example for QEA"
#' )
#'
#' df <- data.frame(
#'   set_id     = c("A", "B"),
#'   set_name   = c("Pathway A", "Pathway B"),
#'   feature_ids = c("x1;x2;x3;x4", "x3;x4;x5;x6")
#' )
#'
#' eset <- EnrichmentSet(df, meta)
#' scores <- c(x1 = 2.0, x2 = 1.5, x3 = -1.0, x4 = 0.5,
#'             x5 = -2.5, x6 = -1.8)
#'
#' res <- qea_test(eset, scores, nperm = 100, quiet = FALSE)
#' head(res$table)
#' }
#'
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

  results  <- list()
  profiles <- list()
  tested   <- 0L
  skipped  <- 0L

  for (set_id in names(sets_list)) {
    set_feats <- intersect(sets_list[[set_id]], all_features)
    Nh <- length(set_feats)
    if (Nh == 0L) { skipped <- skipped + 1L; next }
    if (Nh < min_set_size || Nh > max_set_size) { skipped <- skipped + 1L; next }

    tested <- tested + 1L
    hits_idx <- match(set_feats, all_features)

    # Observed ES
    obs <- es_profile(hits_idx, N, Nh)
    ES  <- obs$ES
    direction <- ifelse(ES >= 0, "positive", "negative")

    # Permutations
    if (nperm > 0) {
      perm_ES <- numeric(nperm)
      for (i in seq_len(nperm)) {
        perm_idx <- sample.int(N, Nh)
        perm_obs <- es_profile(perm_idx, N, Nh)
        perm_ES[i] <- perm_obs$ES
      }

      if (direction == "positive") {
        p_value <- mean(perm_ES >= ES)
      } else {
        p_value <- mean(perm_ES <= ES)
      }
    } else {
      p_value <- NA_real_
    }

    results[[set_id]] <- data.frame(
      set_id    = set_id,
      n_in_set  = Nh,
      ES        = ES,
      direction = direction,
      p_value   = p_value
    )

    if (return_profiles) {
      profiles[[set_id]] <- data.frame(
        set_id        = set_id,
        position      = seq_len(N),
        running_score = obs$running_score,
        hit           = obs$hits
      )
    }
  }

  if (length(results) == 0)
    stop("No sets had overlapping features or met size constraints.")

  res_table <- dplyr::bind_rows(results) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = p_adjust))

  # Add human-readable set names
  set_info <- eset@data %>%
    dplyr::distinct(.data$set_id, .data$set_name)

  res_table <- res_table %>%
    dplyr::left_join(set_info, by = "set_id") %>%
    dplyr::relocate(set_id, set_name, .before = n_in_set) %>%
    dplyr::arrange(p_value)

  if (!quiet) {
    message("QEA tested sets: ", tested,
            " | skipped (no overlap/size): ", skipped,
            " | returned: ", nrow(res_table))
  }

  out <- list(table = res_table)

  if (return_profiles) {
    set_info_named <- set_info$set_name
    names(set_info_named) <- set_info$set_id
    profiles <- lapply(names(profiles), function(sid) {
      df <- profiles[[sid]]
      df$set_name <- set_info_named[[sid]]
      df
    })
    names(profiles) <- names(results)
    out$profiles <- profiles
  }

  class(out) <- "QEAResult"
  out
}
