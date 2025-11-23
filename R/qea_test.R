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

  for (set_id in names(sets_list)) {
    set_feats <- intersect(sets_list[[set_id]], all_features)
    Nh <- length(set_feats)
    if (Nh == 0L) { skipped <- skipped + 1L; next }
    if (Nh < min_set_size || Nh > max_set_size) { skipped <- skipped + 1L; next }

    tested <- tested + 1L
    hits_idx <- match(set_feats, all_features)

    # Observed ES
    obs <- es_profile(hits_idx, N, Nh)
    ES <- obs$ES
    direction <- ifelse(ES >= 0, "positive", "negative")

    # ---- Permutacions (codi original restaurat) --------------------------
    if (nperm > 0) {
      perm_ES <- numeric(nperm)
      for (i in seq_len(nperm)) {
        perm_idx <- sample.int(N, Nh)
        perm_obs <- es_profile(perm_idx, N, Nh)
        perm_ES[i] <- perm_obs$ES
      }

      # p-values segons direcció
      if (direction == "positive") {
        p_value <- mean(perm_ES >= ES)
      } else {
        p_value <- mean(perm_ES <= ES)
      }

    } else {
      p_value <- NA_real_
    }
    # ----------------------------------------------------------------------

    results[[set_id]] <- data.frame(
      set_id = set_id,
      n_in_set = Nh,
      ES = ES,
      direction = direction,
      p_value = p_value
    )

    if (return_profiles) {
      profiles[[set_id]] <- data.frame(
        set_id = set_id,
        position = seq_len(N),
        running_score = obs$running_score,
        hit = obs$hits
      )
    }
  }

  if (length(results) == 0)
    stop("No sets had overlapping features or met size constraints.")

  res_table <- dplyr::bind_rows(results) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = p_adjust))

  # Afegim set_name humà
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
    # també afegim set_name als perfils
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
