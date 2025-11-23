#' Plot QEA Running Score Profiles
#'
#' @description Plot the enrichment running-score curves produced by
#'   \code{qea_test(return_profiles = TRUE)}.
#'
#' @param x A \code{QEAResult} object returned by \code{qea_test()}.
#' @param sets Optional vector of set IDs or set names to plot.
#' @param show One of "all", "selected", or "significant".
#' @param p_cutoff FDR threshold when show = "significant".
#' @param facet Logical; whether to facet plots by set.
#' @param highlight Color for the profile line.
#' @param rug_color Color for tick marks indicating hits.
#'
#' @import ggplot2
#' @return A ggplot object.
#' @export

plot_qea_profile <- function(x,
                             sets = NULL,
                             show = c("all", "selected", "significant"),
                             p_cutoff = 0.05,
                             facet = FALSE,
                             highlight = "#4E79A7",
                             rug_color = "black") {

  show <- match.arg(show)

  if (!inherits(x, "QEAResult"))
    stop("x must be a QEAResult produced by qea_test(return_profiles=TRUE).")

  profiles <- x$profiles
  tbl <- x$table

  # Convertir noms humans a set_id
  if (!is.null(sets)) {
    lookup <- tbl[, c("set_id", "set_name")]
    sets <- sapply(sets, function(s) {
      if (s %in% lookup$set_id) return(s)
      if (s %in% lookup$set_name) return(lookup$set_id[match(s, lookup$set_name)])
      stop(paste("Set", s, "not found."))
    })
  }

  # Filtrar segons show
  if (show == "significant")
    sets <- tbl$set_id[tbl$p_adj <= p_cutoff]
  if (show == "selected" && is.null(sets))
    stop("You must supply 'sets' when show='selected'.")

  df <- dplyr::bind_rows(profiles[sets], .id = "set_id")

  df$set_name <- tbl$set_name[match(df$set_id, tbl$set_id)]

  gg <- ggplot(df, aes(position, running_score, color = set_name)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_rug(data = subset(df, hit == 1),
             aes(position, running_score),
             inherit.aes = FALSE,
             sides = "b",
             color = rug_color) +
    labs(x = "Rank", y = "Running ES", color = "Set",
         title = "QEA Running Score Profiles") +
    theme_minimal(base_size = 12)

  if (facet)
    gg <- gg + facet_wrap(~ set_name, scales = "free_y")

  gg
}
