#' Plot QEA enrichment profiles
#'
#' @param qea_res A \code{QEAResult} from \code{qea_test(return_profiles=TRUE)}.
#' @param sets Character vector of set names to display. If NULL, behavior depends on `show`.
#' @param show One of "significant", "selected", or "all".
#' @param p_cutoff Numeric; significance threshold for `show="significant"` (default 0.05).
#' @param facet Logical; facet by set (TRUE) or overlay curves (FALSE).
#' @param highlight Color for enrichment curves (used if overlaying).
#' @param rug_color Color for tick marks under each panel / curve.
#'
#' @return A ggplot object.
#' @export
plot_qea_profile <- function(qea_res,
                             sets = NULL,
                             show = c("significant", "selected", "all"),
                             p_cutoff = 0.05,
                             facet = TRUE,
                             highlight = "#4E79A7",
                             rug_color = "black") {

  show <- match.arg(show)

  if (!inherits(qea_res, "QEAResult"))
    stop("'qea_res' must be a QEAResult from qea_test().")
  if (is.null(qea_res$profiles))
    stop("This QEAResult has no profiles. Re-run qea_test(return_profiles=TRUE).")

  tbl <- qea_res$table
  prof <- dplyr::bind_rows(qea_res$profiles)

  chosen_sets <- sets
  if (is.null(chosen_sets)) {
    if (show == "significant") {
      chosen_sets <- tbl$set_name[tbl$p_adj < p_cutoff]
    } else if (show == "all") {
      chosen_sets <- tbl$set_name
    } else { # "selected" sense sets -> error clar
      stop("When show='selected', please provide 'sets'.")
    }
  }
  chosen_sets <- unique(chosen_sets)

  plot_df <- dplyr::filter(prof, .data$set_name %in% chosen_sets)
  if (nrow(plot_df) == 0)
    stop("No profiles to plot for the selected sets.")

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$position, y = .data$running_score,
                 color = .data$set_name, group = .data$set_name)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_line(linewidth = 1)

  if (facet) {
    p <- p +
      ggplot2::geom_rug(data = subset(plot_df, hit == 1),
                        sides = "b", color = rug_color, linewidth = 0.3) +
      ggplot2::facet_wrap(~ set_name, scales = "free_y") +
      ggplot2::guides(color = "none")
  } else {
    p <- p +
      ggplot2::geom_rug(data = subset(plot_df, hit == 1),
                        sides = "b", color = rug_color, linewidth = 0.3) +
      ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(title = "Set"))
  }

  p +
    ggplot2::labs(
      x = "Ranked features",
      y = "Running enrichment score",
      title = "QEA enrichment profiles"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
}
