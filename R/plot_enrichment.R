#' Plot enrichment results
#'
#' Visualize results from enrichment analyses such as ORA or QEA.
#'
#' @param results A tibble produced by \code{ora_test()} or \code{qea_test()}.
#' @param type Type of plot: \code{"bar"} (default) or \code{"dot"}.
#' @param top Integer: number of top sets to display (default = 20).
#' @param x_measure Variable for x-axis when type = "bar":
#'   one of \code{"-log10(p_adj)"} (default), \code{"-log10(p_value)"}, or
#'   \code{"fold_enrichment"}.
#' @param color_by Variable for coloring points when \code{type = "dot"}:
#'   one of \code{"p_adj"} (default), \code{"p_value"}, or \code{"fold_enrichment"}.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#' @export
plot_enrichment <- function(results,
                            type = c("bar", "dot"),
                            top = 20,
                            x_measure = "-log10(p_adj)",
                            color_by = "p_adj",
                            title = NULL) {

  #--- 1. Validation ---------------------------------------------------------
  if (!inherits(results, "data.frame")) {
    stop("Input 'results' must be a data frame or tibble (from ora_test).")
  }
  type <- match.arg(type)

  if (!("set_name" %in% names(results))) {
    stop("Results must contain a column 'set_name'.")
  }
  if (!("p_value" %in% names(results))) {
    stop("Results must contain 'p_value'.")
  }

  df <- results |>
    dplyr::mutate(
      neglogP = -log10(p_value),
      neglogPadj = if ("p_adj" %in% names(results)) -log10(p_adj) else NA,
      set_name = forcats::fct_reorder(set_name, -p_value)
    ) |>
    dplyr::slice_head(n = top)

  #--- 2. Select plot type ---------------------------------------------------
  if (type == "bar") {
    xvar <- switch(x_measure,
                   "-log10(p_adj)" = "neglogPadj",
                   "-log10(p_value)" = "neglogP",
                   "fold_enrichment" = "fold_enrichment",
                   stop("Invalid 'x_measure'."))

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[xvar]],
        y = .data$set_name
      )
    ) +
      ggplot2::geom_col(fill = "#4E79A7") +
      ggplot2::labs(
        x = x_measure,
        y = NULL,
        title = title %||% "Top enriched sets"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(face = "bold")
      )

  } else if (type == "dot") {
    if (!("fold_enrichment" %in% names(results))) {
      stop("Results must contain 'fold_enrichment' for dot plot.")
    }

    color_var <- match.arg(color_by, c("p_adj", "p_value", "fold_enrichment"))

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$fold_enrichment,
        y = .data$set_name,
        size = .data$n_in_set,
        color = .data[[color_var]]
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c(direction = -1, option = "C") +
      ggplot2::labs(
        x = "Fold enrichment",
        y = NULL,
        size = "Set size",
        color = color_var,
        title = title %||% "Top enriched sets"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(face = "bold")
      )
  }

  p
}

