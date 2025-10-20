# utilitats internes

#' Null-coalescing operator
#'
#' Returns the left value if not NULL, otherwise the right value.
#'
#' @param a First object.
#' @param b Second object.
#' @return a if it's not NULL, otherwise b.
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
