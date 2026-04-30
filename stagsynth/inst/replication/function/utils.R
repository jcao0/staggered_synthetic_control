#' Convert Long-Format Panel Data to Matrices
#'
#' Transform a data frame in long format (one row per unit-period) into
#' the \eqn{N \times T} matrices \code{Y} and \code{D} expected by
#' \code{\link{ssc}}.
#'
#' @param data A data frame.
#' @param unit Character: name of the unit identifier column.
#' @param time Character: name of the time period column.
#' @param outcome Character: name of the outcome variable column.
#' @param treatment Character: name of the treatment indicator column
#'   (must be 0/1).
#'
#' @return A list with components
#' \describe{
#'   \item{Y}{Numeric \eqn{N \times T} outcome matrix.}
#'   \item{D}{Numeric \eqn{N \times T} treatment matrix.}
#'   \item{units}{Sorted vector of unique unit identifiers.}
#'   \item{times}{Sorted vector of unique time periods.}
#' }
#'
#' @examples
#' df <- data.frame(
#'   id   = rep(1:4, each = 6),
#'   time = rep(1:6, times = 4),
#'   Y    = rnorm(24),
#'   D    = c(rep(0, 12), rep(c(0,0,0,1,1,1), 2))
#' )
#' mat <- panel_to_matrices(df, unit = "id", time = "time",
#'                          outcome = "Y", treatment = "D")
#'
#' @export
panel_to_matrices <- function(data, unit, time, outcome, treatment) {

  stopifnot(
    is.data.frame(data),
    all(c(unit, time, outcome, treatment) %in% names(data))
  )

  data[[outcome]] <- as.numeric(data[[outcome]])
  data <- data[order(data[[unit]], data[[time]]), ]

  units <- sort(unique(data[[unit]]))
  times <- sort(unique(data[[time]]))
  N <- length(units)
  Tt <- length(times)

  Y <- matrix(NA_real_, N, Tt)
  D <- matrix(NA_real_, N, Tt)
  rownames(Y) <- rownames(D) <- as.character(units)
  colnames(Y) <- colnames(D) <- as.character(times)

  for (i in seq_len(N)) {
    idx <- data[[unit]] == units[i]
    sub <- data[idx, ]
    # match times to column positions
    col_pos <- match(sub[[time]], times)
    Y[i, col_pos] <- sub[[outcome]]
    D[i, col_pos] <- sub[[treatment]]
  }

  # check for missing values
  if (anyNA(Y))
    warning("Outcome matrix Y contains NAs (unbalanced panel?).",
            call. = FALSE)
  if (anyNA(D))
    warning("Treatment matrix D contains NAs (unbalanced panel?).",
            call. = FALSE)

  list(Y = Y, D = D, units = units, times = times)
}


#' Compute Smallest Eigenvalue of the SSC Design Matrix
#'
#' A diagnostic function that builds the SSC design matrix
#' \eqn{\sum_s A_s' \hat M A_s} and returns its smallest eigenvalue.
#' This matrix must be positive definite for SSC estimates to exist.
#'
#' @param Y Numeric matrix (\eqn{N \times T_{total}}) of outcomes.
#' @param D Binary matrix (\eqn{N \times T_{total}}) of treatment indicators.
#' @param S Number of post-treatment periods (or \code{NULL} for all).
#'
#' @return A scalar: the smallest eigenvalue.
#'
#' @details
#' A positive value means the SSC estimator is well-defined; a value
#' near zero warns that identification is weak.
#'
#' @examples
#' set.seed(1)
#' N <- 10; Ttot <- 8
#' Y <- matrix(rnorm(N * Ttot), N, Ttot)
#' D <- matrix(0L, N, Ttot)
#' D[1:3, 5:Ttot] <- 1L   # units 1-3 treated from period 5
#' ssc_min_eigenvalue(Y, D, S = 2)
#'
#' @export
ssc_min_eigenvalue <- function(Y, D, S = NULL) {

  res <- ssc(Y, D, S = S, alpha = 0.05)
  res$min_eigenvalue
}
