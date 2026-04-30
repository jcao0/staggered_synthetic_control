#' stagsynth: Staggered Synthetic Control Estimation and Inference
#'
#' Implements the Staggered Synthetic Control (SSC) method of
#' Cao, Lu, and Wu (2020) for estimating treatment effects in
#' panel data with staggered adoption.
#'
#' @section Main function:
#' \code{\link{ssc}}: Estimate event-time ATT, overall ATT, and
#' placebo-in-time confidence intervals.
#'
#' @section Utilities:
#' \itemize{
#'   \item \code{\link{panel_to_matrices}}: Convert long-format panel
#'     data to the \eqn{N \times T} matrices expected by \code{ssc()}.
#'   \item \code{\link{ssc_min_eigenvalue}}: Check the design matrix
#'     invertibility condition.
#'   \item \code{\link{synthetic_control}}: Estimate SC weights for a
#'     single treated unit.
#'   \item \code{\link{synthetic_control_batch}}: Estimate SC weights
#'     for all units.
#' }
#'
#' @docType package
#' @name stagsynth-package
"_PACKAGE"
