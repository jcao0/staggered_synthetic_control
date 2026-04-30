#' Synthetic Control Weights for a Single Treated Unit
#'
#' Estimate synthetic control weights by solving a constrained quadratic
#' program on demeaned pre-treatment outcomes: minimise
#' \eqn{\|\tilde Y_1 - \tilde X b\|^2} subject to
#' \eqn{\sum b_j = 1, \; b_j \ge 0}, where \eqn{\tilde Y_1} and
#' \eqn{\tilde X} are time-demeaned series for the treated unit and controls.
#'
#' @param Y Numeric matrix (\eqn{N \times T}).
#'   The first row is the treated unit; remaining rows are donor (control)
#'   units.  Each column is a pre-treatment time period.
#'
#' @return A list with components
#' \describe{
#'   \item{a_hat}{Scalar intercept
#'     \eqn{\hat a = \bar Y_1 - \bar X' \hat b}.}
#'   \item{b_hat}{Numeric vector of length \eqn{N}.  Entry 1 is 0 (the
#'     treated unit's self-weight); entries \eqn{2, \dots, N} are the
#'     non-negative weights summing to 1.}
#' }
#'
#' @details
#' The QP is solved by \code{\link[quadprog]{solve.QP}}.  A small ridge term
#' (\eqn{10^{-6} I}) is added to the Hessian for numerical stability when
#' \eqn{T} is close to or smaller than \eqn{N-1}.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(5 * 20), 5, 20)  # 5 units, 20 pre-treatment periods
#' res <- synthetic_control(Y)
#' res$b_hat   # SC weights for unit 1
#'
#' @export
synthetic_control <- function(Y) {

  Y <- as.matrix(Y)
  N <- nrow(Y)
  Tp <- ncol(Y)          # number of pre-treatment periods
  n_ctrl <- N - 1L

  # ---- separate treated / controls ----
  y <- Y[1, ]                         # T-vector  (treated)
  X <- t(Y[-1, , drop = FALSE])       # T x (N-1) (controls, cols = units)

  # ---- demean each series ----
  y_dm <- y - mean(y)
  X_dm <- sweep(X, 2L, colMeans(X))   # subtract column means

  # ---- QP: min 1/2 b'Db - d'b  s.t. A'b >= b0 ----
  Dmat <- crossprod(X_dm)                         # (N-1) x (N-1)
  Dmat <- Dmat + Dmat                             # multiply by 2
  Dmat <- Dmat + 1e-6 * diag(n_ctrl)             # ridge for PD
  dvec <- 2 * drop(crossprod(X_dm, y_dm))         # (N-1)

  # constraints: first equality  sum(b) = 1, then b_j >= 0
  Amat <- cbind(rep(1, n_ctrl), diag(n_ctrl))
  bvec <- c(1, rep(0, n_ctrl))

  sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1L)
  b <- sol$solution

  # intercept
  a_hat <- mean(y) - drop(crossprod(colMeans(X), b))

  # full weight vector with 0 for the treated unit (position 1)
  b_hat <- c(0, b)

  list(a_hat = a_hat, b_hat = b_hat)
}


#' Synthetic Control Weights for All Units (Batch)
#'
#' For each unit in turn, treat that unit as the "treated" unit and
#' estimate SC weights from the remaining units.  This produces an
#' \eqn{N \times N} weight matrix \eqn{\hat B} with zeros on the diagonal.
#'
#' @param Y Numeric matrix (\eqn{N \times T}) of pre-treatment outcomes.
#'   Rows are units, columns are time periods.
#'
#' @return A list with components
#' \describe{
#'   \item{a_hat}{Numeric vector of length \eqn{N}: unit-level intercepts.}
#'   \item{B_hat}{Numeric \eqn{N \times N} matrix of SC weights.  Row
#'     \eqn{i} contains the weights used to construct the synthetic control
#'     for unit \eqn{i}; \eqn{B_{ii} = 0}.}
#' }
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(5 * 20), 5, 20)
#' res <- synthetic_control_batch(Y)
#' res$B_hat   # N x N weight matrix
#'
#' @export
synthetic_control_batch <- function(Y) {

  Y <- as.matrix(Y)
  N <- nrow(Y)

  a_hat <- numeric(N)
  B_hat <- matrix(0, N, N)

  for (i in seq_len(N)) {
    # put unit i first, others after
    Y_i <- rbind(Y[i, , drop = FALSE], Y[-i, , drop = FALSE])
    res  <- synthetic_control(Y_i)

    a_hat[i] <- res$a_hat

    # map (N-1) weights back to full N-vector (0 at position i)
    b_full      <- numeric(N)
    b_full[-i]  <- res$b_hat[-1]   # drop the leading 0
    B_hat[i, ]  <- b_full
  }

  list(a_hat = a_hat, B_hat = B_hat)
}
