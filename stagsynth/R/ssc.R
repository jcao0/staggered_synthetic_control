#' Staggered Synthetic Control Estimation
#'
#' Estimate treatment effects in a panel with staggered adoption using
#' the Staggered Synthetic Control (SSC) method of Cao, Lu, and Wu (2020).
#' Returns event-time ATT, overall ATT, heterogeneous treatment effects,
#' and placebo-in-time confidence intervals.
#'
#' @param Y Numeric matrix (\eqn{N \times T_{total}}) of outcomes.
#'   Rows are units, columns are time periods (pre- and post-treatment).
#' @param D Binary matrix (\eqn{N \times T_{total}}) of treatment indicators.
#'   \code{D[i,t] = 1} if unit \eqn{i} is treated at time \eqn{t}.
#'   Treatment must be absorbing (once treated, always treated).
#' @param S Integer or \code{NULL}.  Number of post-treatment periods to use.
#'   If \code{NULL} (default), all available post-treatment periods are used.
#' @param alpha Significance level for confidence intervals (default 0.05).
#'
#' @return An object of class \code{"ssc"}, a list containing:
#' \describe{
#'   \item{att_event}{Numeric vector of length \eqn{S}: event-time ATT
#'     estimates (averaged across units at each event time).}
#'   \item{ci_lower_event, ci_upper_event}{Numeric vectors of length
#'     \eqn{S}: lower and upper bounds of
#'     \eqn{(1 - \alpha)} placebo-in-time confidence intervals.
#'     \code{NA} when \eqn{T < S} (too few pre-treatment periods).}
#'   \item{att_overall}{Scalar: overall ATT (simple average of all
#'     heterogeneous effects).}
#'   \item{ci_lower_overall, ci_upper_overall}{Scalar CI bounds for
#'     the overall ATT.  \code{NA} when \eqn{T < S}.}
#'   \item{p_value}{Two-sided p-value for \eqn{H_0: ATT = 0} based on
#'     the placebo distribution.  \code{NA} when \eqn{T < S}.}
#'   \item{gamma_hat}{Numeric vector of length \eqn{K}: heterogeneous
#'     treatment effects for every treated (unit, post-period) pair.}
#'   \item{te_mat_hat}{Numeric \eqn{N \times S} matrix of
#'     unit-level treatment effects at each post-treatment period.}
#'   \item{B_hat}{Numeric \eqn{N \times N} SC weight matrix.}
#'   \item{a_hat}{Numeric vector of length \eqn{N}: SC intercepts.}
#'   \item{u_hat}{Numeric \eqn{N \times T} matrix of pre-treatment
#'     SC residuals.}
#'   \item{min_eigenvalue}{Smallest eigenvalue of the sample analogue of
#'     the design matrix \eqn{\sum_s A_s' \hat M A_s}.
#'     Must be positive for the estimator to be well-defined.}
#'   \item{index_mat}{Integer \eqn{K \times 3} matrix.  Each row
#'     \eqn{(s, i, e)} records the post-treatment period \eqn{s},
#'     unit \eqn{i}, and event time \eqn{e} for one element of
#'     \eqn{\hat\gamma}.}
#'   \item{N, T, S, K}{Panel dimensions.}
#'   \item{alpha}{Significance level used.}
#' }
#'
#' @details
#' The SSC method proceeds in four steps:
#' \enumerate{
#'   \item \strong{SC weights.}
#'     For every unit, estimate synthetic control weights from
#'     pre-treatment data.
#'   \item \strong{Treatment structure.}
#'     Build the treatment assignment matrices \eqn{A_s} that map
#'     heterogeneous effects \eqn{\gamma} to unit-level outcomes at
#'     each post-treatment period.
#'   \item \strong{Estimation.}
#'     Solve a GLS-type system to recover \eqn{\hat\gamma}, then
#'     aggregate to event-time or overall ATT via a linear map \eqn{L}.
#'   \item \strong{Inference.}
#'     Construct a null distribution by applying the same estimator to
#'     rolling windows of pre-treatment residuals (placebo-in-time).
#'     Confidence intervals are the \eqn{\alpha/2} and
#'     \eqn{1 - \alpha/2} quantiles of this distribution, shifted by
#'     the point estimate.
#' }
#'
#' @references
#' Cao, J., Lu, C., and Wu, Y. (2020).
#' "Synthetic Control Inference for Staggered Adoption."
#'
#' @examples
#' set.seed(1)
#' N <- 5; Ttot <- 15
#' Y <- matrix(rnorm(N * Ttot), N, Ttot)
#' D <- matrix(0L, N, Ttot)
#' D[1, 8:15] <- 1L
#' D[2, 10:15] <- 1L
#' result <- ssc(Y, D)
#' print(result)
#' summary(result)
#'
#' @export
ssc <- function(Y, D, S = NULL, alpha = 0.05) {

  cl <- match.call()
  Y  <- as.matrix(Y)
  D  <- as.matrix(D)

  # ---- validate ----
  stopifnot(
    "Y and D must have the same dimensions" = identical(dim(Y), dim(D)),
    "D must be binary (0/1)" = all(D %in% c(0, 1)),
    "alpha must be in (0, 1)" = alpha > 0 && alpha < 1
  )

  N      <- nrow(Y)
  T_total <- ncol(Y)

  # ---- check absorbing treatment ----
  for (i in seq_len(N)) {
    d_i <- D[i, ]
    if (any(d_i == 1)) {
      first_treat <- which(d_i == 1)[1]
      if (!all(d_i[first_treat:T_total] == 1))
        stop(sprintf(
          "Unit %d: treatment is not absorbing (once treated, must stay treated).",
          i))
    }
  }

  # ---- determine T (last all-zero column before any treatment) ----
  first_treated_col <- which(colSums(D) > 0)[1]
  if (is.na(first_treated_col))
    stop("No treatment found in D.")
  Tp <- first_treated_col - 1L       # number of pre-treatment periods
  if (Tp < 2L)
    stop("Need at least 2 pre-treatment periods.")

  S_max <- T_total - Tp
  if (is.null(S))
    S <- S_max
  stopifnot("S exceeds available post-treatment periods" = S <= S_max)

  # ---- partition ----
  Y_use <- Y[, seq_len(Tp + S), drop = FALSE]
  D_use <- D[, seq_len(Tp + S), drop = FALSE]

  Y_T <- Y_use[, seq_len(Tp),             drop = FALSE]   # N x T
  Y_S <- Y_use[, Tp + seq_len(S),         drop = FALSE]   # N x S
  D_S <- D_use[, Tp + seq_len(S),         drop = FALSE]   # N x S

  # ---- build index matrix (post_time, unit, event_time) ----
  K <- sum(D_use)
  if (K == 0L)
    stop("No treated (unit, period) observations.")

  index_mat <- matrix(0L, K, 3L)
  colnames(index_mat) <- c("post_time", "unit", "event_time")
  ind <- 0L
  for (s in seq_len(S)) {
    for (i in seq_len(N)) {
      if (D_S[i, s] == 1) {
        ind <- ind + 1L
        index_mat[ind, ] <- c(s, i, sum(D_S[i, seq_len(s)]))
      }
    }
  }

  # ---- build A_s treatment structure matrices ----
  A <- array(0, dim = c(N, K, S))
  for (k in seq_len(K)) {
    A[ index_mat[k, 2], k, index_mat[k, 1] ] <- 1
  }

  # ---- SC weights from pre-treatment data ----
  sc <- synthetic_control_batch(Y_T)
  a_hat <- sc$a_hat
  B_hat <- sc$B_hat
  I_N   <- diag(N)
  IB    <- I_N - B_hat                              # N x N
  M_hat <- crossprod(IB)                             # (I-B)'(I-B)

  # ---- design matrix:  sum_s A_s' M A_s ----
  temp1 <- matrix(0, K, K)
  for (s in seq_len(S)) {
    As <- A[, , s, drop = FALSE]
    dim(As) <- c(N, K)
    temp1 <- temp1 + crossprod(As, M_hat %*% As)    # K x K
  }

  min_eig <- min(eigen(temp1, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 1e-12)
    warning(sprintf(
      "Design matrix is near-singular (min eigenvalue = %.2e). ",
      min_eig),
      "SSC estimates may be unreliable.", call. = FALSE)

  # ---- right-hand side:  sum_s A_s' (I-B)' [(I-B) Y_s - a] ----
  temp2 <- rep(0, K)
  for (s in seq_len(S)) {
    As <- A[, , s, drop = FALSE]
    dim(As) <- c(N, K)
    resid_s <- IB %*% Y_S[, s] - a_hat              # N-vector
    temp2   <- temp2 + drop(crossprod(As, t(IB)) %*% resid_s)
  }

  # ---- solve for gamma_hat ----
  gamma_hat <- drop(solve(temp1, temp2))

  # ---- event-time ATT ----
  L_event <- matrix(0, S, K)
  for (s in seq_len(S)) {
    mask <- index_mat[, 3] == s
    if (any(mask))
      L_event[s, mask] <- 1 / sum(mask)
  }
  att_event <- drop(L_event %*% gamma_hat)

  # ---- overall ATT ----
  att_overall <- mean(gamma_hat)

  # ---- pre-treatment residuals ----
  u_hat <- Y_T - (outer(a_hat, rep(1, Tp)) + B_hat %*% Y_T)   # N x T

  # ---- placebo-in-time inference ----
  n_placebo <- Tp - S
  has_inference <- n_placebo > 0L

  if (has_inference) {
    V_mat <- matrix(0, K, n_placebo)
    for (t in seq_len(n_placebo)) {
      v_temp <- rep(0, K)
      for (s in seq_len(S)) {
        As <- A[, , s, drop = FALSE]
        dim(As) <- c(N, K)
        v_temp <- v_temp + drop(crossprod(As, t(IB)) %*% u_hat[, t + s])
      }
      V_mat[, t] <- drop(solve(temp1, v_temp))
    }

    # event-time CI
    null_event  <- L_event %*% V_mat          # S x n_placebo
    ci_lower_ev <- numeric(S)
    ci_upper_ev <- numeric(S)
    for (s in seq_len(S)) {
      q_lo <- stats::quantile(null_event[s, ], probs = alpha / 2)
      q_hi <- stats::quantile(null_event[s, ], probs = 1 - alpha / 2)
      ci_lower_ev[s] <- att_event[s] - q_hi
      ci_upper_ev[s] <- att_event[s] - q_lo
    }

    # overall CI
    L_overall   <- matrix(1 / K, 1, K)
    null_overall <- drop(L_overall %*% V_mat)   # n_placebo vector
    q_lo <- stats::quantile(null_overall, probs = alpha / 2)
    q_hi <- stats::quantile(null_overall, probs = 1 - alpha / 2)
    ci_lower_ov <- att_overall - q_hi
    ci_upper_ov <- att_overall - q_lo

    # two-sided p-value
    p_value <- mean(abs(null_overall) > abs(att_overall))

  } else {
    ci_lower_ev <- rep(NA_real_, S)
    ci_upper_ev <- rep(NA_real_, S)
    ci_lower_ov <- NA_real_
    ci_upper_ov <- NA_real_
    p_value     <- NA_real_
    null_event  <- NULL
    null_overall <- NULL
    message(sprintf(
      "T = %d <= S = %d: not enough pre-treatment periods for ",
      Tp, S),
      "placebo-in-time inference.  Point estimates are available ",
      "but confidence intervals are NA.")
  }

  # ---- unit-level treatment effects ----
  te_mat_hat <- matrix(0, N, S)
  for (s in seq_len(S)) {
    As <- A[, , s, drop = FALSE]
    dim(As) <- c(N, K)
    te_mat_hat[, s] <- drop(As %*% gamma_hat)
  }

  # ---- assemble output ----
  out <- list(
    att_event        = att_event,
    ci_lower_event   = ci_lower_ev,
    ci_upper_event   = ci_upper_ev,
    att_overall      = att_overall,
    ci_lower_overall = ci_lower_ov,
    ci_upper_overall = ci_upper_ov,
    p_value          = p_value,
    gamma_hat        = gamma_hat,
    te_mat_hat       = te_mat_hat,
    B_hat            = B_hat,
    a_hat            = a_hat,
    u_hat            = u_hat,
    min_eigenvalue   = min_eig,
    index_mat        = index_mat,
    N = N, T = Tp, S = S, K = K,
    alpha = alpha,
    call = cl
  )
  class(out) <- "ssc"
  out
}
