## ---- synthetic_control tests ----

test_that("synthetic_control returns correct structure", {
  set.seed(42)
  N <- 5; Tp <- 20
  Y <- matrix(rnorm(N * Tp), N, Tp)

  res <- synthetic_control(Y)

  expect_type(res, "list")
  expect_named(res, c("a_hat", "b_hat"))
  expect_length(res$b_hat, N)
  expect_equal(res$b_hat[1], 0)              # self-weight is 0
  expect_true(all(res$b_hat >= -1e-8))       # non-negative
  expect_equal(sum(res$b_hat), 1, tolerance = 1e-6)  # sum to 1
})


test_that("synthetic_control_batch returns N x N weight matrix", {
  set.seed(42)
  N <- 5; Tp <- 20
  Y <- matrix(rnorm(N * Tp), N, Tp)

  res <- synthetic_control_batch(Y)

  expect_type(res, "list")
  expect_named(res, c("a_hat", "B_hat"))
  expect_length(res$a_hat, N)
  expect_equal(dim(res$B_hat), c(N, N))

  # diagonal should be zero
  expect_equal(diag(res$B_hat), rep(0, N))

  # each row sums to 1
  row_sums <- rowSums(res$B_hat)
  expect_equal(row_sums, rep(1, N), tolerance = 1e-6)

  # all weights non-negative
  expect_true(all(res$B_hat >= -1e-8))
})


## ---- panel_to_matrices tests ----

test_that("panel_to_matrices correctly reshapes long data", {
  df <- data.frame(
    id   = rep(1:3, each = 5),
    time = rep(1:5, 3),
    Y    = rnorm(15),
    D    = c(rep(0, 5), rep(c(0,0,0,1,1), 2))
  )

  mat <- panel_to_matrices(df, unit = "id", time = "time",
                           outcome = "Y", treatment = "D")

  expect_equal(dim(mat$Y), c(3, 5))
  expect_equal(dim(mat$D), c(3, 5))
  expect_equal(mat$units, 1:3)
  expect_equal(mat$times, 1:5)

  # unit 1 is never treated
  expect_equal(unname(mat$D[1, ]), rep(0, 5))

  # units 2-3 treated from period 4
  expect_equal(unname(mat$D[2, ]), c(0, 0, 0, 1, 1))
  expect_equal(unname(mat$D[3, ]), c(0, 0, 0, 1, 1))
})


## ---- ssc main function tests ----

test_that("ssc runs on simple simulated data", {
  set.seed(123)
  N <- 10; Tp <- 30; S <- 5

  # simulate factor model
  f <- cumsum(rnorm(Tp + S))
  lambda <- runif(N, 0.5, 1.5)
  Y0 <- outer(lambda, f) + matrix(rnorm(N * (Tp + S), sd = 0.5), N)

  # treatment: units 1-4 treated at staggered times
  D <- matrix(0, N, Tp + S)
  D[1, (Tp + 1):(Tp + S)] <- 1          # unit 1: treated from T+1
  D[2, (Tp + 3):(Tp + S)] <- 1          # unit 2: treated from T+3
  D[3, (Tp + 4):(Tp + S)] <- 1          # unit 3: treated from T+4
  D[4, (Tp + 5):(Tp + S)] <- 1          # unit 4: treated from T+5

  # add treatment effect = 2
  Y <- Y0 + 2 * D

  result <- ssc(Y, D, alpha = 0.05)

  # class and structure
  expect_s3_class(result, "ssc")
  expect_equal(result$N, N)
  expect_equal(result$T, Tp)
  expect_equal(result$S, S)

  # event-time ATT should be close to 2
  expect_length(result$att_event, S)
  expect_true(all(abs(result$att_event - 2) < 2))  # rough check

  # overall ATT should be close to 2
  expect_true(abs(result$att_overall - 2) < 2)

  # min eigenvalue should be positive
  expect_true(result$min_eigenvalue > 0)

  # inference available (T > S)
  expect_true(!is.na(result$p_value))
  expect_true(!anyNA(result$ci_lower_event))
})


test_that("ssc handles T < S (no inference) gracefully", {
  set.seed(456)
  N <- 10; Tp <- 5; S <- 8

  Y0 <- matrix(rnorm(N * (Tp + S)), N)
  D  <- matrix(0, N, Tp + S)
  D[1, (Tp + 1):(Tp + S)] <- 1
  D[2, (Tp + 3):(Tp + S)] <- 1

  Y <- Y0 + D

  expect_message(result <- ssc(Y, D), "not enough pre-treatment")

  expect_s3_class(result, "ssc")
  expect_true(is.na(result$p_value))
  expect_true(all(is.na(result$ci_lower_event)))
  expect_length(result$att_event, S)  # point estimates still available
})


test_that("ssc errors on invalid input", {
  # mismatched dimensions
  expect_error(ssc(matrix(1, 3, 5), matrix(1, 4, 5)), "same dimensions")

  # no treatment
  expect_error(ssc(matrix(1, 3, 5), matrix(0, 3, 5)), "No treatment")

  # non-binary D
  expect_error(ssc(matrix(1, 3, 5), matrix(2, 3, 5)), "binary")
})


## ---- print / summary / plot methods ----

test_that("print.ssc runs without error", {
  set.seed(789)
  N <- 8; Tp <- 20; S <- 3
  Y <- matrix(rnorm(N * (Tp + S)), N)
  D <- matrix(0, N, Tp + S)
  D[1, (Tp + 1):(Tp + S)] <- 1
  D[2, (Tp + 2):(Tp + S)] <- 1
  Y <- Y + D

  result <- ssc(Y, D)
  expect_output(print(result), "Staggered Synthetic Control")
})


test_that("summary.ssc runs without error", {
  set.seed(789)
  N <- 8; Tp <- 20; S <- 3
  Y <- matrix(rnorm(N * (Tp + S)), N)
  D <- matrix(0, N, Tp + S)
  D[1, (Tp + 1):(Tp + S)] <- 1
  D[2, (Tp + 2):(Tp + S)] <- 1
  Y <- Y + D

  result <- ssc(Y, D)
  expect_output(summary(result), "Summary")
})
