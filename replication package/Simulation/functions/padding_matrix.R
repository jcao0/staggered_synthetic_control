# Two padding functions to unify matrix sizes when storing simulation results, 
# Matrix sizes may vary due to inconsistent event-time lengths across simulations.

# padding function for to unify matrix sizes: due to inconsistent event-time lengths
pad_matrix <- function(mat, nrow, ncol) {
  padded_mat <- matrix(NA, nrow = nrow, ncol = ncol)
  padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
  return(padded_mat)
}

pad_vector <- function(vec, nrow) c(vec, rep(NA, max(0, nrow - length(vec))))