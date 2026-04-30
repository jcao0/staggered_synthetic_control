## Run all scripts to generate Figure 2, Figure 3, and Table 1.
## Set the working directory to the project root before running.
##
## Outputs are written to tempdir() by default.  To save to a custom
## location, set REPLICATION_OUT before sourcing this script:
##   Sys.setenv(REPLICATION_OUT = "path/to/output")

## Set working directory to this folder before running, e.g.:
## setwd("path/to/replication")

timing <- list()

# Step 1: clean data
cat("== Step 1: Data cleaning ==\n")
t0 <- proc.time()
source("data_clean.R")
timing[["data_clean.R"]] <- (proc.time() - t0)[["elapsed"]]

# Step 2: Figure 2 -- treatment scheme
cat("\n== Step 2: Figure 2 -- treatment scheme ==\n")
t0 <- proc.time()
source("figure2_treatment_scheme.R")
timing[["figure2_treatment_scheme.R"]] <- (proc.time() - t0)[["elapsed"]]

# Step 3: Figure 3 -- ATT estimates
cat("\n== Step 3: Figure 3 -- ATT estimates ==\n")
t0 <- proc.time()
source("figure3_att_estimates.R")
timing[["figure3_att_estimates.R"]] <- (proc.time() - t0)[["elapsed"]]

# Step 4: Table 1 -- smallest eigenvalue
cat("\n== Step 4: Table 1 -- smallest eigenvalue ==\n")
t0 <- proc.time()
source("table1_eigenvalue.R")
timing[["table1_eigenvalue.R"]] <- (proc.time() - t0)[["elapsed"]]

out_dir <- Sys.getenv("REPLICATION_OUT", unset = tempdir())
cat("\n== Done. Outputs saved in", out_dir, "==\n")

# Write running times to txt
lines <- c("Running times (seconds)", "========================")
for (nm in names(timing)) {
  val <- timing[[nm]]
  if (is.na(val)) {
    lines <- c(lines, sprintf("%-35s  skipped (output already exists)", nm))
  } else {
    lines <- c(lines, sprintf("%-35s  %.2f s", nm, val))
  }
}
times_path <- file.path(out_dir, "running_times.txt")
writeLines(lines, times_path)
cat("Running times saved to:", times_path, "\n")
