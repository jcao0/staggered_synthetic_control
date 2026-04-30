## Table 1: SSC estimation summary (incl. smallest eigenvalue) for all outcomes.
## Run from the project root directory.
##
## Input:  file.path(tempdir(), "cleaned_data", "psrm_crime_data.csv")
##         file.path(tempdir(), "cleaned_data", "psrm_cartel_data.csv")
## Output: file.path(tempdir(), "Table1_eigenvalue.csv")
##
## Set REPLICATION_OUT to redirect input/output, e.g.:
##   Sys.setenv(REPLICATION_OUT = "path/to/output")

source("function/synthetic_control.R")
source("function/ssc.R")
source("function/utils.R")
source("function/methods.R")

out_dir  <- Sys.getenv("REPLICATION_OUT", unset = tempdir())
data_dir <- file.path(out_dir, "cleaned_data")

crime_data  <- read.csv(file.path(data_dir, "psrm_crime_data.csv"))
cartel_data <- read.csv(file.path(data_dir, "psrm_cartel_data.csv"))

run_ssc <- function(data, unit_col, time_col, treatment_col,
                    outcome_col, time_filter = NULL, outcome_label) {
  if (!is.null(time_filter))
    data <- data[eval(time_filter, envir = data), ]

  mat    <- panel_to_matrices(data,
                              unit      = unit_col,
                              time      = time_col,
                              outcome   = outcome_col,
                              treatment = treatment_col)
  result <- ssc(mat$Y, mat$D, alpha = 0.05)

  data.frame(
    Outcome          = outcome_label,
    N                = result$N,
    T_pre            = result$T,
    S_post           = result$S,
    Min_Eigenvalue   = round(result$min_eigenvalue, 6),
    ATT              = round(result$att_overall, 4),
    CI_lower         = round(result$ci_lower_overall, 4),
    CI_upper         = round(result$ci_upper_overall, 4),
    p_value          = round(result$p_value, 4),
    stringsAsFactors = FALSE
  )
}

rows <- list(
  run_ssc(crime_data,  "idunico", "time", "Policial", "hom_all_rate",
          expression(time < 253),  "(a) Homicide rate"),
  run_ssc(crime_data,  "idunico", "time", "Policial", "hom_ym_rate",
          expression(time < 253),  "(b) Cartel-related homicide rate"),
  run_ssc(crime_data,  "idunico", "time", "Policial", "theft_violent_rate",
          expression(time >= 133), "(c) Violent theft rate"),
  run_ssc(crime_data,  "idunico", "time", "Policial", "theft_nonviolent_rate",
          expression(time >= 133), "(d) Nonviolent theft rate"),
  run_ssc(cartel_data, "idunico", "Year", "policial", "presence_strength",
          NULL,                    "(e) Cartel strength"),
  run_ssc(cartel_data, "idunico", "Year", "policial", "co_num",
          NULL,                    "(f) Number of cartels"),
  run_ssc(cartel_data, "idunico", "Year", "policial", "war",
          NULL,                    "(g) Cartel war")
)

tab <- do.call(rbind, rows)
tab <- tab[, c("Outcome", "Min_Eigenvalue")]

tab1_path <- file.path(out_dir, "Table1_eigenvalue.csv")
write.csv(tab, tab1_path, row.names = FALSE)

cat("Table 1:\n")
print(tab, row.names = FALSE)
cat("\nTable saved to:", tab1_path, "\n")
