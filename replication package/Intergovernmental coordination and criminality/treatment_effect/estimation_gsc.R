## Empirical Study
## Replication Materials

# This R script replicates key results in Alcocer M (2025):
# "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico."

# It implements the Generalized Synthetic Control method (Xu, 2017) to estimate treatment effects on crime and cartel outcomes.
# The analysis is based on the original authors' code and data, available at:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FPBTTDM

## %% Set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls())

# start time
begin.time <- Sys.time()

# Set seed
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(1234)

# Set working directory to the folder where this script is located
get_script_dir <- function() {
  # 1) Run by RStudio / interactive IDE
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) {
      return(dirname(normalizePath(path)))
    }
  }
  
  # 2) Fallback
  return(NULL)
}

script_dir <- get_script_dir()

if (!is.null(script_dir)) {
  setwd(script_dir)
  cat("Working directory:", script_dir, "\n")
} else {
  message("Could not determine script directory.")
}


# Helper function to extract event-time average treatment effects
extract_event_time_att <- function(gsc_output, event_times, outcome) {
  est_att <- gsc_output$est.att
  target_rows <- as.character(seq_len(event_times))

  if (all(target_rows %in% rownames(est_att))) {
    return(est_att[target_rows, , drop = FALSE])
  }

  warning(
    "Could not identify event-time rows by row name for outcome ", outcome,
    ". Falling back to the final ", event_times, " rows of gsynth$est.att.",
    call. = FALSE
  )
  utils::tail(est_att, event_times)
}


## %% Main analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Load data
psrm_crime_data  <- read.csv("cleaned_data/psrm_crime_data.csv")
psrm_cartel_data <- read.csv("cleaned_data/psrm_cartel_data.csv")


# Model specifications
# Each list element defines:
#   - data      : dataset
#   - outcome   : dependent variable
#   - covars    : right-hand-side covariates (excluding outcome)
#   - subset    : row filter
#   - index     : panel identifiers

model_specs <- list(

  # ---------------- Crime outcomes ----------------
  theft.viol.out = list(
    data    = psrm_crime_data,
    outcome = "theft_violent_rate",
    covars  = c(
      "Policial", "n.cell", "n.weak", "n.strong", "war",
      "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(Year > 2010),
    index  = c("idunico", "time")
  ),

  theft.nonviol.out = list(
    data    = psrm_crime_data,
    outcome = "theft_nonviolent_rate",
    covars  = c(
      "Policial", "n.cell", "n.weak", "n.strong", "war",
      "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(Year > 2010),
    index  = c("idunico", "time")
  ),

  hom.out = list(
    data    = psrm_crime_data,
    outcome = "hom_all_rate",
    covars  = c(
      "Policial", "n.cell", "n.weak", "n.strong", "war",
      "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(Year < 2021),
    index  = c("idunico", "time")
  ),

  hom.ym.out = list(
    data    = psrm_crime_data,
    outcome = "hom_ym_rate",
    covars  = c(
      "Policial", "n.cell", "n.weak", "n.strong", "war",
      "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(Year < 2021),
    index  = c("idunico", "time")
  ),

  # ---------------- Cartel outcomes ----------------
  out_strength = list(
    data    = psrm_cartel_data,
    outcome = "presence_strength",
    covars  = c(
      "policial", "log_pop", "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(TRUE),
    index  = c("idunico", "Year")
  ),

  out_number = list(
    data    = psrm_cartel_data,
    outcome = "co_num",
    covars  = c(
      "policial", "log_pop", "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(TRUE),
    index  = c("idunico", "Year")
  ),

  out_war = list(
    data    = psrm_cartel_data,
    outcome = "war",
    covars  = c(
      "policial", "log_pop", "log_econin",
      "state_political_vulnerability",
      "federal_political_vulnerability",
      "state_federal_political_vulnerability",
      "TT_PERSO"
    ),
    subset = quote(TRUE),
    index  = c("idunico", "Year")
  )
)

# Common gsynth function arguments
common_args <- list(
  EM        = TRUE,
  force     = "two-way",
  CV        = TRUE,
  r         = c(0, 5),
  se        = TRUE,
  cl        = "idunico",
  inference = "parametric",
  nboots    = 1000,
  parallel  = TRUE
)

# Helper function to run gsynth model for different outcomes
run_gsynth <- function(spec) {

  # Build formula: outcome ~ x1 + x2 + ...
  fmla <- as.formula(
    paste(spec$outcome, "~", paste(spec$covars, collapse = " + "))
  )

  # Subset data
  data_sub <- subset(spec$data, eval(spec$subset, spec$data))

  # Combine arguments
  args <- c(
    list(
      formula = fmla,
      data    = data_sub,
      index   = spec$index
    ),
    common_args
  )

  # Run gsynth
  do.call(gsynth, args)
}


# Run gsynth models
results <- setNames({
  set.seed(1234) # Reset seed to produce identical results when running this section independently
  lapply(names(model_specs), function(model_name) {
    cat("Running:", model_name, "\n")
    run_gsynth(model_specs[[model_name]])
  })
}, names(model_specs))


# %% save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# collect ATT estimates for post-treatment periods

vtheft_att1 <- extract_event_time_att(results$theft.viol.out, 90, "theft_violent_rate")
vtheft_att1 = data.frame(vtheft_att1)
vtheft_att1 = vtheft_att1 |> mutate(outcome='theft_violent_rate')

nvtheft_att1 <- extract_event_time_att(results$theft.nonviol.out, 90, "theft_nonviolent_rate")
nvtheft_att1 = data.frame(nvtheft_att1)
nvtheft_att1 = nvtheft_att1 |> mutate(outcome='theft_nonviolent_rate')

hom_att1 <- extract_event_time_att(results$hom.out, 78, "hom_all_rate")
hom_att1 = data.frame(hom_att1)
hom_att1 = hom_att1 |> mutate(outcome='hom_all_rate')

homym_att1 <- extract_event_time_att(results$hom.ym.out, 78, "hom_ym_rate")
homym_att1 = data.frame(homym_att1)
homym_att1 = homym_att1 |> mutate(outcome='hom_ym_rate')


strength_att <- extract_event_time_att(results$out_strength, 7, "presence_strength")
strength_att = data.frame(strength_att)
strength_att = strength_att |> mutate(outcome='presence_strength')

number_att <- extract_event_time_att(results$out_number, 7, "co_num")
number_att = data.frame(number_att)
number_att = number_att |> mutate(outcome='co_num')

war_att <- extract_event_time_att(results$out_war, 7, "war")
war_att = data.frame(war_att)
war_att = war_att |> mutate(outcome='war')

# combine all results and save
results_gsc = rbind( hom_att1, homym_att1, vtheft_att1, nvtheft_att1, strength_att, number_att, war_att)
results_gsc = results_gsc |> group_by(outcome) |> mutate(event_time = row_number())
if (!dir.exists("output")) {
  dir.create("output")
}
write.csv(results_gsc, "output/results_gsc.csv", row.names = FALSE)

# save running time
time <- Sys.time() - begin.time
time <- as.numeric(time, units = "secs") # convert to seconds
write.table(time, file = "output/running_time_estimation_gsc.txt", row.names = FALSE, col.names = "Running Time (seconds)")
