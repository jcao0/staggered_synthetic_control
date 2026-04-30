## Figure 3: ATT estimates -- SSC and GSC overlay.
## Run from the project root directory.
##
## Input:  file.path(tempdir(), "cleaned_data", "psrm_crime_data.csv")
##         file.path(tempdir(), "cleaned_data", "psrm_cartel_data.csv")
##         data/results_gsc.csv
## Output: file.path(tempdir(), "Figure3_att_estimates.png")
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
gsc_results <- read.csv("data/results_gsc.csv")

COL_SSC <- "#E8593C"
COL_GSC <- "#3CBFBF"

build_panel <- function(data, unit_col, time_col, treatment_col,
                        outcome_col, time_filter = NULL,
                        outcome_name, label) {

  if (!is.null(time_filter))
    data <- data[eval(time_filter, envir = data), ]

  mat <- panel_to_matrices(data,
                           unit      = unit_col,
                           time      = time_col,
                           outcome   = outcome_col,
                           treatment = treatment_col)

  result <- ssc(mat$Y, mat$D, alpha = 0.05)
  cat(sprintf("%s:  N=%d  T=%d  S=%d  ATT=%.4f",
              outcome_name, result$N, result$T, result$S, result$att_overall))
  if (!is.na(result$ci_lower_overall))
    cat(sprintf("  [%.4f, %.4f]  p=%.4f",
                result$ci_lower_overall, result$ci_upper_overall,
                result$p_value))
  cat("\n")

  p <- plot(result, main = "")

  gsc_sub     <- gsc_results[gsc_results$outcome == outcome_name, ]
  gsc_sub$et  <- gsc_sub$event_time - 1L

  p <- p +
    ggplot2::geom_ribbon(
      data = gsc_sub,
      ggplot2::aes(x = et, ymin = CI.lower, ymax = CI.upper, fill = "GSC"),
      alpha = 0.15, inherit.aes = FALSE) +
    ggplot2::geom_line(
      data = gsc_sub,
      ggplot2::aes(x = et, y = ATT, colour = "GSC"),
      linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(name = NULL,
      values = c("SSC" = COL_SSC, "GSC" = COL_GSC),
      guide = ggplot2::guide_legend(
        override.aes = list(linetype = c("dashed", "solid")))) +
    ggplot2::scale_fill_manual(name = NULL,
      values = c("SSC" = COL_SSC, "GSC" = COL_GSC),
      guide = "none") +
    ggplot2::labs(x = "event time", y = "ATT estimates", caption = label) +
    ggplot2::theme(
      plot.caption         = ggplot2::element_text(hjust = 0.5, size = 11,
                               margin = ggplot2::margin(t = 4)),
      legend.position      = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank(),
      legend.key.width     = ggplot2::unit(3, "cm"),
      legend.text          = ggplot2::element_text(size = 12)
    )

  p
}

p_a <- build_panel(crime_data, "idunico", "time", "Policial", "hom_all_rate",
                   expression(time < 253), "hom_all_rate", "(a) Homicide rate")

p_b <- build_panel(crime_data, "idunico", "time", "Policial", "hom_ym_rate",
                   expression(time < 253), "hom_ym_rate", "(b) Cartel-related homicide rate")

p_c <- build_panel(crime_data, "idunico", "time", "Policial", "theft_violent_rate",
                   expression(time >= 133), "theft_violent_rate", "(c) Violent theft rate")

p_d <- build_panel(crime_data, "idunico", "time", "Policial", "theft_nonviolent_rate",
                   expression(time >= 133), "theft_nonviolent_rate", "(d) Nonviolent theft rate")

p_e <- build_panel(cartel_data, "idunico", "Year", "policial", "presence_strength",
                   NULL, "presence_strength", "(e) Cartel strength")

p_f <- build_panel(cartel_data, "idunico", "Year", "policial", "co_num",
                   NULL, "co_num", "(f) Number of cartels")

p_g <- build_panel(cartel_data, "idunico", "Year", "policial", "war",
                   NULL, "war", "(g) Cartel war")

blank <- grid::nullGrob()

fig <- gridExtra::arrangeGrob(
  p_a, p_b,
  p_c, p_d,
  p_e, p_f,
  p_g, blank,
  ncol = 2
)

fig3_path <- file.path(out_dir, "Figure3_att_estimates.png")
ggplot2::ggsave(fig3_path, fig, width = 10, height = 18, dpi = 150)

cat("\nFigure saved to:", fig3_path, "\n")
