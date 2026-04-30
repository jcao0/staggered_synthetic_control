## Figure 2: Treatment scheme plots.
## Run from the project root directory.
##
## Input:  file.path(tempdir(), "cleaned_data", "psrm_crime_data.csv")
##         file.path(tempdir(), "cleaned_data", "psrm_cartel_data.csv")
## Output: file.path(tempdir(), "Figure2_treatment_scheme.png")
##
## Set REPLICATION_OUT to redirect both input and output, e.g.:
##   Sys.setenv(REPLICATION_OUT = "path/to/output")

library(panelView)
library(ggplot2)
library(cowplot)
library(dplyr)

out_dir  <- Sys.getenv("REPLICATION_OUT", unset = tempdir())
data_dir <- file.path(out_dir, "cleaned_data")

# homicide
psrm_crime_data <- read.csv(file.path(data_dir, "psrm_crime_data.csv"))
psrm_crime_data <- psrm_crime_data[(psrm_crime_data$Year < 2021), ]
psrm_crime_data$time <- psrm_crime_data$time - min(psrm_crime_data$time) + 1
p_1 <- panelview(hom_all_rate ~ Policial, data = psrm_crime_data,
          index = c("idunico", "time"),
          by.timing = TRUE, main = "", ylab = "municipalities", xlab = "time",
          axis.lab = c("time"), legend.labs = c("untreated", "treated"),
          color = c("#0099991A", "#009999"), background = "white",
          cex.legend = 22, cex.lab = 22, cex.axis = 22, axis.lab.gap = c(12))
n_units <- length(unique(psrm_crime_data$idunico))
n_times <- length(unique(psrm_crime_data$time))
cat("\nOutcomes: homicide | T/N ratio:", n_times / n_units, "\n")

# theft
psrm_crime_data <- read.csv(file.path(data_dir, "psrm_crime_data.csv"))
psrm_crime_data <- psrm_crime_data[(psrm_crime_data$Year > 2010), ]
psrm_crime_data$time <- psrm_crime_data$time - min(psrm_crime_data$time) + 1
p_2 <- panelview(theft_violent_rate ~ Policial, data = psrm_crime_data,
          index = c("idunico", "time"),
          by.timing = TRUE, main = "", ylab = "municipalities", xlab = "time",
          axis.lab = c("time"), legend.labs = c("untreated", "treated"),
          color = c("#0099991A", "#009999"), background = "white",
          cex.legend = 22, cex.lab = 22, cex.axis = 22, axis.lab.gap = c(12))
n_units <- length(unique(psrm_crime_data$idunico))
n_times <- length(unique(psrm_crime_data$time))
cat("\nOutcomes: theft | T/N ratio:", n_times / n_units, "\n")

# cartel
psrm_cartel_data <- read.csv(file.path(data_dir, "psrm_cartel_data.csv"))
psrm_cartel_data$Year <- psrm_cartel_data$Year - min(psrm_cartel_data$Year) + 1
p_3 <- panelview(presence_strength ~ policial, data = psrm_cartel_data,
          index = c("idunico", "Year"),
          by.timing = TRUE, main = "", ylab = "municipalities", xlab = "time",
          axis.lab = c("time"), legend.labs = c("untreated", "treated"),
          color = c("#0099991A", "#009999"), background = "white",
          cex.legend = 22, cex.lab = 22, cex.axis = 22, axis.lab.gap = c(3))
n_units <- length(unique(psrm_cartel_data$idunico_num))
n_times <- length(unique(psrm_cartel_data$Year))
cat("\nOutcomes: cartel | T/N ratio:", n_times / n_units, "\n")

plots    <- list(p_1, p_2, p_3)
captions <- c("(a) Homicide and cartel-related homicide rate",
              "(b) Violent and nonviolent theft rate",
              "(c) Cartel strength, number and war")

x_ranges     <- sapply(plots, function(p) diff(range(p$data$period, na.rm = TRUE)))
scaled_widths <- sqrt(x_ranges / max(x_ranges))
scaled_widths <- scaled_widths[2:3] / sum(scaled_widths[2:3]) * 0.95

plots_rescaled <- Map(function(p, cap) {
  p <- p +
    theme(
      legend.box.spacing = unit(1, "pt"),
      legend.box.margin  = margin(t = 4, r = 0, b = 15, l = 0),
      legend.margin      = margin(t = 2, r = 0, b = 6,  l = 0),
      legend.key.size    = unit(0.5, "cm"),
      legend.text        = element_text(margin = margin(l = 5, r = 5))
    )
  ggdraw() +
    draw_plot(p, width = 1, x = 0) +
    draw_label(cap, x = 0.5, y = 0, hjust = 0.47, vjust = 1, size = 22) +
    theme(plot.margin = margin(t = 40, r = 25, b = 40, l = 25))
}, plots, captions)

row1 <- ggdraw() + draw_plot(plots_rescaled[[1]], width = 1, x = 0)
row2 <- ggdraw() +
  draw_plot(plots_rescaled[[2]], x = 0,                      width = scaled_widths[1]) +
  draw_plot(plots_rescaled[[3]], x = 0.95 - scaled_widths[2], width = scaled_widths[2])
combined <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

fig2_path <- file.path(out_dir, "Figure2_treatment_scheme.png")
ggsave(fig2_path, combined, width = 15, height = 10, dpi = 300, bg = "white")

cat("\nFigure saved to:", fig2_path, "\n")
