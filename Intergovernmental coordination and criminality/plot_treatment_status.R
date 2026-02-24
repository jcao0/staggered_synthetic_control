## %% Set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))
library(panelView)
library(ggplot2)
library(cowplot)
library(dplyr)

## %% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Treated status: homicide
psrm_crime_data <- read.csv("data/processed/psrm_crime_data.csv")
psrm_crime_data <- psrm_crime_data[(psrm_crime_data$Year<2021),]
psrm_crime_data$time <- psrm_crime_data$time - min(psrm_crime_data$time) + 1
p_1 = panelview(hom_all_rate ~ Policial, data = psrm_crime_data, 
          index = c("idunico","time"), 
          by.timing = TRUE, main = "",ylab = "municipalities", xlab = "time", axis.lab = c("time"), legend.labs = c("untreated", "treated"), color =  c("#0099991A","#009999"), background = "white",cex.legend = 22,cex.lab = 22,cex.axis= 22,axis.lab.gap = c(12)) 
n_units <- length(unique(psrm_crime_data$idunico))
n_times <- length(unique(psrm_crime_data$time))
cat("\n Outcomes: homicide | The ratio of time periods to units is:", n_times/n_units, "\n")


#Treated status: theft
psrm_crime_data <- read.csv("data/processed/psrm_crime_data.csv")
psrm_crime_data <- psrm_crime_data[(psrm_crime_data$Year>2010),]
psrm_crime_data$time <- psrm_crime_data$time - min(psrm_crime_data$time) + 1
p_2 = panelview(theft_violent_rate ~ Policial, data = psrm_crime_data, 
          index = c("idunico","time"), 
          by.timing = TRUE, main = "",ylab = "municipalities", xlab = "time", axis.lab = c("time"), legend.labs = c("untreated", "treated"), color =  c("#0099991A","#009999"), background = "white",cex.legend = 22,cex.lab = 22,cex.axis= 22,axis.lab.gap = c(12)) 
n_units <- length(unique(psrm_crime_data$idunico))
n_times <- length(unique(psrm_crime_data$time))
cat("\n Outcomes: theft | The ratio of time periods to units is:", n_times/n_units, "\n")


#Treated status: cartel
psrm_cartel_data <- read.csv("data/processed/psrm_cartel_data.csv")
psrm_cartel_data$Year <- psrm_cartel_data$Year - min(psrm_cartel_data$Year) + 1
p_3 = panelview(presence_strength ~ policial, data = psrm_cartel_data, 
          index = c("idunico","Year"), 
          by.timing = TRUE, main = "", ylab = "municipalities", xlab = "time", axis.lab = c("time"), legend.labs = c("untreated", "treated"), color =  c("#0099991A","#009999"), background = "white", cex.legend = 22,cex.lab = 22,cex.axis= 22, axis.lab.gap = c(3))
n_units <- length(unique(psrm_cartel_data$idunico_num))
n_times <- length(unique(psrm_cartel_data$Year))
cat("\n Outcomes: cartel | The ratio of time periods to units is:", n_times/n_units, "\n")



## %% Combine plots and save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plots <- list(p_1, p_2, p_3)
captions <- c("(a) Homicide and cartel-related homicide rate",
              "(b) Violent and nonviolent theft rate",
              "(c) Cartel strength, number and war")


x_ranges <- sapply(plots, function(p) {
  diff(range(p$data$period, na.rm = TRUE))
})
scaled_widths <- sqrt(x_ranges / max(x_ranges))
scaled_widths <- scaled_widths[2:3]/sum(scaled_widths[2:3]) * 0.95 # scale to fit within 95% of the width for the second row

# Rescale individual plots with your theme + captions
plots_rescaled <- Map(function(p, cap) {
  p <- p + 
    theme(
      legend.box.spacing = unit(1, "pt"), # margin between legend and plot
      legend.box.margin = margin(t = 4, r = 0, b = 15, l = 0), # margin outside legend box
      legend.margin = margin(t = 2, r = 0, b = 6, l = 0),  # margin inside legend box
      legend.key.size = unit(0.5, "cm"),
       legend.text = element_text(margin = margin(l = 5,r = 5))   # ← space between key and text
    )
  
  ggdraw() +
    draw_plot(p, width = 1, x = 0) +   # full width
    draw_label(
      cap,
      x = 0.5,  
      y = 0,   
      hjust = 0.47,
      vjust = 1,
      size = 22
    ) +
    theme(plot.margin = margin(t = 40, r = 25, b = 40, l = 25)) # margin around the subplot to accommodate caption
}, plots, captions)

# First row: p1 (full width)
row1 <- ggdraw() + draw_plot(plots_rescaled[[1]], width = 1, x = 0)

# Second row: p2 + p3 side by side, but full width combined
row2 <- ggdraw() +
  draw_plot(plots_rescaled[[2]], x = 0, width = scaled_widths[1]) +
  draw_plot(plots_rescaled[[3]], x = 0.95 - scaled_widths[2], width = scaled_widths[2])

# Combine vertically
combined <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

ggsave("results/Figure2_treatment_status.png", width = 15, height = 10, dpi = 300,  bg = "white")

