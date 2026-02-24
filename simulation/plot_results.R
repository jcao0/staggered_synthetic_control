## Simulation Study
## Replication Materials


## This scripts plots the simulation results from different methods
## Main steps:
## 1. Compute summary statistics (RMSE of ATT estimates, elapsed time)
## 2. Plot RMSE


# %% Set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))

# install.packages(c("ggplot2", "patchwork"))
library(ggplot2)
library(patchwork)
source('functions/plot.R') # load plotting functions


# %% Summary results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load simulation results
files <- list.files("output/", pattern = "\\.csv$")
files <- files[order(as.numeric(sub(".*_r([0-9]+)_.*", "\\1", files)), as.numeric(sub(".*T([0-9]+)_.*", "\\1", files)))]

case_names <- unique(sapply(strsplit(files, "_"), function(x) paste(x[3], x[4], x[5], sep = "_")))
methods <- unique(sapply(strsplit(files, "_"), function(x) { fname <- x[length(x)]; sub("\\.csv$", "", fname)}))

# create a nested list to store summary results: case -> method -> summary
summary = list()
for (case in case_names) {
  summary[[case]] <- list()
  for (method in methods) {
    cat("Processing case:", case, "method:", method, "\n")
    # load sim results
    filename <- paste0("output/sim_results_", case, "_", method , ".csv")
    result <- read.csv(filename, check.names = FALSE)
    
    # compute summary statistics
    att.e <- as.matrix(result[, !names(result) %in% "time_sec"]) # exclude time column
    att.e <- att.e[,1:7] #v2
    att_true <- c(1:ncol(att.e)) # true ATT values equal to event times + 1


    rmse_att_e <- sapply(att_true, function(n) if (any(!is.na(att.e[,n]))) sqrt(mean((att.e[,n] - n)^2, na.rm = TRUE)) else NA)
    elapsed_time <- sum(result$`time_sec`, na.rm = TRUE); elapsed_time <- rep(elapsed_time, times = nrow(att.e)); elapsed_time[is.na(att.e)] <- NA

    n_sims <- nrow(att.e)

    # store summary results
    summary[[case]][[method]] <- list(
      rmse = rmse_att_e,
      time = elapsed_time,
      n_sims = n_sims
    )
  }
}

# %% Plot rmse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare data for plotting
plot_data <- data.frame()
for (case in case_names) {
  for (method in methods) {
    rmse_values <- summary[[case]][[method]]$rmse
    temp_df <- data.frame(
      case = case,
      method = ifelse (method == "asy", "ppsc", method), # rename "asy" to "ppsc" for clarity in plots
      event_time = seq_along(rmse_values) -1, # event times start from 0
      rmse = rmse_values
    )
    plot_data <- rbind(plot_data, temp_df)
  }
}

# plot rmse for each case
plot_list <- lapply(case_names, function(cname) {
  plot_rmse_case(plot_data, cname) # plot function from plot.R
})
names(plot_list) <- case_names

save_file = 'output/Figure1_simulation_results.png'
combined = combine_plots(plot_list, save_file = save_file)  # plot combining function from plot.R
cat("\nSimulation figures saved in 'output' folder.\n")
print(combined)


# print elapsed time for simulation
cat("\nElapsed time of each method:\n")
for (case in case_names) {
  cat("Case:", case, "\n")
  for (method in methods) {
    elapsed_time <- unique(na.omit(summary[[case]][[method]]$time))
    cat("  Method:", method, "- Elapsed time (secs):", elapsed_time, "\n")
  }
}


