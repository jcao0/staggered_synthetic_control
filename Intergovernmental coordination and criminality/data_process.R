######################
# This R script processes data for the empirical replication of:
# "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico"

# Raw data is provided by the original authors at:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FPBTTDM

# Processing steps:
# 1. Exclude observations removed in the original study.
# 2. Save processed data as CSV files for use in both R and Matlab implementations.

# Note: All data processing steps strictly follow the procedures described in the original paper.
######################
rm(list=ls())
# Set working directory ----
setwd("/Users/wuhang/Desktop/metric/sc/replication_files")


# Packages -----
library(pacman)
# p_load(gsynth, dplyr, ggplot2, stargazer)
p_load(dplyr)

## %% monthly data for outcome related to theft and homicide %%%%%%%%%%%%%%%%%%%%%%%%%%
load("data/raw/psrm_crime_data.RData")
psrm_crime_data <- mu_policial_monthly  |> ungroup()
# remove units that implement then remove reform
psrm_crime_data <- psrm_crime_data[!(psrm_crime_data$idunico %in% c("11016","11010","11027","11033","11005")),] 
# remove Leon, an outlier
psrm_crime_data <- psrm_crime_data[psrm_crime_data$idunico!=11020,]
# store as csv
write.csv(psrm_crime_data, 'data/processed/psrm_crime_data.csv', row.names = FALSE)

## %% annual data for outcome related to cartel presence %%%%%%%%%%%%%%%%%%%%%%%%%%
load("data/raw/psrm_cartel_data.RData")
# remove units that get rid of reform
psrm_cartel_data <- mu_policial_annual |> ungroup()
psrm_cartel_data <- psrm_cartel_data[!(psrm_cartel_data$idunico %in% c("11027","11010","11016","11033")),] 
# remove Leon, an outlier
psrm_cartel_data <- psrm_cartel_data[psrm_cartel_data$idunico!=11020,]
# store as csv
write.csv(psrm_cartel_data, 'data/processed/psrm_cartel_data.csv', row.names = FALSE)

