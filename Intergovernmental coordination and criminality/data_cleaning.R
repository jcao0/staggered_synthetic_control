## Empirical Study
## Replication Materials

# This R script processes data for the empirical replication of:
# "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico"

# Raw data is provided by the original authors at:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FPBTTDM

# Processing steps:
# 1. Exclude observations removed in the original study.
# 2. Save processed data as CSV files.

# Note: All data processing steps strictly follow the procedures described in the original paper.


## %% set environment %%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls())

# record start time
begin.time <- Sys.time()

# set working directory to where this script is located
setwd("/Users/replication_files/Intergovernmental coordination and criminality")

# load packages
# install.packages("dplyr")
library(dplyr)

## %% process data %%%%%%%%%%%%%%%%%%%%%%%%%%

# monthly data for outcome related to theft and homicide
load("raw_data/psrm_crime_data.RData")
psrm_crime_data <- mu_policial_monthly  |> ungroup()
# remove units that implement then remove reform
psrm_crime_data <- psrm_crime_data[!(psrm_crime_data$idunico %in% c("11016","11010","11027","11033","11005")),] 
# remove Leon, an outlier
psrm_crime_data <- psrm_crime_data[psrm_crime_data$idunico!=11020,]



# annual data for outcome related to cartel presence
load("raw_data/psrm_cartel_data.RData")
# remove units that get rid of reform
psrm_cartel_data <- mu_policial_annual |> ungroup()
psrm_cartel_data <- psrm_cartel_data[!(psrm_cartel_data$idunico %in% c("11027","11010","11016","11033")),] 
# remove Leon, an outlier
psrm_cartel_data <- psrm_cartel_data[psrm_cartel_data$idunico!=11020,]


## %% save data %%%%%%%%%%%%%%%%%%%%%%%%%%
## save the processed data in each of the folders producing results

# find folders producing results in the current directory
folders <- list.dirs(path = ".", recursive = FALSE)
folders <- folders[!grepl("data", folders)] # exclude "raw_data" folder

# save the processed data in a new "cleaned_data" folder within each of those folders
for (folder in folders) {
  dir.create(file.path(folder, "cleaned_data"), showWarnings = FALSE)
  write.csv(psrm_crime_data, file.path(folder, "cleaned_data/psrm_crime_data.csv"), row.names = FALSE)
  write.csv(psrm_cartel_data, file.path(folder, "cleaned_data/psrm_cartel_data.csv"), row.names = FALSE)
}

# running time
time <- Sys.time() - begin.time
cat("Data cleaning completed in ", time, " ", attr(time, "units"), ".\n", sep = "")


