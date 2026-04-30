## Data cleaning for Cao, Lu, and Wu (2020) replication.
## Run from the project root directory.
##
## Input:  data/rawdata/psrm_crime_data.RData
##         data/rawdata/psrm_cartel_data.RData
## Output: file.path(tempdir(), "cleaned_data", "psrm_crime_data.csv")
##         file.path(tempdir(), "cleaned_data", "psrm_cartel_data.csv")
##
## To save outputs to a custom location, set the environment variable
## REPLICATION_OUT before sourcing this script, e.g.:
##   Sys.setenv(REPLICATION_OUT = "path/to/output")

library(dplyr)

# monthly crime data
load("data/rawdata/psrm_crime_data.RData")
psrm_crime_data <- mu_policial_monthly |> ungroup()
psrm_crime_data <- psrm_crime_data[!(psrm_crime_data$idunico %in%
                     c("11016", "11010", "11027", "11033", "11005")), ]
psrm_crime_data <- psrm_crime_data[psrm_crime_data$idunico != 11020, ]

# annual cartel data
load("data/rawdata/psrm_cartel_data.RData")
psrm_cartel_data <- mu_policial_annual |> ungroup()
psrm_cartel_data <- psrm_cartel_data[!(psrm_cartel_data$idunico %in%
                      c("11027", "11010", "11016", "11033")), ]
psrm_cartel_data <- psrm_cartel_data[psrm_cartel_data$idunico != 11020, ]

out_dir <- Sys.getenv("REPLICATION_OUT", unset = tempdir())
data_out <- file.path(out_dir, "cleaned_data")
dir.create(data_out, showWarnings = FALSE, recursive = TRUE)
write.csv(psrm_crime_data,  file.path(data_out, "psrm_crime_data.csv"),  row.names = FALSE)
write.csv(psrm_cartel_data, file.path(data_out, "psrm_cartel_data.csv"), row.names = FALSE)

cat("Cleaned data saved to:", data_out, "\n")
