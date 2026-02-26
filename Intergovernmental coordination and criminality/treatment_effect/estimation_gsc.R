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
setwd("/Users/replication_files/Intergovernmental coordination and criminality/treatment_effect")


# Packages 
# install.packages(c("gsynth", "dplyr"))
library(gsynth)
library(dplyr)




# %% Main analysis: outcome related to theft and homicide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Note: different time periods available for different outcomes

# load data
psrm_crime_data <- read.csv("cleaned_data/psrm_crime_data.csv")

## Theft
#theft violent
theft.viol.out <- gsynth(theft_violent_rate ~ Policial + n.cell + n.weak + n.strong + war + log_econin  + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO, 
                          data = psrm_crime_data[(psrm_crime_data$Year>2010),], EM = TRUE, 
                          index = c("idunico","time"), force = "two-way", 
                          CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                          inference = "parametric", nboots = 1000, 
                          parallel = T) 

#theft nonviolent
theft.nonviol.out <- gsynth(theft_nonviolent_rate ~ Policial + n.cell + n.weak + n.strong + war + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO, 
                             data = psrm_crime_data[(psrm_crime_data$Year>2010),], EM = TRUE, 
                             index = c("idunico","time"), force = "two-way", 
                             CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                             inference = "parametric", nboots = 1000, 
                             parallel = T)

## Homicides
#all homs
hom.out <- gsynth(hom_all_rate ~ Policial + n.cell + n.weak + n.strong + war + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO, 
                   data = psrm_crime_data[(psrm_crime_data$Year<2021),], EM = TRUE, 
                   index = c("idunico","time"), force = "two-way", 
                   CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                   inference = "parametric", nboots = 1000, 
                   parallel = T)
#homicides of young men
hom.ym.out <- gsynth(hom_ym_rate ~ Policial + n.cell + n.weak + n.strong + war + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO, 
                      data = psrm_crime_data[(psrm_crime_data$Year<2021),], EM = TRUE, 
                      index = c("idunico","time"), force = "two-way", 
                      CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                      inference = "parametric", nboots = 1000, 
                      parallel = T)



# %% Main analysis: outcome related to cartel presence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load data
psrm_cartel_data <- read.csv("cleaned_data/psrm_cartel_data.csv")


# cartel presence strength
out_strength <- gsynth(presence_strength ~ policial + log_pop + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO,# + mun_election_year, 
                       data = psrm_cartel_data, EM = TRUE, 
                       index = c("idunico","Year"), force = "two-way", 
                       CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                       inference = "parametric", nboots = 1000, 
                       parallel = T)

# number of cartels
out_number <- gsynth(co_num ~ policial + log_pop + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO,# + mun_election_year, 
                     data = psrm_cartel_data, EM = TRUE, 
                     index = c("idunico","Year"), force = "two-way", 
                     CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                     inference = "parametric", nboots = 1000, 
                     parallel = T)

# cartel war
out_war <- gsynth(war ~ policial + log_pop + log_econin + state_political_vulnerability + federal_political_vulnerability + state_federal_political_vulnerability + TT_PERSO,# + mun_election_year, 
                  data = psrm_cartel_data, EM = TRUE, 
                  index = c("idunico","Year"), force = "two-way", 
                  CV = TRUE, r = c(0,5), se = TRUE, cl = "idunico",
                  inference = "parametric", nboots = 1000, 
                  parallel = T)


# %% save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# collect ATT estimates for post-treatment periods

vtheft_att1 <- theft.viol.out$est.att[130:219,]
vtheft_att1 = data.frame(vtheft_att1)
vtheft_att1 = vtheft_att1 |> mutate(outcome='theft_violent_rate')

nvtheft_att1 <- theft.nonviol.out$est.att[130:219,]
nvtheft_att1 = data.frame(nvtheft_att1)
nvtheft_att1 = nvtheft_att1 |> mutate(outcome='theft_nonviolent_rate')

hom_att1 <- hom.out$est.att[195:272,]
hom_att1 = data.frame(hom_att1)
hom_att1 = hom_att1 |> mutate(outcome='hom_all_rate')

homym_att1 <- hom.ym.out$est.att[195:272,]
homym_att1 = data.frame(homym_att1)
homym_att1 = homym_att1 |> mutate(outcome='hom_ym_rate')


strength_att <- out_strength$est.att[17:23,] # extract ATT for post-treatment periods
strength_att = data.frame(strength_att)
strength_att = strength_att |> mutate(outcome='presence_strength')

number_att <- out_number$est.att[17:23,]
number_att = data.frame(number_att)
number_att = number_att |> mutate(outcome='co_num')

war_att <- out_war$est.att[17:23,]
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

