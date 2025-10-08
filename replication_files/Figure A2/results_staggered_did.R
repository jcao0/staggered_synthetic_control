# install.packages("pacman")
rm(list=ls())
library(pacman)
p_load(did, dplyr)
# Data cleaning ------------------------------------------------

data <- read.csv("data_boardgendereige.csv")
treat_year <- data   |>  filter(treat == 1) |> group_by(unit) |> summarise(first.treat = min(time))
data <- data |> left_join(treat_year, by = c("unit")) |> mutate(first.treat = ifelse(is.na(first.treat), 0, first.treat)) |> arrange(unit, time)

# Event-study / dynamic effects ------------------------------------------------
out <- att_gt(
  yname = "fratio",
  gname = "first.treat",
  idname = "id",
  tname = "time",
  xformla = ~1,
  data = data,
  est_method = "reg"
)
summary(out)
es <- aggte(out, type = "dynamic")

# Export results ------------------------------------------------
# Dynamic effect estimates are in es$egt, es$att.egt, es$se.egt
es_df <- data.frame(
  event_time = es$egt,
  att = es$att.egt,
  se = es$se.egt
)

# Compute 95% confidence intervals
es_df$ci_lower <- es_df$att - 1.96 * es_df$se
es_df$ci_upper <- es_df$att + 1.96 * es_df$se

# View the data frame
print(es_df)


write.csv(es_df, "results_staggered_did.csv")
