# load in packages
library(tidyverse)
library(cowplot)
library(here)

# load in functions 
source(here("code", "processing", "calculation_funs.R"))
source(here("code", "processing", "data_import_funs.R"))

# build imprinting HA groups: 1 is H1N1 and H2N2, 2 is H3N2
get_p_unprot <- function(Wmat, H, group) {
  if (group == 1)
    p_imprint <- Wmat$H1N1_probs + Wmat$H2N2_probs
  else if (group == 2)
    p_imprint <-  Wmat$H3N2_probs
  else
    stop("please select valid group (1 or 2)")
  
  p_imprint*H + (1 - p_imprint)
}

# run imprinting probs
setwd(here::here())
getwd()
Wmat <- get_imprinting_probabilities(
  observation_years = 2024
)


# add in the failure rate parameter
p_get_anyway <-  # the H parameter
  
p_unprot_H5 <- get_p_unprot(Wmat, H=0.25, group=1)
p_unprot_H7 <- get_p_unprot(Wmat, H=0.24, group=2)

# load in the county risk (age and comorbs) data
age_risk <- read_csv(here("data", "processed-data", "county_risk_weighted_with_dups.csv")) %>%
  rename(county_fips = county_fips.x)

# checking
age_risk |>
  filter(county_fips == "12119")


### age_risk_long <- age_risk |> do not use
#  pivot_longer(
#    cols = starts_with("avg_"),
#    names_to = "age_group",
#    values_to = "prob_severe_inf"
#  ) |>
#  mutate(
#    age_min = parse_number(str_extract(age_group, "(?<=_)[\\d]+")),
#    age_max = parse_number(str_extract(age_group, "(?<=_)[\\d]+$"))
#  )
  
age_risk_long <- age_risk %>%
    rename(prob_severe_inf = risk) %>%
    mutate(
      age_label = str_remove(age_group, "^avg_"),
      
      age_min = parse_number(
        str_extract(age_label, "^[0-9]+")),
      
      age_max = case_when(
        str_detect(age_label, "\\+$") ~ 107,
        TRUE ~ parse_number(
          str_extract(age_label, "[0-9]+$"))))


# checking
age_risk_long %>%
    filter(county_fips == "12119") %>%
    select(
      age_group,
      age_min,
      age_max,
      prob_severe_inf
    )
  

# takes a bit longer to run
age_risk_comp <- age_risk_long |> 
  pivot_longer(age_min:age_max, names_to = "name", values_to = "value") |> 
  mutate(
    name = "age", 
    age_group = fct_inorder(age_group),
    value = case_match(value, 
                       0.5 ~ 1, 
                       NA ~ 107, 
                       .default = value)
  ) |> 
  group_by(county_fips, age_group) |> 
  complete(value = first(value):max(value)) |> 
  ungroup()

# checking
age_risk_comp |>
  filter(
    county_fips == "12119",
    value >= 60) |>
  select(
    age_group,
    value,
    prob_severe_inf) |>
  print(n = Inf)


#age_risk_comp <- age_risk_comp |> 
#  fill(prob_severe_inf, name, .direction="up") |> 
#  select(fips=county_fips, age=value, prob_severe_inf) |> 
#  distinct(fips, age, .keep_all=TRUE) # removes the few counties with two entries per fips

age_risk_comp_filled <- age_risk_comp |>
  group_by(county_fips, age_group) |>
  arrange(value, .by_group = TRUE) |>
  fill(
    prob_severe_inf,
    name,
    .direction = "up") |>
  ungroup()

# maps for imprinting distribution
load(here("data", "processed-data", "county_fine_age_pop.rda")) # read age distributions by county

county_age_dist <- county_fine_age_pop |>
  ungroup() |> 
  # add 0 probs for ages over the cutoff: OK because of sum and since all these ages imprinted same
  complete(
    nesting(birthyear=1924:1918, age=106:100), nesting(fips, county_name), 
    fill=list(pop=0, frac_pop=0)
  ) |>
  arrange(fips, desc(birthyear))

county_risk_for_plot <- county_age_dist |> 
  mutate(fips = as.character(fips)) |> 
  left_join(
    age_risk_comp |> mutate(fips = as.character(fips)),
    by = c("fips", "age")
  ) |> 
  group_by(fips) |> 
  # weighted average of prob_unprotected over the age distribution
  summarise(
    risk_unprot_H5 = sum(frac_pop * p_unprot_H5),
    risk_severe_H5 = sum(frac_pop * p_unprot_H5 * prob_severe_inf),
    risk_comorb = sum(frac_pop * prob_severe_inf),
    risk_unprot_H7 = sum(frac_pop * p_unprot_H7),
    risk_severe_H7 = sum(frac_pop * p_unprot_H7 * prob_severe_inf),
  )

county_risk_for_plot |>
  filter(fips == "12119")

# save data to use in plots
write_csv(county_risk_for_plot, here("data", "processed-data", "map_data_numbers.csv"))
