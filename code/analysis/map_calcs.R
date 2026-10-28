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
  
age_risk_long <- age_risk |>
  rename(prob_severe_inf = risk) |>
  mutate(
    age_label = str_remove(age_group, "^avg_"),
    
    age_min = case_when(
      age_label == "0_0.5" ~ 0,
      age_label == "0.5_4" ~ 1,
      TRUE ~ parse_number(
        str_extract(age_label, "^[0-9]+"))),
    
    age_max = case_when(
      age_label == "0_0.5" ~ 0,
      age_label == "0.5_4" ~ 4,
      str_detect(age_label, "\\+$") ~ 107,
      TRUE ~ parse_number(
        str_extract(age_label, "[0-9]+$")))
  )

# takes a bit longer to run
age_risk_comp <- age_risk_long |>
  pivot_longer(
    age_min:age_max,
    names_to = "name",
    values_to = "value") |>
  mutate(
    name = "age",
    age_group = fct_inorder(age_group)) |>
  group_by(county_fips, age_group) |>
  complete(
    value = first(value):max(value)) |>
  ungroup()

# fill in
age_risk_comp_filled <- age_risk_comp |>
  group_by(county_fips, age_group) |>
  arrange(value, .by_group = TRUE) |>
  fill(
    prob_severe_inf,
    name,
    .direction = "up") |>
  ungroup()

# checking
age_risk_comp_filled |>
  filter(
    county_fips == "12119",
    value >= 60) |>
  select(
    age_group,
    value,
    prob_severe_inf) |>
  print(n = Inf)


# maps for imprinting distribution
load(here("data", "processed-data", "county_fine_age_pop.rda")) # read age distributions by county

county_age_dist <- county_fine_age_pop |>
  ungroup() |> 
  # add 0 probs for ages over the cutoff: OK because of sum and since all these ages imprinted same
  complete(
    nesting(birthyear=1924:1918, age=106:100), nesting(fips, county_name), 
    fill=list(pop=0, frac_pop=0)) |>
  arrange(fips, desc(birthyear))

county_risk_joined <- county_age_dist |>
  mutate(
    fips = sprintf("%05s", as.character(fips))
  ) |>
  filter(!str_starts(fips, "72")) |>
  group_by(fips) |>
  arrange(desc(birthyear), .by_group = TRUE) |>
  mutate(
    p_unprot_H5 = .env$p_unprot_H5,
    p_unprot_H7 = .env$p_unprot_H7
  ) |>
  ungroup() |>
  left_join(
    age_risk_comp_filled |>
      transmute(
        fips = sprintf("%05s", as.character(county_fips)),
        age = value,
        prob_severe_inf
      ),
    by = c("fips", "age")
  )


county_risk_for_plot <- county_risk_joined |>
  group_by(fips) |>
  summarise(
    risk_unprot_H5 = sum(frac_pop * p_unprot_H5),
    
    risk_severe_H5 = sum(
      frac_pop * p_unprot_H5 * prob_severe_inf),
    
    risk_comorb = sum(
      frac_pop * prob_severe_inf),
    
    risk_unprot_H7 = sum(
      frac_pop * p_unprot_H7),
    
    risk_severe_H7 = sum(
      frac_pop * p_unprot_H7 * prob_severe_inf),
    
    .groups = "drop")

# checking to make sure none missing
county_risk_joined |>
  summarise(
    n_counties = n_distinct(fips),
    n_missing_risk = sum(is.na(prob_severe_inf))
  )

# save data to use in plots
write_csv(county_risk_for_plot, here("data", "processed-data", "map_data_numbers.csv"))





#### Create dataset of ages and comorbs for later plotting #####################

# county level risk results we made above
here()
map_numbers <- read_csv(
  here("data", "processed-data", "map_data_numbers.csv"),
  show_col_types = FALSE) %>%
  mutate(fips = str_pad(as.character(fips), width = 5, pad = "0")
  )

# load age and comorbidity data and create a county lookup
comorbidities <- read_csv(here("data","processed-data","county_risk_weighted_with_dups.csv"),
  show_col_types = FALSE) %>%
  mutate(fips = str_pad(as.character(county_fips.x), width = 5, pad = "0")
  )

# make only one row per county so the join does not duplicate counties
county_lookup <- comorbidities %>%
  transmute(fips, county_name, state_fips = str_pad(as.character(state_fips), width = 2, pad = "0")) %>%
  distinct(fips, .keep_all = TRUE
  )


# median-age data
# Read in data
age2 <- read.csv(here("data", "processed-data", "ACSST5Y2023.S0101-Data.csv"))

# take off first row
age2 <- age2[-1, ]

# Clean and prepare median age
age2_median <- age2 %>%
  select(GEO_ID, NAME, median_age = S0101_C01_032E) %>%
  mutate(
    median_age = as.numeric(median_age),
    fips = stringr::str_pad(substr(GEO_ID, 10, 14), width = 5, pad = "0")
  )

age2_median <- age2_median %>%
  mutate(fips = str_pad(as.character(fips), width = 5, pad = "0")) %>%
  select(fips, median_age) %>%
  distinct(fips,.keep_all = TRUE
  )

# join
age_comorb_risksev <- map_numbers %>%
  left_join(county_lookup, by = "fips") %>%
  left_join(age2_median, by = "fips"
  )


# checking
names(age_comorb_risksev)

str(age_comorb_risksev$median_age)

age_comorb_risksev %>%
  summarise(
    rows = n(),
    counties = n_distinct(fips),
    missing_county_names = sum(is.na(county_name)),
    missing_median_age = sum(is.na(median_age)),
    missing_risk_comorb = sum(is.na(risk_comorb)),
    missing_risk_severe_H5 = sum(is.na(risk_severe_H5)),
    missing_risk_severe_H7 = sum(is.na(risk_severe_H7))
  )

age_comorb_risksev %>%
  select(
    fips,
    county_name,
    state_fips,
    median_age,
    risk_comorb,
    risk_unprot_H5,
    risk_severe_H5,
    risk_unprot_H7,
    risk_severe_H7
  ) %>%
  print(n = 20)


# save data to use in plots
write_csv(age_comorb_risksev, here("data", "processed-data", "age_comorb_risksev.csv"))


