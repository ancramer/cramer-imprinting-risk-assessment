# Load data
#   age_comorb_risksev = county-level risk_comorb values
age_comorb_risksev <- read.csv(here("data", "processed-data", "age_comorb_risksev.csv"))

#   county_avg_both_with_ages = county-level avg_* age-risk columns
#   this is county_risk_weighted_with_dups
county_avg_both_with_ages <- read.csv(here("data", "processed-data", "county_risk_weighted_with_dups.csv"))

# load packages
library(tidyverse)
library(here)

#### Functions #################################################################
# calculate a weighted mean after removing observations for which either the value or weight is missing
weighted_mean_complete <- function(x, w) {
  keep <- !is.na(x) & !is.na(w)
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w[keep])
}

# state references
state_lookup <- tibble(
  state_fips = c(
    "01", "02", "04", "05", "06", "08", "09", "10", "11", "12",
    "13", "15", "16", "17", "18", "19", "20", "21", "22", "23",
    "24", "25", "26", "27", "28", "29", "30", "31", "32", "33",
    "34", "35", "36", "37", "38", "39", "40", "41", "42", "44",
    "45", "46", "47", "48", "49", "50", "51", "53", "54", "55",
    "56"
  ),
  state_abbr = c(
    "AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL",
    "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME",
    "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH",
    "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI",
    "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI",
    "WY"
  )
)

# counties used in results
target_fips <- c(
  "12119", # Sumter County, FL
  "12015", # Charlotte County, FL
  "26001", # Alcona County, MI
  "48243", # Jeff Davis County, TX
  "48263", # Kent County, TX
  "49049", # Utah County, UT
  "49005", # Cache County, UT
  "48041", # Brazos County, TX
  "37133", # Onslow County, NC
  "13059"  # Clarke County, GA
)


#### Load data #################################################################
here()
map_data <- read_csv(here("data", "processed-data","map_data_numbers.csv"),
  show_col_types = FALSE) %>%
  mutate(fips = str_pad(as.character(fips), width = 5, pad = "0"))

county_pop_totals <- read_csv(here("data", "processed-data", "county_pop_totals.csv"),
  show_col_types = FALSE) %>%
  transmute(fips = str_pad(as.character(county_fips.x), width = 5, pad = "0"),
    county_pop = county_total)

map_data_pop <- map_data %>%
  left_join(county_pop_totals, by = "fips") %>%
  mutate(state_fips = str_sub(fips, 1, 2))


#### Vulnerability / lacking imprinting results ################################
#### state and county levels

state_weighted_vulnerability <- map_data_pop %>%
  group_by(state_fips) %>%
  summarise(
    H5_vulnerability = weighted_mean_complete(risk_unprot_H5, county_pop),
    H7_vulnerability = weighted_mean_complete(risk_unprot_H7, county_pop),
    .groups = "drop") %>%
  left_join(state_lookup, by = "state_fips") %>%
  select(state_abbr, state_fips, everything())

state_vulnerability_rank_H5 <- state_weighted_vulnerability %>%
  arrange(desc(H5_vulnerability))

state_vulnerability_rank_H7 <- state_weighted_vulnerability %>%
  arrange(desc(H7_vulnerability))

county_vulnerability_extremes <- list(
  H5_highest = map_data_pop %>%
    filter(!is.na(risk_unprot_H5)) %>%
    slice_max(risk_unprot_H5, n = 10, with_ties = FALSE) %>%
    select(fips, county_pop, risk_unprot_H5),
  
  H5_lowest = map_data_pop %>%
    filter(!is.na(risk_unprot_H5)) %>%
    slice_min(risk_unprot_H5, n = 10, with_ties = FALSE) %>%
    select(fips, county_pop, risk_unprot_H5),
  
  H7_highest = map_data_pop %>%
    filter(!is.na(risk_unprot_H7)) %>%
    slice_max(risk_unprot_H7, n = 10, with_ties = FALSE) %>%
    select(fips, county_pop, risk_unprot_H7),
  
  H7_lowest = map_data_pop %>%
    filter(!is.na(risk_unprot_H7)) %>%
    slice_min(risk_unprot_H7, n = 10, with_ties = FALSE) %>%
    select(fips, county_pop, risk_unprot_H7)
)


#### Severe disease risk #######################################################
#### lacking imprinting vulernarbility * (age effects + comorbs)

state_weighted_severe <- map_data_pop %>%
  group_by(state_fips) %>%
  summarise(
    H5_severe = weighted_mean_complete(risk_severe_H5, county_pop),
    H7_severe = weighted_mean_complete(risk_severe_H7, county_pop),
    .groups = "drop") %>%
  left_join(state_lookup, by = "state_fips") %>%
  select(state_abbr, state_fips, everything())

state_severe_rank_H5 <- state_weighted_severe %>%
  arrange(desc(H5_severe))

state_severe_rank_H7 <- state_weighted_severe %>%
  arrange(desc(H7_severe))

national_weighted_severe <- map_data_pop %>%
  summarise(
    H5_severe = weighted_mean_complete(risk_severe_H5, county_pop),
    H7_severe = weighted_mean_complete(risk_severe_H7, county_pop)
  )

# value for counties at a specific threshold
severe_threshold <- 0.35

county_pct_above_threshold <- map_data_pop %>%
  summarise(
    threshold = severe_threshold,
    H5_pct_counties = mean(
      risk_severe_H5 >= severe_threshold,
      na.rm = TRUE
    ) * 100,
    H7_pct_counties = mean(
      risk_severe_H7 >= severe_threshold,
      na.rm = TRUE
    ) * 100
  )

# unweighted proportion: each county contributes equally, regardless of size
prop_counties_H7_higher <- map_data_pop %>%
  summarise(
    proportion = mean(risk_severe_H7 > risk_severe_H5, na.rm = TRUE)
  )


# severe-risk threshold curve

threshold_results <- tibble(threshold = seq(0, 1, by = 0.01)) %>%
  mutate(
    H5_pct_counties = map_dbl(
      threshold,
      ~ mean(map_data_pop$risk_severe_H5 >= .x, na.rm = TRUE) * 100
    ),
    H7_pct_counties = map_dbl(
      threshold,
      ~ mean(map_data_pop$risk_severe_H7 >= .x, na.rm = TRUE) * 100
    )
  ) %>%
  pivot_longer(
    cols = ends_with("pct_counties"),
    names_to = "group",
    values_to = "percent_counties"
  ) %>%
  mutate(
    group = recode(
      group,
      H5_pct_counties = "H5",
      H7_pct_counties = "H7"
    )
  )

threshold_plot <- ggplot(
  threshold_results,
  aes(
    x = threshold * 100,
    y = percent_counties,
    color = group
  )
) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Severe-risk threshold (%)",
    y = "Counties at or above threshold (%)",
    color = "Group"
  ) +
  theme_minimal()


#### Comorbidity results #######################################################

# age_comorb_risksev must already exist before this section is run.
stopifnot(exists("age_comorb_risksev"))

comorbidity_data <- age_comorb_risksev %>%
  mutate(
    fips = str_pad(as.character(fips), width = 5, pad = "0"),
    state_fips = str_sub(fips, 1, 2)
  ) %>%
  left_join(county_pop_totals, by = "fips")

highest_comorbidity_county <- comorbidity_data %>%
  filter(!is.na(risk_comorb)) %>%
  slice_max(risk_comorb, n = 1, with_ties = TRUE)

lowest_comorbidity_county <- comorbidity_data %>%
  filter(!is.na(risk_comorb)) %>%
  slice_min(risk_comorb, n = 1, with_ties = TRUE)

state_weighted_comorbidity <- comorbidity_data %>%
  group_by(state_fips) %>%
  summarise(
    comorbidity_risk = weighted_mean_complete(risk_comorb, county_pop),
    .groups = "drop"
  ) %>%
  left_join(state_lookup, by = "state_fips") %>%
  select(state_abbr, state_fips, everything()) %>%
  arrange(desc(comorbidity_risk))


# age-risk data for age summaries

# county_avg_both_with_ages must still be in wide formay
stopifnot(exists("county_avg_both_with_ages"))

age_levels <- c(
  "avg_0_2", "avg_2_4", "avg_5_9", "avg_10_14", "avg_15_19",
  "avg_20_24", "avg_25_29", "avg_30_34", "avg_35_39", "avg_40_44",
  "avg_45_49", "avg_50_54", "avg_55_59", "avg_60_64",
  "avg_65_69", "avg_70_74", "avg_75+"
)

age_group_lookup <- tibble(
  age_group = age_levels,
  broad_age_group = c(
    "0-2", "2-4", "5-9", "10-14", "15-19",
    rep("20-44", 5),
    rep("45-64", 4),
    rep("65+", 3)
  )
)

county_age_long <- county_avg_both_with_ages %>%
  pivot_longer(
    cols = any_of(age_levels),
    names_to = "age_group",
    values_to = "risk"
  ) %>%
  mutate(
    age_group = factor(age_group, levels = age_levels),
    fips = str_pad(as.character(county_fips.x), width = 5, pad = "0")
  ) %>%
  left_join(age_group_lookup, by = "age_group") %>%
  left_join(county_pop_totals, by = "fips")

broad_age_county_weighted_risk <- county_age_long %>%
  group_by(broad_age_group) %>%
  summarise(
    risk = weighted_mean_complete(risk, county_pop),
    .groups = "drop"
  ) %>%
  mutate(
    broad_age_group = factor(
      broad_age_group,
      levels = c(
        "0-2", "2-4", "5-9", "10-14", "15-19",
        "20-44", "45-64", "65+"
      )
    )
  ) %>%
  arrange(broad_age_group)

# values for those counties

target_county_results <- map_data %>%
  filter(fips %in% target_fips) %>%
  select(
    fips,
    risk_unprot_H5,
    risk_unprot_H7,
    risk_severe_H5,
    risk_severe_H7
  ) %>%
  arrange(match(fips, target_fips))

# main result objects

state_vulnerability_rank_H5
state_vulnerability_rank_H7
county_vulnerability_extremes
state_severe_rank_H5
state_severe_rank_H7
national_weighted_severe
county_pct_above_threshold
prop_counties_H7_higher
threshold_plot
highest_comorbidity_county
lowest_comorbidity_county
state_weighted_comorbidity
broad_age_county_weighted_risk
target_county_results
