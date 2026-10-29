# load in packages
library(tidyverse)
library(usmap)
library(cowplot)
library(here)

# load in data
county_risk_plots <- read.csv(here("data", "processed-data", "map_data_numbers.csv"))

# county shapefiles for plotting
counties <- us_map(regions="counties")
states <- us_map(regions="states")

# Standardize fips for joins
county_risk_plots <- county_risk_plots |>
  mutate(fips = as.character(fips))

counties <- counties |>
  mutate(fips = as.character(fips))

states <- states |>
  mutate(fips = as.character(fips))

# imprinting only maps
p1 <- county_risk_plots |> 
  pivot_longer(contains("unprot")) |> 
  mutate(
    fips = stringr::str_pad(as.character(fips), width = 5, pad = "0"),
    label = stringr::str_extract(name, "H\\d$")
  ) |> 
  left_join(
    counties |> 
      mutate(fips = stringr::str_pad(as.character(fips), width = 5, pad = "0")),
    by = "fips"
  ) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80", alpha = 0.3) +
  facet_wrap(~label, nrow = 2) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = NULL, title = "Proportion without Imprinting Protection") +
  theme_map()
p1


severe_limits <- range(
  county_risk_plots$risk_severe_H5,
  county_risk_plots$risk_severe_H7,
  na.rm = TRUE
)

p2 <- county_risk_plots |> 
  pivot_longer(contains("severe")) |> 
  mutate(
    fips = stringr::str_pad(as.character(fips), width = 5, pad = "0"),
    label = stringr::str_extract(name, "H\\d$")
  ) |> 
  left_join(
    counties |> 
      mutate(fips = stringr::str_pad(as.character(fips), width = 5, pad = "0")),
    by = "fips"
  ) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80", alpha = 0.3) +
  facet_wrap(~label, nrow = 2) +
  scale_fill_viridis_c(
    option = "inferno",
    limits = severe_limits,
    oob = scales::squish
  ) +
  labs(
    fill = NULL,
    title = "Limited imprinting protection and severe-illness risk factors"
  ) +
  theme_map()

p2


source(file.path(project_dir, "R", "calculation_funs.R"))
source(file.path(project_dir, "R", "data_import_funs.R"))


#===============================================================================
# Calculate imprinting probabilities
#===============================================================================

get_p_unprot <- function(Wmat, H, group) {
  if (group == 1) {
    p_imprint <- Wmat$H1N1_probs + Wmat$H2N2_probs
  } else if (group == 2) {
    p_imprint <- Wmat$H3N2_probs
  } else {
    stop("Please select group 1 or group 2.")
  }
  
  p_imprint * H + (1 - p_imprint)
}

# This function expects to run from the project's R folder
original_wd <- getwd()
setwd(file.path(project_dir, "R"))

Wmat <- get_imprinting_probabilities(
  observation_years = 2024
)

setwd(original_wd)

p_unprot_H5 <- get_p_unprot(
  Wmat,
  H = 0.25,
  group = 1
)

p_unprot_H7 <- get_p_unprot(
  Wmat,
  H = 0.24,
  group = 2
)


#===============================================================================
# Load county comorbidity/age risk
#===============================================================================

age_risk <- read_csv(
  file.path(
    project_dir,
    "processed-data",
    "county_risk_weighted_with_dups.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(county_fips = county_fips.x)


#===============================================================================
# Convert age groups into individual ages
#===============================================================================

age_risk_comp <- age_risk %>%
  pivot_longer(
    cols = starts_with("avg_"),
    names_to = "age_group",
    values_to = "prob_severe_inf"
  ) %>%
  mutate(
    age_min = parse_number(
      str_extract(age_group, "(?<=_)[\\d]+")
    ),
    age_max = parse_number(
      str_extract(age_group, "(?<=_)[\\d]+$")
    )
  ) %>%
  pivot_longer(
    cols = age_min:age_max,
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    name = "age",
    age_group = fct_inorder(age_group),
    value = case_match(
      value,
      0.5 ~ 1,
      NA ~ 107,
      .default = value
    )
  ) %>%
  group_by(county_fips, age_group) %>%
  complete(value = first(value):max(value)) %>%
  ungroup() %>%
  fill(
    prob_severe_inf,
    name,
    .direction = "up"
  ) %>%
  select(
    fips = county_fips,
    age = value,
    prob_severe_inf
  ) %>%
  distinct(
    fips,
    age,
    .keep_all = TRUE
  )


#===============================================================================
# Load county population distributions
#===============================================================================

load(
  file.path(
    project_dir,
    "processed-data",
    "county_fine_age_pop.rda"
  )
)

county_age_dist <- county_fine_age_pop %>%
  ungroup() %>%
  complete(
    nesting(
      birthyear = 1924:1918,
      age = 106:100
    ),
    nesting(fips, county_name),
    fill = list(
      pop = 0,
      frac_pop = 0
    )
  ) %>%
  arrange(fips, desc(birthyear))


#===============================================================================
# Calculate severe-illness risk by county
#===============================================================================

county_risk <- county_age_dist %>%
  left_join(
    age_risk_comp,
    by = c("fips", "age")
  ) %>%
  group_by(fips) %>%
  summarise(
    risk_severe_H5 = sum(
      frac_pop * p_unprot_H5 * prob_severe_inf,
      na.rm = TRUE
    ),
    risk_severe_H7 = sum(
      frac_pop * p_unprot_H7 * prob_severe_inf,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


#===============================================================================
# Load national map boundaries
#===============================================================================

counties <- us_map(regions = "counties")
states   <- us_map(regions = "states")


#===============================================================================
# Create national p2 — original scale, not state-specific
#===============================================================================

p2 <- county_risk %>%
  pivot_longer(
    cols = contains("severe"),
    names_to = "risk_type",
    values_to = "value"
  ) %>%
  mutate(
    label = str_extract(risk_type, "H\\d$"),
    label = recode(
      label,
      "H5" = "Group 1",
      "H7" = "Group 2"
    )
  ) %>%
  left_join(counties, by = "fips") %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(
    aes(fill = value),
    color = NA
  ) +
  geom_sf(
    data = states,
    fill = NA,
    color = "gray80"
  ) +
  facet_wrap(
    ~ label,
    nrow = 1
  ) +
  scale_fill_viridis_c(
    option = "inferno"
  ) +
  labs(
    fill = "Proportion\nat risk"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    strip.background = element_blank(),
    strip.text = element_text(
      margin = margin(b = 1, t = 1)
    ),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

p2



