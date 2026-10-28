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
