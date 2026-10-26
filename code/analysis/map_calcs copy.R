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






p1 <- county_risk_plots |> 
  pivot_longer(contains("unprot")) |> 
  mutate(label = str_c(str_extract(name, "H\\d$"))) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80") +
  facet_wrap(~label, nrow = 1) +
  scale_fill_viridis_c(
    option = "inferno",
    trans = "log10",                      # 👈 apply log10 scale
    labels = scales::label_number(format = "e") # optional, nice formatting
  ) +
  labs(
    fill = "Proportion\nunprotected (log scale)",
    title = "Proportion without Imprinting Protection"
  ) +
  theme_map() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

p1

p2 <- county_risk_plots |> 
  pivot_longer(contains("severe")) |> 
  mutate(
    label = str_c(str_extract(name, "H\\d$")),
    value_log = log10(value + 1e-6)  #  add small constant to avoid log(0)
  ) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value_log), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80", alpha = 0.3) +
  facet_wrap(~label, ncol = 2) +
  scale_fill_viridis_c(
    option = "inferno",
    direction = 1,
    name = expression(log[10]~"(Proportion at risk)"),  # nice math label
    labels = scales::label_number(format = "e")          # show exponential-style numbers
  ) +
  labs(
    title = "Proportion at risk of severe illness (log10 scale)"
  ) +
  theme_map() +
  theme(strip.text = element_text(angle = 0))

p2

# save figure
ggsave(
  here("reconstructions-us-counties-main", "figs", "second_plots_logscale.png"),
  plot = p2,
  width = 8, height = 6, dpi = 300
)

library(dplyr)
library(ggplot2)
library(sf)
library(viridis)
library(patchwork)

# --- Prepare data ---
p2_data <- county_risk_plots |>
  pivot_longer(contains("severe")) |>
  mutate(
    label = str_c(str_extract(name, "H\\d$")),
    value_log = log10(value + 1e-6)  # avoid log(0)
  ) |>
  left_join(counties) |>
  sf::st_as_sf()  # ensure it's sf before plotting

# --- H5 plot ---
p2_H5 <- p2_data |>
  filter(label == "H5") |>
  ggplot() +
  geom_sf(aes(fill = value_log), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray80", alpha = 0.3) +
  scale_fill_viridis_c(
    option = "inferno",
    direction = 1,      # default
    name = "H5 (log10)"
  ) +
  labs(title = "H5 severe risk (log10 scale)") +
  theme_map()

# --- H7 plot (flipped) ---
p2_H7 <- p2_data |>
  filter(label == "H7") |>
  ggplot() +
  geom_sf(aes(fill = value_log), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray80", alpha = 0.3) +
  scale_fill_viridis_c(
    option = "inferno",
    direction = 1,    # flipped
    name = "H7 (log10)"
  ) +
  labs(title = "H7 severe risk (log10 scale)") +
  theme_map()

# Combine
library(patchwork)
p2 <- p2_H5 + p2_H7
p2






# Create the full dataset once
p2_data <- county_risk_plots |> 
  pivot_longer(contains("severe")) |> 
  mutate(label = str_c(str_extract(name, "H\\d$"), ": H=", rep(c(0.25, 0.24), n()/2))) |> 
  left_join(counties) |> 
  sf::st_as_sf()

#View(p2_data)
here()
saveRDS(p2_data, here("reconstructions-us-counties-main", "processed-data", "combined_map_data.rds"))

# Unique labels
labels <- unique(p2_data$label)

# First plot
p2_plot1 <- ggplot(p2_data |> filter(label == labels[1])) +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80", alpha = 0.3) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = NULL, title = paste("probability of severe infection -", labels[1])) +
  theme_map()

# Second plot
p2_plot2 <- ggplot(p2_data |> filter(label == labels[2])) +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80", alpha = 0.3) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = NULL, title = paste("probability of severe infection -", labels[2])) +
  theme_map()

# Display separately
p2_plot1

p2_plot2


ggsave(
  paste0("figs/HPAI-risk-map-p2_plot1-", today(), ".pdf"),
  plot = p2_plot1,
  width = 14.5,
  height = 11
)

ggsave(
  paste0("figs/HPAI-risk-map-p2_plot2-", today(), ".pdf"),
  plot = p2_plot2,
  width = 14.5,
  height = 11
)

## Find risk for a specific county
county_risk |> 
  filter(fips == "12119") |> 
  glimpse()


## Find number of counties that are at a risk of 0.30 or greate
county_risk |> 
  filter(risk_severe_H5 > 0.35) |> 
  nrow()/3222

county_risk |> 
  filter(risk_severe_H7 > 0.35) |> 
  nrow()/3222



# Get the mapping of counties in Florida
fl_counties <- us_map("counties") %>%
  filter(state == "FL") %>%
  select(fips, county)

# List of counties you want
target_counties <- c("Sumter", "Charlotte", "Sarasota", "Citrus", "Indian River", "Highlands", "Hernando")

# Filter for these counties and join with county_risk
county_risk %>%
  inner_join(fl_counties, by = "fips") %>%
  filter(county %in% target_counties) %>%
  select(county, fips, risk_severe_H5, risk_severe_H7)












# p_unprot <- plot_risk_map(county_risk, "risk_unprot", "unprotected", "inferno")
# p_severe <- plot_risk_map(county_risk, "risk_severe", "of severe infection", "viridis")
# plot_grid(p_unprot, p_severe, nrow=1)
# ggsave(paste0("figs/H5-risk-map-", today(), ".pdf"), width=10.5, height=6.5)

# Figure 2: association between county risk and social vulnerability index--------
library(ggpmisc) # for adding R^2 to regression plot

svi <- read_csv("../SVI_By_County.csv") |> 
  rename(fips=FIPS)

county_risk |> 
  left_join(svi) |> 
  ggplot(aes(SPL_THEMES, risk)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2"))) +
  geom_point(shape=1, alpha=0.7) +
  theme_half_open()

# Fig 3: compare the pred. case distribution between high and low risk counties---
hi_lo_counties <- filter(county_risk, risk == min(risk) | risk == max(risk))

hi_lo_risk <- county_age_dist |> 
  filter(fips %in% hi_lo_counties$fips) |> 
  mutate(rdist=frac_pop * prob_unprotected / sum(frac_pop * prob_unprotected)) # normalize to 1

ggplot(hi_lo_risk, aes(birthyear, rdist)) +
  geom_line() +
  facet_wrap(~county_name) +
  labs(y="predicted case distribution") +
  theme_half_open()

# the median age of these counties are respectively 23.9 vs. 68.3 according to Google
