#### Packages #################################################################
library(tidyverse)
library(here)
library(usmap)
library(cowplot)
library(viridis)
library(sf)
library(ggnewscale)

source(here("code", "processing", "calculation_funs.R"))
source(here("code", "processing", "data_import_funs.R"))

#### Figure 1 #################################################################
# Further edits made by Daniel Humphrey, graphic designer at the Odum School
# of Ecology, University of Georgia, Athens, GA.
# Contact atc10260@uga.edu with questions.

# takes a min to run
US <- get_imprinting_probabilities(observation_years = 2010:2024)

unitedstates <- bind_rows(
  as.data.frame(US$H1N1_probs) |> mutate(obs_year = rownames(US$H1N1_probs), strain = "H1N1"),
  as.data.frame(US$H2N2_probs) |> mutate(obs_year = rownames(US$H2N2_probs), strain = "H2N2"),
  as.data.frame(US$H3N2_probs) |> mutate(obs_year = rownames(US$H3N2_probs), strain = "H3N2"),
  as.data.frame(US$naive_probs) |> mutate(obs_year = rownames(US$naive_probs), strain = "naive")
)

year_cols <- grep("^[0-9]{4}$", names(unitedstates), value = TRUE)

unitedstates_long <- unitedstates |>
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "birth_year",
    values_to = "probability"
  ) |>
  mutate(
    birth_year = as.numeric(birth_year),
    obs_year = as.numeric(obs_year),
    strain = factor(strain, levels = c("naive", "H1N1", "H2N2", "H3N2"))
  )

unitedstates_plot <- unitedstates_long |>
  filter(obs_year == 2024) |>
  ggplot(aes(birth_year, probability, fill = strain)) +
  geom_col(width = 1) +
  scale_fill_manual(
    values = c(
      "H1N1" = "#e08214",
      "H2N2" = "#fdb863",
      "H3N2" = "#8073ac",
      "naive" = "grey80"
    )
  ) +
  scale_x_continuous(
    breaks = c(1918, 1957, 1968, 1977, 2009, 2024),
    minor_breaks = NULL
  ) +
  labs(
    x = "Birth Year",
    y = "Proportion of Population",
    title = "Imprinting Probabilities in the US",
    fill = "Strain"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

unitedstates_plot

ggsave(
  here("code", "figs", "US_birth_cohorts.png"),
  unitedstates_plot,
  width = 8,
  height = 4,
  dpi = 300
)

################################################################################
#### Shared formatting for panel Figures 2, 3, and 4 ###########################
################################################################################

group_colors <- c("Group 1" = "#e08214", "Group 2" = "#8073ac")
line_colors <- c("Group 1" = "#ad6511", "Group 2" = "#574d78")
panel_label_size <- 14

single_map_theme <- theme(
  panel.spacing = unit(0, "lines"),
  plot.margin = margin(0, 0, 0, 0),
  strip.background = element_blank(),
  strip.text = element_blank(),
  legend.position = "right",
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0)
)

faceted_map_theme <- theme(
  panel.spacing = unit(0, "lines"),
  strip.background = element_blank(),
  strip.text = element_text(margin = margin(b = 1, t = 1)),
  legend.position = "right"
)

line_plot_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
)

violin_axis_theme <- theme(
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title = element_text(size = 14)
)

figure4_theme <- theme(
  plot.title = element_blank(),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.position = "bottom",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.6, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0)
)

add_panel_title <- function(plot, title, height = 0.95, vjust = 1.5,
                            lineheight = 1.1) {
  ggdraw() +
    draw_label(
      title,
      x = 0.5,
      y = 1,
      hjust = 0.5,
      vjust = vjust,
      size = 18,
      lineheight = lineheight
    ) +
    draw_plot(plot, y = 0, height = height)
}

assemble_three_panel <- function(a, b, c, titles) {
  a_title <- add_panel_title(a, titles[[1]])
  b_title <- add_panel_title(b, titles[[2]], height = 0.85, vjust = 1.2)
  c_title <- add_panel_title(c, titles[[3]])
  
  top_row <- plot_grid(
    a_title,
    b_title,
    ncol = 2,
    rel_widths = c(1.5, 1.2),
    labels = c("A", "B"),
    label_size = 16,
    label_fontface = "bold",
    align = "v"
  )
  
  bottom_row <- plot_grid(
    c_title,
    ncol = 1,
    labels = "C",
    label_size = 16,
    label_fontface = "bold"
  )
  
  plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
}

pad_fips <- function(x) str_pad(as.character(x), width = 5, pad = "0")

r2_label <- function(model, digits = 3) {
  paste0("R² = ", formatC(summary(model)$r.squared, digits = digits, format = "f"))
}

counties <- us_map(regions = "counties") |> mutate(fips = pad_fips(fips))
states <- us_map(regions = "states")

#### Panel Figure 2 ############################################################
#### A: Median-age map #########################################################

age2 <- read.csv(here("data", "processed-data", "ACSST5Y2023.S0101-Data.csv"))[-1, ]

age2_median <- age2 |>
  select(GEO_ID, NAME, median_age = S0101_C01_032E) |>
  mutate(
    median_age = as.numeric(median_age),
    fips = pad_fips(substr(GEO_ID, 10, 14))
  )

plot_data <- left_join(counties, age2_median, by = "fips")

us_states <- counties |>
  group_by(abbr) |>
  summarise(geometry = st_union(geom), .groups = "drop")

p_age_formatted <- plot_data |>
  mutate(dummy = "Map A") |>
  ggplot() +
  geom_sf(aes(fill = median_age), color = NA) +
  geom_sf(data = counties, fill = NA, color = NA, alpha = 0.6) +
  geom_sf(data = us_states, fill = NA, color = "gray80", linewidth = 0.25) +
  facet_wrap(~dummy, nrow = 1) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "gray90",
    name = "Median \nage"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  single_map_theme

p_age_formatted

#### B: Vulnerability by birth year ############################################

get_p_unprot <- function(Wmat, H, group) {
  p_imprint <- switch(
    as.character(group),
    "1" = Wmat$H1N1_probs + Wmat$H2N2_probs,
    "2" = Wmat$H3N2_probs,
    stop("Please select group 1 or 2.")
  )
  p_imprint * H + (1 - p_imprint)
}

Wmat <- get_imprinting_probabilities(observation_years = 2024)
p_unprot_H5 <- get_p_unprot(Wmat, H = 0.25, group = 1)
p_unprot_H7 <- get_p_unprot(Wmat, H = 0.24, group = 2)

plot_df <- tibble(
  birthyear = as.integer(names(p_unprot_H5)),
  unprot_H5 = as.numeric(p_unprot_H5),
  unprot_H7 = as.numeric(p_unprot_H7)
) |>
  pivot_longer(starts_with("unprot_"), names_to = "strain", values_to = "prop_unprotected") |>
  mutate(strain = recode(strain, unprot_H5 = "Group 1", unprot_H7 = "Group 2")) |>
  arrange(birthyear)

p_vul <- ggplot(plot_df, aes(birthyear, prop_unprotected, color = strain)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    limits = c(0.2, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  scale_x_continuous(
    limits = c(1918, 2024),
    breaks = c(1918, 1957, 1968, 1977, 2009, 2015, 2024),
    minor_breaks = seq(1918, 2024, 5)
  ) +
  labs(
    x = "Birth year",
    y = "Proportion vulnerable (unprotected)",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  line_plot_theme

p_vul

#### C: Vulnerability maps #####################################################

county_risk_plots <- read.csv(
  here("data", "processed-data", "map_data_numbers.csv")
) |>
  mutate(fips = pad_fips(fips))

p1_data <- county_risk_plots |>
  pivot_longer(contains("unprot"), names_to = "risk_type", values_to = "value") |>
  mutate(
    label = str_extract(risk_type, "H\\d$"),
    label = recode(label, H5 = "Group 1", H7 = "Group 2")
  ) |>
  left_join(counties, by = "fips") |>
  st_as_sf()

p1 <- ggplot(p1_data) +
  geom_sf(aes(fill = value), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray80", alpha = 0.3) +
  facet_wrap(~label, nrow = 1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "Proportion\nvulnerable") +
  coord_sf(expand = FALSE) +
  theme_map() +
  faceted_map_theme

p1

figure2_panel <- assemble_three_panel(
  p_age_formatted,
  p_vul,
  p1,
  c(
    "Age distribution by county across the United States",
    "Proportion unprotected by imprint by birth year",
    "Proportion of population unprotected by imprint in the United States"
  )
)

figure2_panel

ggsave(
  here("code", "figs", "fig_2.png"),
  figure2_panel,
  width = 13,
  height = 8,
  dpi = 300
)

#### Panel Figure 3 ############################################################
#### A: Comorbidity map ########################################################

age_comorb_risksev <- read.csv(
  here("data", "processed-data", "age_comorb_risksev.csv")
) |>
  mutate(fips = pad_fips(fips))

county_sf <- counties |>
  left_join(age_comorb_risksev |> select(fips, risk_comorb), by = "fips") |>
  st_as_sf(sf_column_name = "geom")

p_comorbs <- county_sf |>
  mutate(dummy = "Map A") |>
  ggplot() +
  geom_sf(aes(fill = risk_comorb), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray80", linewidth = 0.3) +
  facet_wrap(~dummy, nrow = 1) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "gray90",
    name = "Comorbidity\nrate"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  single_map_theme

p_comorbs

#### B: Risk-factor violin plot ################################################

county_avg_both2 <- read_csv(
  here ("data", "processed-data", "county_avg_both2.csv"),
  show_col_types = FALSE
)

age_levels <- c(
  "avg_0_2", "avg_2_4", "avg_5_9", "avg_10_14", "avg_15_19",
  "avg_20_24", "avg_25_29", "avg_30_34", "avg_35_39", "avg_40_44",
  "avg_45_49", "avg_50_54", "avg_55_59", "avg_60_64", "avg_65_69",
  "avg_70_74", "avg_75+"
)

county_avg_both_with_ages <- county_avg_both2 |>
  select(-starts_with("avg_X")) |>
  mutate(
    avg_0_2 = 1,
    avg_2_4 = avg_0.5_4,
    avg_65_69 = 1,
    avg_70_74 = 1,
    `avg_75+` = 1
  ) |>
  select(-avg_0_0.5, -avg_0.5_4) |>
  pivot_longer(starts_with("avg_"), names_to = "age_group", values_to = "risk") |>
  mutate(age_group = factor(age_group, levels = age_levels, ordered = TRUE))

q <- ggplot(county_avg_both_with_ages, aes(age_group, risk)) +
  geom_jitter(width = 0.15, alpha = 0.3, color = "#a6317d", size = 1.5) +
  geom_violin(fill = "#fcfdbf", color = "black") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(
    labels = function(x) str_replace_all(str_remove(x, "^avg_"), "_", "-")
  ) +
  labs(x = "Age group", y = "Proportion of population\nwith a risk factor") +
  theme_minimal(base_size = 14) +
  violin_axis_theme

q

#### C: Severe-influenza-risk maps ############################################

p2_data <- county_risk_plots |>
  pivot_longer(contains("severe"), names_to = "risk_type", values_to = "value") |>
  mutate(
    label = str_extract(risk_type, "H\\d$"),
    label = recode(label, H5 = "Group 1", H7 = "Group 2")
  ) |>
  left_join(counties, by = "fips") |>
  st_as_sf()

p2 <- ggplot(p2_data) +
  geom_sf(aes(fill = value), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray80", alpha = 0.3) +
  facet_wrap(~label, nrow = 1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "Proportion\nat risk") +
  coord_sf(expand = FALSE) +
  theme_map() +
  faceted_map_theme

p2

figure3_panel <- assemble_three_panel(
  p_comorbs,
  q,
  p2,
  c(
    "Comorbidity distribution by county across the United States",
    "Prevalence of demographic vulnerability factors\n(age and comorbidities)",
    "Proportion of population at high risk of severe influenza in the United States"
  )
)

figure3_panel

ggsave(
  here("code", "figs", "fig_3.png"),
  figure3_panel,
  width = 13,
  height = 8,
  dpi = 300
)

#### Panel Figure 4 ############################################################
#### A-C: County-level relationships ###########################################

age_comorb_risksev$median_age <- as.numeric(age_comorb_risksev$median_age)

lm_age_H5 <- lm(risk_severe_H5 ~ median_age, data = age_comorb_risksev)
lm_age_H7 <- lm(risk_severe_H7 ~ median_age, data = age_comorb_risksev)

plot_age2 <- ggplot(age_comorb_risksev) +
  geom_point(aes(median_age, risk_severe_H7, color = "Group 2"), alpha = 0.15, size = 2) +
  geom_point(aes(median_age, risk_severe_H5, color = "Group 1"), alpha = 0.15, size = 2) +
  geom_abline(
    slope = coef(lm_age_H5)[["median_age"]],
    intercept = coef(lm_age_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]], linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_age_H7)[["median_age"]],
    intercept = coef(lm_age_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]], linewidth = 1.2
  ) +
  annotate("text", x = 50, y = 0.60, label = r2_label(lm_age_H5),
           color = line_colors[["Group 1"]], hjust = 1, size = 7) +
  annotate("text", x = 50, y = 0.70, label = r2_label(lm_age_H7),
           color = line_colors[["Group 2"]], hjust = 1, size = 7) +
  scale_color_manual(values = group_colors, name = "Group") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  labs(x = "Median age", y = "Proportion at high risk\nfor severe influenza") +
  theme_minimal(base_size = 18) +
  figure4_theme

lm_risk_H5 <- lm(risk_severe_H5 ~ risk_comorb, data = age_comorb_risksev)
lm_risk_H7 <- lm(risk_severe_H7 ~ risk_comorb, data = age_comorb_risksev)

plot_risk2 <- ggplot(age_comorb_risksev) +
  geom_point(aes(risk_comorb, risk_severe_H7, color = "Group 2"), alpha = 0.15, size = 2) +
  geom_point(aes(risk_comorb, risk_severe_H5, color = "Group 1"), alpha = 0.15, size = 2) +
  geom_abline(
    slope = coef(lm_risk_H5)[["risk_comorb"]],
    intercept = coef(lm_risk_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]], linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_risk_H7)[["risk_comorb"]],
    intercept = coef(lm_risk_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]], linewidth = 1.2
  ) +
  annotate("text", x = 0.485, y = 0.60, label = r2_label(lm_risk_H5),
           color = line_colors[["Group 1"]], hjust = 1, size = 7) +
  annotate("text", x = 0.50, y = 0.70, label = r2_label(lm_risk_H7),
           color = line_colors[["Group 2"]], hjust = 1, size = 7) +
  scale_color_manual(values = group_colors, name = "Group") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  labs(x = "Proportion with comorbidity", y = "Proportion at high risk\nfor severe influenza") +
  coord_cartesian(xlim = c(0.2, 0.7)) +
  theme_minimal(base_size = 18) +
  figure4_theme

lm_unprotected_H5 <- lm(risk_severe_H5 ~ risk_unprot_H5, data = age_comorb_risksev)
lm_unprotected_H7 <- lm(risk_severe_H7 ~ risk_unprot_H7, data = age_comorb_risksev)

x_H5 <- quantile(age_comorb_risksev$risk_unprot_H5, 0.01, na.rm = TRUE)
x_H7 <- quantile(age_comorb_risksev$risk_unprot_H7, 0.01, na.rm = TRUE)
y_H5 <- quantile(age_comorb_risksev$risk_severe_H5, 0.97, na.rm = TRUE)
y_H7 <- quantile(age_comorb_risksev$risk_severe_H7, 0.97, na.rm = TRUE)

plot_unimprinted <- ggplot(age_comorb_risksev) +
  geom_point(aes(risk_unprot_H7, risk_severe_H7, color = median_age), alpha = 0.6, size = 2) +
  scale_color_gradient(low = "#d8daeb", high = "#2d004b", name = "Median age") +
  ggnewscale::new_scale_color() +
  geom_point(aes(risk_unprot_H5, risk_severe_H5, color = median_age), alpha = 0.6, size = 2) +
  scale_color_gradient(low = "#fee0b6", high = "#7f3b08", name = "Median age") +
  geom_abline(
    slope = coef(lm_unprotected_H7)[["risk_unprot_H7"]],
    intercept = coef(lm_unprotected_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]], linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_unprotected_H5)[["risk_unprot_H5"]],
    intercept = coef(lm_unprotected_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]], linewidth = 1.2
  ) +
  annotate("text", x = x_H5, y = y_H5, label = r2_label(lm_unprotected_H5),
           color = line_colors[["Group 1"]], hjust = 0, size = 7) +
  annotate("text", x = x_H7, y = y_H7, label = r2_label(lm_unprotected_H7),
           color = line_colors[["Group 2"]], hjust = 0, size = 7) +
  labs(x = "Proportion unimprinted/unprotected",
       y = "Proportion at high risk\nfor severe influenza") +
  theme_minimal(base_size = 18) +
  figure4_theme

#### D-E: State and census-region results ######################################

state_lookup_df <- tibble(
  state_fips = c("01", "02", "04", "05", "06", "08", "09", "10", "12", "13", "15", "16",
                 "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28",
                 "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                 "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
                 "55", "56"),
  state = state.name
)

map_numbers <- county_risk_plots |>
  mutate(state_fips = substr(fips, 1, 2)) |>
  filter(!state_fips %in% c("11", "72")) |>
  left_join(state_lookup_df, by = "state_fips") |>
  left_join(tibble(state = state.name, state_abb = state.abb), by = "state") |>
  transmute(
    across(everything()),
    state = state_abb
  ) |>
  select(-state_abb) |>
  mutate(
    census_region_main = case_when(
      state %in% c("ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA") ~ "Northeast",
      state %in% c("OH", "IN", "IL", "MI", "WI", "MN", "IA", "MO", "ND", "SD", "NE", "KS") ~ "Midwest",
      state %in% c("DE", "MD", "VA", "WV", "NC", "SC", "GA", "FL", "KY", "TN", "MS", "AL", "AR", "LA", "OK", "TX") ~ "South",
      state %in% c("MT", "ID", "WY", "CO", "NM", "AZ", "UT", "NV", "WA", "OR", "CA", "AK", "HI") ~ "West",
      TRUE ~ NA_character_
    )
  )

county_pop_totals <- read.csv(
  here("data", "processed-data", "county_pop_totals.csv")
) |>
  mutate(county_fips.x = sprintf("%05d", as.integer(county_fips.x)))

data_for_regions_states <- map_numbers |>
  left_join(county_pop_totals, by = c("fips" = "county_fips.x")) |>
  rename(county_pop = county_total)

state_risks <- data_for_regions_states |>
  group_by(state, state_fips, census_region_main) |>
  summarise(
    across(starts_with("risk"), ~ weighted.mean(.x, county_pop, na.rm = TRUE)),
    total_pop = sum(county_pop, na.rm = TRUE),
    .groups = "drop"
  )

region_colors <- c(
  "Northeast" = "#99DDFFFF",
  "West" = "#DCB0F2FF",
  "South" = "#F89C74FF",
  "Midwest" = "#C9DB74FF"
)

make_state_plot <- function(risk_var, group_number) {
  ordered_states <- state_risks |>
    mutate(state = fct_reorder(state, .data[[risk_var]]))
  
  county_data <- data_for_regions_states |>
    mutate(state = factor(state, levels = levels(ordered_states$state)))
  
  ggplot() +
    geom_violin(
      data = county_data,
      aes(y = state, x = .data[[risk_var]], fill = census_region_main),
      color = NA, alpha = 0.9, scale = "width"
    ) +
    geom_jitter(
      data = county_data,
      aes(y = state, x = .data[[risk_var]]),
      color = "grey40", alpha = 0.3, height = 0.25, size = 1
    ) +
    geom_point(
      data = ordered_states,
      aes(y = state, x = .data[[risk_var]]),
      color = "grey30", alpha = 0.8, size = 2
    ) +
    scale_fill_manual(values = region_colors) +
    labs(
      y = "State",
      x = paste0("Proportion at high risk \nfor severe group ", group_number, " disease"),
      fill = "U.S. Census Region"
    ) +
    theme_minimal(base_size = 12) +
    figure4_theme +
    theme(axis.text.y = element_text(size = 7))
}

G1_results <- make_state_plot("risk_severe_H5", 1)
G2_results <- make_state_plot("risk_severe_H7", 2)

#### Assemble Figure 4 #########################################################

a <- plot_age2 + labs(y = NULL)
b <- plot_risk2 + labs(y = NULL)
c <- plot_unimprinted + labs(y = NULL)
d <- G1_results
e <- G2_results

legend_groups <- get_legend(a)
legend_age <- get_legend(c)
legend_bottom <- get_legend(d)

a <- a + theme(legend.position = "none")
b <- b + theme(legend.position = "none")
c <- c + theme(legend.position = "none")
d <- d + theme(legend.position = "none")
e <- e + theme(legend.position = "none")

shared_y_axis <- ggdraw() +
  draw_label("Proportion at high risk\nfor severe influenza", angle = 90, size = 14)

top_panels <- plot_grid(
  a, b, c,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = panel_label_size,
  label_fontface = "bold",
  label_x = -0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "hv",
  axis = "tblr"
)

top_main <- plot_grid(shared_y_axis, top_panels, ncol = 2, rel_widths = c(0.06, 1))

top_legend_row <- plot_grid(
  legend_groups,
  NULL,
  plot_grid(NULL, legend_age, ncol = 2, rel_widths = c(0.15, 1)),
  ncol = 3,
  rel_widths = c(1, 0.2, 1)
)

top_row <- plot_grid(top_main, top_legend_row, ncol = 1, rel_heights = c(1, 0.14))

bottom_main <- plot_grid(
  d, e,
  ncol = 2,
  labels = c("D", "E"),
  label_size = panel_label_size,
  label_fontface = "bold",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "h"
)

bottom_row <- plot_grid(
  bottom_main,
  legend_bottom,
  ncol = 1,
  rel_heights = c(1, 0.06)
)

figure4_panel <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.78, 1.9)
)

figure4_panel

ggsave(
  here("code", "figs", "fig_4.png"),
  figure4_panel,
  width = 12,
  height = 10,
  dpi = 300,
  limitsize = FALSE
)

