# load in packages
library(tidyverse)
library(usmap)
library(cowplot)
library(here)

# load in functions 
here()
source(here("code", "processing", "calculation_funs.R"))
source(here("code", "processing", "data_import_funs.R"))

# build imrpinting HA groups: 1 is H1N1 and H2N2, 2 is H3N2
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

#here::here()
#setwd("/Users/annalisecramer/Documents/H5 OAS Project/reconstructions-us-counties-main/R")
#getwd()
#Wmat <- get_imprinting_probabilities(observation_years=2024)
#setwd("..")
#getwd()

# add in the failure rate parameter
p_get_anyway <-  # the H parameter
  
p_unprot_H5 <- get_p_unprot(Wmat, H=0.25, group=1)
p_unprot_H7 <- get_p_unprot(Wmat, H=0.24, group=2)

# load in the county risk (age and comorbs) data
age_risk <- read_csv(here("data", "processed-data", "county_risk_weighted_with_dups.csv")) %>%
  rename(county_fips = county_fips.x)


age_risk_long <- age_risk |>
  pivot_longer(
    cols = starts_with("avg_"),
    names_to = "age_group",
    values_to = "prob_severe_inf"
  ) |>
  mutate(
    age_min = parse_number(str_extract(age_group, "(?<=_)[\\d]+")),
    age_max = parse_number(str_extract(age_group, "(?<=_)[\\d]+$"))
  )

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

age_risk_comp <- age_risk_comp |> 
  fill(prob_severe_inf, name, .direction="up") |> 
  select(fips=county_fips, age=value, prob_severe_inf) |> 
  distinct(fips, age, .keep_all=TRUE) # removes the few counties with two entries per fips

# 10/16/26 pick up here!!! 
# Figure 1: county level risk-----------------------------------------------------
load("processed-data/county_fine_age_pop.rda") # read Spencer's age distributions by county
glimpse(county_fine_age_pop)

county_age_dist <- county_fine_age_pop |>
  ungroup() |> 
  # add 0 probs for ages over the cutoff: OK because of sum and since all these ages imprinted same
  complete(
    nesting(birthyear=1924:1918, age=106:100), nesting(fips, county_name), 
    fill=list(pop=0, frac_pop=0)
  ) |>
  arrange(fips, desc(birthyear))

county_risk <- county_age_dist |> 
  left_join(age_risk_comp, by=c("fips", "age")) |> 
  # select(-county_name.y, county_name=county_name.x) |> 
  # filter(is.na(prob_severe_inf)) |>
  # View()
  group_by(fips) |> 
  # weighted average of prob_unprotected over the age distribution
  summarise(
    risk_unprot_H5=sum(frac_pop * p_unprot_H5),
    risk_severe_H5=sum(frac_pop * p_unprot_H5 * prob_severe_inf),
    risk_unprot_H7=sum(frac_pop * p_unprot_H7),
    risk_severe_H7=sum(frac_pop * p_unprot_H7 * prob_severe_inf),
  )


county_risk_for_plot <- county_age_dist |> 
  left_join(age_risk_comp, by=c("fips", "age")) |> 
  # select(-county_name.y, county_name=county_name.x) |> 
  # filter(is.na(prob_severe_inf)) |>
  # View()
  group_by(fips) |> 
  # weighted average of prob_unprotected over the age distribution
  summarise(
    risk_unprot_H5=sum(frac_pop * p_unprot_H5),
    risk_severe_H5=sum(frac_pop * p_unprot_H5 * prob_severe_inf),
    risk_comorb=sum(frac_pop * prob_severe_inf),
    risk_unprot_H7=sum(frac_pop * p_unprot_H7),
    risk_severe_H7=sum(frac_pop * p_unprot_H7 * prob_severe_inf),
  )

glimpse(county_risk_for_plot)
str(county_risk_for_plot)
#View(county_age_dist)
fl <- county_risk_for_plot %>%
  filter(fips == "12119")
str(fl)

here()
write_csv(county_risk_for_plot, here("reconstructions-us-counties-main", "processed-data", "map_data_numbers.csv"))

counties <- us_map(regions="counties") # county shapefiles for plotting
states <- us_map(regions="states")

p1 <- county_risk |> 
  pivot_longer(contains("unprot")) |> 
  mutate(label=str_c(str_extract(name, "H\\d$"))) |> #, ": H=", rep(c(0.25, 0.24), n()/2))) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill=value), col=NA) +
  geom_sf(data=states, fill=NA, col="gray80", alpha=0.3) +
  facet_wrap(~label, nrow=2) +
  scale_fill_viridis_c(option="inferno") +
  labs(fill=NULL, title="Proportion without Imprinting Protection") +
  # labs(fill=NULL, title=paste0("probability unprotected with H=", p_get_anyway)) +
  theme_map()
p1

p2 <- county_risk |> 
  pivot_longer(contains("severe")) |> 
  mutate(label=str_c(str_extract(name, "H\\d$"))) |> #, ": H=", rep(c(0.25, 0.24), n()/2))) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill=value), col=NA) +
  geom_sf(data=states, fill=NA, col="gray80", alpha=0.3) +
  facet_wrap(~label, nrow=2) +
  scale_fill_viridis_c(
    option = "inferno",
    limits = c(0.15, 0.71),       # tighten range
    oob = scales::squish          # values outside get squished to edges
  ) +
  labs(fill=NULL, title="Probability of severe illness (Lack imprinting protection, with age risks and pre-existing conditions)") +
  # labs(fill=NULL, title=paste0("probability unprotected with H=", p_get_anyway)) +
  theme_map()
p2

################################################################################
#### Georgia only p2

names(county_risk)
names(counties)
names(states)

p2GA <- county_risk |>
  pivot_longer(contains("severe")) |>
  mutate(
    fips  = stringr::str_pad(as.character(fips), 5, pad = "0"),
    label = stringr::str_extract(name, "H\\d$")
  ) |>
  left_join(counties, by = "fips") |>
  sf::st_as_sf() |>
  dplyr::filter(stringr::str_starts(fips, "13")) |>   # Georgia counties
  ggplot() +
  geom_sf(aes(fill = value), color = NA) +
  geom_sf(
    data = states |> dplyr::filter(fips == "13"),     # Georgia outline
    fill = NA, color = "gray80", linewidth = 0.4
  ) +
  facet_wrap(~label, nrow = 2) +
  scale_fill_viridis_c(
    option = "inferno",
    limits = c(0.15, 0.71),
    oob = scales::squish
  ) +
  labs(
    fill = NULL,
    title = "Proportion at risk of severe flu"
  ) +
  theme_map()

p2GA

ggsave(here("reconstructions-us-counties-main", "figs", "GA_only.png"), 
       plot = p2GA, width = 8, height = 6, dpi = 300)


#### get map numbers for GA only
ga_risk <- county_risk |>
  mutate(fips = stringr::str_pad(as.character(fips), 5, pad = "0")) |>
  dplyr::filter(stringr::str_starts(fips, "13")) |>
  dplyr::select(fips, risk_severe_H5, risk_severe_H7)

View(ga_risk)
names(ga_risk)

ga_risk <- ga_risk |>
  left_join(counties |> dplyr::select(fips, county), by = "fips")

View(ga_risk)

ga_risk_df <- ga_risk |>
  select(-geom)

View(ga_risk_df)
names(ga_risk_df)


# save data for use in python
write.csv(ga_risk_df, here("reconstructions-us-counties-main", "processed-data", 
                           "ga_data.csv"), row.names = FALSE)










summary(county_risk |> select(contains("severe")))

p1 <- county_risk |> 
  pivot_longer(contains("unprot")) |> 
  mutate(label = str_extract(name, "H\\d$"),
         label = recode(label, "H5" = "Group 1", "H7" = "Group 2")) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80") +
  facet_wrap(~label, nrow = 1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "Proportion\nvulnerable") +
  coord_sf(expand = FALSE) +                  # removes border padding
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),         # facets closer together
    plot.margin = margin(0, 0, 0, 0),         # more room for maps
    strip.background = element_blank(),       # less wasted space
    strip.text = element_text(
      margin = margin(b = 1, t = 1)
    ),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
# Show plot
p1
ggsave(here("reconstructions-us-counties-main", "figs", "first_plots_4.png"), plot = p1, width = 8, height = 6, dpi = 300)

county_risk2 <- county_risk %>%
  filter(fips != "12119")


p1 <- county_risk |> 
  pivot_longer(contains("unprot")) |> 
  mutate(label = str_extract(name, "H\\d$"),
         label = recode(label, "H5" = "Group 1", "H7" = "Group 2")) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80") +
  facet_wrap(~label, nrow = 1) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "Proportion\nvulnerable") +
  coord_sf(expand = FALSE) +                  # removes border padding
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),         # facets closer together
    plot.margin = margin(0, 0, 0, 0),         # more room for maps
    strip.background = element_blank(),       # less wasted space
    strip.text = element_text(
      margin = margin(b = 0.5, t = 1)
    ),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
# Show plot
p1
ggsave(here("reconstructions-us-counties-main", "figs", "first_plots_4.png"), plot = p1, width = 8, height = 6, dpi = 300)

county_risk2 <- county_risk %>%
  filter(fips != "12119")


p2 <- county_risk |> 
  pivot_longer(contains("severe")) |> 
  mutate(label = str_extract(name, "H\\d$"),
         label = recode(label, "H5" = "Group 1", "H7" = "Group 2")) |> 
  left_join(counties) |> 
  sf::st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = states, fill = NA, col = "gray80") +
  facet_wrap(~label, nrow = 1) +                                # match p1 layout
  scale_fill_viridis_c(option = "inferno") +
  labs(
    fill = "Proportion\nat risk"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),         # facets closer together
    plot.margin = margin(0, 0, 0, 0),         # reduce outside whitespace
    strip.background = element_blank(),
    strip.text = element_text(
      margin = margin(b = 1, t = 1)
    ),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

# Show
p2
ggsave(here("reconstructions-us-counties-main", "figs", "second_plots_5.png"), plot = p2, width = 8, height = 6, dpi = 600)


p1
p2

plot_grid(p1, p2, nrow=1)
ggsave(paste0("figs/HPAI-risk-map-", today(), ".pdf"), width=14.5, height=11)















p1 <- county_risk |> 
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

p2 <- county_risk |> 
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
p2_data <- county_risk |>
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
p2 <- p2_H5 + p2_H7
p2

################################################################################
# Plot FL and UT
# Filter states
fl_ut_states <- states |> 
  filter(abbr %in% c("FL", "UT"))

# Calculate bounding boxes for each state
fl_bbox <- sf::st_bbox(fl_ut_states |> filter(abbr == "FL"))
ut_bbox <- sf::st_bbox(fl_ut_states |> filter(abbr == "UT"))

# -------------------
# First plot: unprotected
# -------------------
fl_ut_unprot <- county_risk |> 
  pivot_longer(contains("unprot")) |> 
  mutate(label = str_c(str_extract(name, "H\\d$"), ": H=",
                       rep(c(0.25, 0.24), n()/2))) |> 
  left_join(counties, by = "fips") |> 
  filter(abbr %in% c("FL", "UT")) |> 
  sf::st_as_sf()

# Create separate plots for FL and UT
p1_fl <- fl_ut_unprot |> 
  filter(abbr == "FL") |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = fl_ut_states |> filter(abbr == "FL"), fill = NA, col = "gray80", linewidth = 0.4) +
  facet_wrap(~label, ncol = 2) +
  coord_sf(xlim = c(fl_bbox["xmin"], fl_bbox["xmax"]), 
           ylim = c(fl_bbox["ymin"], fl_bbox["ymax"]), 
           expand = FALSE) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = NULL, title = "FL - Probability Unprotected") +
  theme_map()

p1_ut <- fl_ut_unprot |> 
  filter(abbr == "UT") |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = fl_ut_states |> filter(abbr == "UT"), fill = NA, col = "gray80", linewidth = 0.4) +
  facet_wrap(~label, ncol = 2) +
  coord_sf(xlim = c(ut_bbox["xmin"], ut_bbox["xmax"]), 
           ylim = c(ut_bbox["ymin"], ut_bbox["ymax"]), 
           expand = FALSE) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = NULL, title = "UT - Probability Unprotected") +
  theme_map()

# -------------------
# Second plot: severe infection
# -------------------
fl_ut_severe <- county_risk |> 
  pivot_longer(contains("severe")) |> 
  mutate(label = str_c(str_extract(name, "H\\d$"), ": H=",
                       rep(c(0.25, 0.24), n()/2))) |> 
  left_join(counties, by = "fips") |> 
  filter(abbr %in% c("FL", "UT")) |> 
  sf::st_as_sf()

p2_fl <- fl_ut_severe |> 
  filter(abbr == "FL") |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = fl_ut_states |> filter(abbr == "FL"), fill = NA, col = "gray80", linewidth = 0.4) +
  facet_wrap(~label, ncol = 2) +
  coord_sf(xlim = c(fl_bbox["xmin"], fl_bbox["xmax"]), 
           ylim = c(fl_bbox["ymin"], fl_bbox["ymax"]), 
           expand = FALSE) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0.15, 0.71),
    oob = scales::squish
  ) +
  labs(fill = NULL, title = "FL - Probability of Severe Infection") +
  theme_map()

p2_ut <- fl_ut_severe |> 
  filter(abbr == "UT") |> 
  ggplot() +
  geom_sf(aes(fill = value), col = NA) +
  geom_sf(data = fl_ut_states |> filter(abbr == "UT"), fill = NA, col = "gray80", linewidth = 0.4) +
  facet_wrap(~label, ncol = 2) +
  coord_sf(xlim = c(ut_bbox["xmin"], ut_bbox["xmax"]), 
           ylim = c(ut_bbox["ymin"], ut_bbox["ymax"]), 
           expand = FALSE) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0.15, 0.71),
    oob = scales::squish
  ) +
  labs(fill = NULL, title = "UT - Probability of Severe Infection") +
  theme_map()

# -------------------
# Arrange in a grid: FL unprotected | UT unprotected
#                    FL severe      | UT severe
# -------------------
plot_grid(p1_fl, p1_ut, p2_fl, p2_ut, ncol = 2)


ggsave(
  filename = here("reconstructions-us-counties-main", "figs", "FL_UT_plots.PNG"), 
  width = 10, height = 8, dpi = 300
)







library(patchwork)

# Helper function: FL + UT plots without removing legend
make_fl_ut_plot <- function(data, outcome, H, fill_limits = NULL, fill_option = "inferno") {
  data_filtered <- data |> 
    filter(str_detect(label, H)) |> 
    filter(abbr %in% c("FL", "UT"))
  
  fl_plot <- data_filtered |> 
    filter(abbr == "FL") |> 
    ggplot() +
    geom_sf(aes(fill = value), col = NA) +
    geom_sf(data = fl_ut_states |> filter(abbr == "FL"), fill = NA, col = "gray80", linewidth = 0.4) +
    facet_wrap(~label, ncol = 1) +
    coord_sf(xlim = c(fl_bbox["xmin"], fl_bbox["xmax"]), 
             ylim = c(fl_bbox["ymin"], fl_bbox["ymax"]), expand = FALSE) +
    scale_fill_viridis_c(option = fill_option, limits = fill_limits, oob = scales::squish) +
    labs(fill = NULL) +
    theme_map() +
    theme(legend.position = "right")
  
  ut_plot <- data_filtered |> 
    filter(abbr == "UT") |> 
    ggplot() +
    geom_sf(aes(fill = value), col = NA) +
    geom_sf(data = fl_ut_states |> filter(abbr == "UT"), fill = NA, col = "gray80", linewidth = 0.4) +
    facet_wrap(~label, ncol = 1) +
    coord_sf(xlim = c(ut_bbox["xmin"], ut_bbox["xmax"]), 
             ylim = c(ut_bbox["ymin"], ut_bbox["ymax"]), expand = FALSE) +
    scale_fill_viridis_c(option = fill_option, limits = fill_limits, oob = scales::squish) +
    labs(fill = NULL) +
    theme_map() +
    theme(legend.position = "right")
  
  fl_plot + ut_plot
}

# --------------------------
# Fill limits per outcome
# --------------------------
fill_limits_unprot <- range(fl_ut_unprot$value, na.rm = TRUE)
fill_limits_severe <- range(fl_ut_severe$value, na.rm = TRUE)

library(patchwork)

# --------------------------
# Top row: unprotected (inferno)
# --------------------------
top_left  <- make_fl_ut_plot(fl_ut_unprot, " ", "H5", fill_limits_unprot, "inferno")
top_right <- make_fl_ut_plot(fl_ut_unprot, " ", "H7", fill_limits_unprot, "inferno")

top_row <- (top_left | top_right) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.direction = "vertical")

# --------------------------
# Bottom row: severe infection (viridis)
# --------------------------
bottom_left  <- make_fl_ut_plot(fl_ut_severe, " ", "H5", fill_limits_severe, "viridis")
bottom_right <- make_fl_ut_plot(fl_ut_severe, " ", "H7", fill_limits_severe, "viridis")

bottom_row <- (bottom_left | bottom_right) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.direction = "vertical")

library(patchwork)

library(patchwork)
library(ggplot2)
library(grid)

# Assume top_row and bottom_row are your 2×1 or 1×2 patchwork objects

# Create text grobs for row labels
top_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Probability Unprotected", size = 6, fontface = "bold") +
  theme_void()

bottom_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Probability of Severe Infection", size = 6, fontface = "bold") +
  theme_void()

# Stack them with your plot rows
final_plot <- top_label / top_row / bottom_label / bottom_row + 
  plot_layout(heights = c(0.1, 1, 0.1, 1), guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

final_plot


# --------------------------
# Save
# --------------------------
ggsave(
  filename = here("reconstructions-us-counties-main", "figs", "FL_UT_plots2.png"),
  plot = final_plot,
  width = 12,
  height = 12,
  dpi = 300
)
















# Create the full dataset once
p2_data <- county_risk |> 
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
