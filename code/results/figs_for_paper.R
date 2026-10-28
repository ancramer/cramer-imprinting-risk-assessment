#### load packages #############################################################
library(tidyverse)
library(here)

#### Figure 1 ##################################################################
# further edits made by Daniel Humphrey, graphic designer at the Odum School of Ecology
# at University of Georgia, Athens GA

# set working directory and load functions
here()
source(here("code", "processing", "calculation_funs.R"))
source(here("code", "processing", "data_import_funs.R"))

# get imprinting probabilities for US using Maya's updated code
US <- get_imprinting_probabilities(observation_years = 2010:2024)

# extract strain-specific data
US_H1N1 <- US$H1N1_probs
US_H2N2 <- US$H2N2_probs
US_H3N2 <- US$H3N2_probs
US_naive <- US$naive_probs

# add strain labels, change format from matrix to df, add year label
# H1N1
US_H1N1 <- as.data.frame(US$H1N1_probs)
US_H1N1$obs_year <- rownames(US$H1N1_probs)
US_H1N1$strain <- "H1N1"

# H2N2
US_H2N2 <- as.data.frame(US$H2N2_probs)
US_H2N2$obs_year <- rownames(US$H2N2_probs)
US_H2N2$strain <- "H2N2"

# H3N2
US_H3N2 <- as.data.frame(US$H3N2_probs)
US_H3N2$obs_year <- rownames(US$H3N2_probs)
US_H3N2$strain <- "H3N2"

# Naive
US_naive <- as.data.frame(US$naive_probs)
US_naive$obs_year <- rownames(US$naive_probs)
US_naive$strain <- "naive"

# combine all into one data frame
unitedstates <- bind_rows(US_H1N1, US_H2N2, US_H3N2, US_naive)

# long format, use only columns with 4-digit years bc there was issues
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

# Reorder strains so naive is on top
unitedstates_long$strain <- factor(unitedstates_long$strain, levels = c("naive", "H1N1", "H2N2", "H3N2"))

# Plot
unitedstates_plot <- unitedstates_long |>
  filter(obs_year == 2024) |>
  ggplot(aes(x = birth_year, y = probability, fill = strain)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Birth Year",
    y = "Proportion of Population",
    title = "Imprinting Probabilities in the US",
    fill = "Strain"
  ) +
  scale_fill_manual(values = c(
    "H1N1" = "#e08214",
    "H2N2" = "#fdb863",
    "H3N2" = "#8073ac",
    "naive" = "grey80"
  )) +
  scale_x_continuous(
    breaks = c(1918, 1957, 1968, 1977, 2009, 2024),
    labels = c("1918", "1957", "1968", "1977", "2009", "2024"),
    minor_breaks = NULL
  ) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),  # bottom labels, angled
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

# Show plot
unitedstates_plot

# Save plot
ggsave(
  filename = here::here("code", "figs", "US_birth_cohorts.png"),
  plot = unitedstates_plot,
  width = 8,
  height = 4,
  dpi = 300)

################################################################################
#### Panel Figure 2: A is median age map, B is vulnerability age line plot, C is vulnerability map set
################################################################################

#### A: making median age map ##################################################
here()
library(usmap)
library(ggplot2)
library(dplyr)
library(here)
library(viridis)
library(sf)

# Load data
age2 <- read.csv(here("data", "processed-data", "ACSST5Y2023.S0101-Data.csv"))

age2 <- age2[-1, ]

# Clean and prepare median age
age2_median <- age2 %>%
  select(GEO_ID, NAME, median_age = S0101_C01_032E) %>%
  mutate(
    median_age = as.numeric(median_age),
    fips = stringr::str_pad(substr(GEO_ID, 10, 14), width = 5, pad = "0")
  )

# Get properly formatted county data
us_counties <- us_map(regions="counties") # county shapefiles for plotting

# Join median age to counties
plot_data <- left_join(us_counties, age2_median, by = "fips")

# Plot
# If using sf package
# Create state boundaries by dissolving counties
us_states <- us_counties %>% 
  group_by(abbr) %>%  # using the state abbreviation column
  summarise(geometry = st_union(geom), .groups = 'drop')  # note: your geometry column is named 'geom'

# Then use in your plot
p_age <- ggplot(plot_data) +
  geom_sf(aes(fill = median_age), col = NA) +
  geom_sf(data = us_counties, fill = NA, col = "NA", alpha = 0.6) +
  geom_sf(data = us_states, fill = NA, col = "grey80", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Median Age", na.value = "gray90") +
  labs(title = "") +
  theme_void()
p_age

p_age_formatted <- plot_data |>
  mutate(dummy = "Map A") |>  # Add dummy variable for faceting
  ggplot() +
  geom_sf(aes(fill = median_age), col = NA) +                          # fill counties
  geom_sf(data = us_counties, fill = NA, col = NA, alpha = 0.6) +     # county boundaries
  geom_sf(data = us_states, fill = NA, col = "gray80", linewidth = 0.25) + # thinner state borders
  facet_wrap(~dummy, nrow = 1) +                                      # faceting like p2
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "gray90",
    name = "Median \nage"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    strip.background = element_blank(),
    strip.text = element_blank(),   # hide facet label
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

# Show the plot
p_age_formatted

#### B: making vulnerability age line plot #####################################

#### C: making vulnerability maps ##############################################


#### making panel ##############################################################
# Apply existing themes
a <- a + map_theme
b <- b + line_plot_theme
c <- c + map_theme_c  # Use theme with larger left margin

# Add titles in the panel using ggdraw
a_with_title <- ggdraw() +
  draw_label("Age distribution by county across the United States", x = 0.5, y = 1, 
             hjust = 0.5, vjust = 1.5, size = 18) +
  draw_plot(a, y = 0, height = 0.95)

b_with_title <- ggdraw() +
  draw_label(
    "Proportion unprotected by imprint by birth year", 
    x = 0.5, y = 1, 
    hjust = 0.5, vjust = 1.2,  # move slightly above top
    size = 18
  ) +
  draw_plot(b, y = 0, height = 0.9)  # reduce plot height to give more room for title

c_with_title <- ggdraw() +
  draw_label("Proportion of population unprotected by imprint in the United States", 
             x = 0.5, y = 1, hjust = 0.5, vjust = 1.5, size = 18) +
  draw_plot(c, y = 0, height = 0.95)

# Top row: A + B with labels
top_row <- plot_grid(
  a_with_title, b_with_title,
  ncol = 2,
  rel_widths = c(1.5, 1.2),  # Give A more space, B slightly wider
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  align = "v"
)

# Bottom row: C with label
bottom_row <- plot_grid(
  c_with_title,
  ncol = 1,
  labels = "C",
  label_size = 16,
  label_fontface = "bold"
)

# Combine rows
final_panel <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)  # Equal heights
)

# Display and save
final_panel

ggsave(
  filename = here("reconstructions-us-counties-main", "figs", "fig_2_2.png"),
  plot = final_panel,
  width = 13,  # Increased width
  height = 8,
  dpi = 300
)

################################################################################
#### Panel Figure 3: A is comorb map, B is risk factor violin, C is high risk sev flu map set
################################################################################

#### A: making comorb map ######################################################

library(usmap)
here()
age_comorb_risksev <- read.csv(here("data", "processed-data", "age_comorb_risksev.csv"))
# make fips 5-digit strings
age_comorb_risksev <- age_comorb_risksev %>%
  mutate(fips = sprintf("%05s", fips))

# get county geometries (sf)
county_map <- usmap::us_map(regions = "counties")

# merge comorbidity data
county_map <- county_map %>%
  left_join(age_comorb_risksev %>% select(fips, risk_comorb), by = "fips")

# chsnge geometry for plotting
county_sf <- sf::st_as_sf(county_map, sf_column_name = "geom")

# in the original p_comorbs creation add
p_comorbs <- county_sf |>
  mutate(dummy = "Map A") |>  # add a dummy variable
  ggplot() +
  geom_sf(aes(fill = risk_comorb), col = NA) +
  geom_sf(
    data = usmap::us_map(regions = "states") |> sf::st_as_sf(),
    fill = NA,
    col = "gray80",
    linewidth = 0.3
  ) +
  facet_wrap(~dummy, nrow = 1) +  # Add facet like p2
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "gray90"
  ) +
  labs(
    fill = "Comorbidity\nrate"
  ) +
  coord_sf(expand = FALSE) +
  theme_map() +
  theme(
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    strip.background = element_blank(),
    strip.text = element_blank(),  # Hide the facet label
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
p_comorbs

#### B: making risk factor violin ##############################################
#### C: making high risk for severe flu maps ###################################

#### making these into a panel #################################################
a <- p_comorbs # in Comorb_explorations.R
b <- q         # in county_age_dists.R in the risk_by_county folder
c <- p2        # in county-riskmap.R

# Apply existing themes
a <- a + map_theme
b <- b + violin_axis_theme
c <- c + map_theme_c  # Use theme with larger left margin

# Add titles in the panel using ggdraw
a_with_title <- ggdraw() +
  draw_label(
    "Comorbidity distribution by county across the United States",
    x = 0.5, y = 1, 
    hjust = 0.5, vjust = 1.5,
    size = 18
  ) +
  draw_plot(a, y = 0, height = 0.95)

b_with_title <- ggdraw() +
  draw_label(
    "Prevalence of demographic vulnerability factors\n(age and comorbidities)", 
    x = 0.5, y = 1, 
    hjust = 0.5, vjust = 1.2,  # slightly above the top of the plot
    size = 18,
    lineheight = 1.1
  ) +
  draw_plot(b, y = 0, height = 0.85)  # shrink plot slightly to give more room for 2-line title

c_with_title <- ggdraw() +
  draw_label(
    "Proportion of population at high risk of severe influenza in the United States",
    x = 0.5, y = 1, 
    hjust = 0.5, vjust = 1.5,
    size = 18
  ) +
  draw_plot(c, y = 0, height = 0.95)

# Top row: A + B with labels
top_row <- plot_grid(
  a_with_title, b_with_title,
  ncol = 2,
  rel_widths = c(1.5, 1.2),  # Give A more space, B slightly wider
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  align = "v"
)

# Bottom row: C with label
bottom_row <- plot_grid(
  c_with_title,
  ncol = 1,
  labels = "C",
  label_size = 16,
  label_fontface = "bold"
)

# Combine rows
final_panel <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)  # Equal heights
)

# Display and save
final_panel

ggsave(
  filename = here("reconstructions-us-counties-main", "figs", "fig_3_5.png"),
  plot = final_panel,
  width = 13,  # Increased width
  height = 8,
  dpi = 300
)

################################################################################
#### Panel Figure 4: A & B & C are scatterplots, D & E are state/census region results
################################################################################

#### A, B, C: making scatterplots ##############################################
library(ggplot2)
library(ggnewscale)
age_comorb_risksev <- read.csv(here("data", "processed-data", "age_comorb_risksev.csv"))

# Ensure median age is numeric
age_comorb_risksev$median_age <-
  as.numeric(age_comorb_risksev$median_age)

# Shared colors
group_colors <- c(
  "Group 1" = "#e08214",
  "Group 2" = "#8073ac"
)

line_colors <- c(
  "Group 1" = "#ad6511",
  "Group 2" = "#574d78"
)

# Function for R-squared labels
r2_label <- function(model, digits = 3) {
  paste0(
    "R² = ",
    formatC(
      summary(model)$r.squared,
      digits = digits,
      format = "f"
    )
  )
}

########## A ###########
# Fit models in the same direction as the plotted axes
lm_age_H5 <- lm(
  risk_severe_H5 ~ median_age,
  data = age_comorb_risksev
)

lm_age_H7 <- lm(
  risk_severe_H7 ~ median_age,
  data = age_comorb_risksev
)

plot_age2 <- ggplot(age_comorb_risksev) +
  geom_point(
    aes(
      x = median_age,
      y = risk_severe_H7,
      color = "Group 2"
    ),
    alpha = 0.15,
    size = 2
  ) +
  geom_point(
    aes(
      x = median_age,
      y = risk_severe_H5,
      color = "Group 1"
    ),
    alpha = 0.15,
    size = 2
  ) +
  geom_abline(
    slope = coef(lm_age_H5)[["median_age"]],
    intercept = coef(lm_age_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]],
    linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_age_H7)[["median_age"]],
    intercept = coef(lm_age_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]],
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 50,
    y = 0.60,
    label = r2_label(lm_age_H5),
    color = line_colors[["Group 1"]],
    hjust = 1,
    size = 7
  ) +
  annotate(
    "text",
    x = 50,
    y = 0.70,
    label = r2_label(lm_age_H7),
    color = line_colors[["Group 2"]],
    hjust = 1,
    size = 7
  ) +
  scale_color_manual(
    values = group_colors,
    name = "Group"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 4)
    )
  ) +
  labs(
    x = "Median age",
    y = "Proportion at high risk\nfor severe influenza"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    legend.position = "bottom"
  )

a <- plot_age2
a

########## B ##########
lm_risk_H5 <- lm(
  risk_severe_H5 ~ risk_comorb,
  data = age_comorb_risksev
)

lm_risk_H7 <- lm(
  risk_severe_H7 ~ risk_comorb,
  data = age_comorb_risksev
)

plot_risk2 <- ggplot(age_comorb_risksev) +
  geom_point(
    aes(
      x = risk_comorb,
      y = risk_severe_H7,
      color = "Group 2"
    ),
    alpha = 0.15,
    size = 2
  ) +
  geom_point(
    aes(
      x = risk_comorb,
      y = risk_severe_H5,
      color = "Group 1"
    ),
    alpha = 0.15,
    size = 2
  ) +
  geom_abline(
    slope = coef(lm_risk_H5)[["risk_comorb"]],
    intercept = coef(lm_risk_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]],
    linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_risk_H7)[["risk_comorb"]],
    intercept = coef(lm_risk_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]],
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 0.485,
    y = 0.60,
    label = r2_label(lm_risk_H5),
    color = line_colors[["Group 1"]],
    hjust = 1,
    size = 7
  ) +
  annotate(
    "text",
    x = 0.50,
    y = 0.70,
    label = r2_label(lm_risk_H7),
    color = line_colors[["Group 2"]],
    hjust = 1,
    size = 7
  ) +
  scale_color_manual(
    values = group_colors,
    name = "Group"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 4)
    )
  ) +
  labs(
    x = "Proportion with comorbidity",
    y = "Proportion at high risk\nfor severe influenza"
  ) +
  coord_cartesian(xlim = c(0.2, 0.7)) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = grid::unit(1.2, "cm"),
    legend.position = "bottom"
  )

b <- plot_risk2
b
########## C ##########
lm_unprotected_H5 <- lm(
  risk_severe_H5 ~ risk_unprot_H5,
  data = age_comorb_risksev
)

lm_unprotected_H7 <- lm(
  risk_severe_H7 ~ risk_unprot_H7,
  data = age_comorb_risksev
)

# Data-based label positions
x_H5 <- quantile(
  age_comorb_risksev$risk_unprot_H5,
  probs = 0.01,
  na.rm = TRUE
)

x_H7 <- quantile(
  age_comorb_risksev$risk_unprot_H7,
  probs = 0.01,
  na.rm = TRUE
)

y_H5 <- quantile(
  age_comorb_risksev$risk_severe_H5,
  probs = 0.97,
  na.rm = TRUE
)

y_H7 <- quantile(
  age_comorb_risksev$risk_severe_H7,
  probs = 0.97,
  na.rm = TRUE
)

plot_unimprinted <- ggplot(age_comorb_risksev) +
  # Group 2/H7 layer
  geom_point(
    aes(
      x = risk_unprot_H7,
      y = risk_severe_H7,
      color = median_age
    ),
    alpha = 0.6,
    size = 2
  ) +
  scale_color_gradient(
    low = "#d8daeb",
    high = "#2d004b",
    name = "Median age"
  ) +
  ggnewscale::new_scale_color() +
  # Group 1/H5 layer
  geom_point(
    aes(
      x = risk_unprot_H5,
      y = risk_severe_H5,
      color = median_age
    ),
    alpha = 0.6,
    size = 2
  ) +
  scale_color_gradient(
    low = "#fee0b6",
    high = "#7f3b08",
    name = "Median age"
  ) +
  geom_abline(
    slope = coef(lm_unprotected_H7)[["risk_unprot_H7"]],
    intercept = coef(lm_unprotected_H7)[["(Intercept)"]],
    color = line_colors[["Group 2"]],
    linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(lm_unprotected_H5)[["risk_unprot_H5"]],
    intercept = coef(lm_unprotected_H5)[["(Intercept)"]],
    color = line_colors[["Group 1"]],
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = x_H5,
    y = y_H5,
    label = r2_label(lm_unprotected_H5),
    color = line_colors[["Group 1"]],
    hjust = 0,
    size = 7
  ) +
  annotate(
    "text",
    x = x_H7,
    y = y_H7,
    label = r2_label(lm_unprotected_H7),
    color = line_colors[["Group 2"]],
    hjust = 0,
    size = 7
  ) +
  labs(
    x = "Proportion unimprinted/unprotected",
    y = "Proportion at high risk\nfor severe influenza"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

c <- plot_unimprinted
c


#### D, E: making state and census region results ##############################
########## D ###########
# Read in data from maps
here()
map_numbers <- read.csv(here("data", "processed-data", "map_data_numbers.csv"))

# Make sure FIPS is character and padded to 5 digits
map_numbers <- map_numbers %>%
  mutate(fips = str_pad(as.character(fips), width = 5, pad = "0"))

# Get state FIPS
map_numbers <- map_numbers %>%
  mutate(state_fips = substr(fips, 1, 2))

length(unique(map_numbers$state_fips))

# State fips to state names
state_lookup <- c(
  "01" = "Alabama", "02" = "Alaska", "04" = "Arizona", "05" = "Arkansas",
  "06" = "California", "08" = "Colorado", "09" = "Connecticut",
  "10" = "Delaware", "12" = "Florida", "13" = "Georgia", "15" = "Hawaii",
  "16" = "Idaho", "17" = "Illinois", "18" = "Indiana", "19" = "Iowa",
  "20" = "Kansas", "21" = "Kentucky", "22" = "Louisiana", "23" = "Maine",
  "24" = "Maryland", "25" = "Massachusetts", "26" = "Michigan",
  "27" = "Minnesota", "28" = "Mississippi", "29" = "Missouri",
  "30" = "Montana", "31" = "Nebraska", "32" = "Nevada",
  "33" = "New Hampshire", "34" = "New Jersey", "35" = "New Mexico",
  "36" = "New York", "37" = "North Carolina", "38" = "North Dakota",
  "39" = "Ohio", "40" = "Oklahoma", "41" = "Oregon",
  "42" = "Pennsylvania", "44" = "Rhode Island",
  "45" = "South Carolina", "46" = "South Dakota", "47" = "Tennessee",
  "48" = "Texas", "49" = "Utah", "50" = "Vermont",
  "51" = "Virginia", "53" = "Washington", "54" = "West Virginia",
  "55" = "Wisconsin", "56" = "Wyoming"
)

state_lookup_df <- tibble(
  state_fips = names(state_lookup),
  state = unname(state_lookup)
)

# Join state names
map_numbers <- map_numbers %>%
  left_join(state_lookup_df, by = "state_fips")

# Remove non-states sorry Puerto Rico
map_numbers <- map_numbers %>%
  filter(!state_fips %in% c("72", "11"))

# State abbreviations
state_abbrev_lookup <- tibble(
  state = state.name,
  state_abb = state.abb
)
names(map_numbers)
map_numbers <- map_numbers %>%
  left_join(state_abbrev_lookup, by = c("state" = "state")) %>%
  mutate(state = state_abb) %>%
  select(-state_abb)

# Assign census regions

map_numbers <- map_numbers %>%
  mutate(
    census_region = case_when(
      state %in% c("ME","NH","VT","MA","RI","CT") ~ "New England",
      state %in% c("NY","NJ","PA") ~ "Middle Atlantic",
      state %in% c("DE","MD","VA","WV","NC","SC","GA","FL") ~ "South Atlantic",
      state %in% c("KY","TN","MS","AL") ~ "East South Central",
      state %in% c("AR","LA","OK","TX") ~ "West South Central",
      state %in% c("OH","IN","IL","MI","WI") ~ "East North Central",
      state %in% c("MN","IA","MO","ND","SD","NE","KS") ~ "West North Central",
      state %in% c("MT","ID","WY","CO","NM","AZ","UT","NV") ~ "Mountain",
      state %in% c("WA","OR","CA","AK","HI") ~ "Pacific",
      TRUE ~ NA_character_
    ),
    census_region_main = case_when(
      census_region %in% c("New England","Middle Atlantic") ~ "Northeast",
      census_region %in% c("South Atlantic","East South Central","West South Central") ~ "South",
      census_region %in% c("East North Central","West North Central") ~ "Midwest",
      census_region %in% c("Mountain","Pacific") ~ "West",
      TRUE ~ NA_character_
    )
  )

# merge population data in

county_pop_totals <- read.csv(here("data", "processed-data", "county_pop_totals.csv"))

county_pop_totals <- county_pop_totals %>%
  mutate(county_fips.x = sprintf("%05d", as.integer(county_fips.x)))

data_for_regions_states <- map_numbers %>%
  left_join(county_pop_totals, by = c("fips" = "county_fips.x")) %>%
  rename(county_pop = county_total)

# State level weighted risk

state_risks <- data_for_regions_states %>%
  group_by(state, state_fips, census_region_main) %>%
  summarise(
    across(starts_with("risk"),
           ~ weighted.mean(.x, county_pop, na.rm = TRUE)),
    total_pop = sum(county_pop, na.rm = TRUE),
    .groups = "drop"
  )

# Plot for group 1

state_risksH5 <- state_risks %>%
  mutate(state = fct_reorder(state, risk_severe_H5))

data_for_regions_statesH5 <- data_for_regions_states %>%
  mutate(state = factor(state, levels = levels(state_risksH5$state)))

region_colors <- c(
  "Northeast" = "#99DDFFFF",
  "West"      = "#DCB0F2FF",
  "South"     = "#F89C74FF",
  "Midwest"   = "#C9DB74FF"
)

G1_results <- ggplot() +
  geom_violin(
    data = data_for_regions_statesH5,
    aes(y = state, x = risk_severe_H5, fill = census_region_main),
    color = NA,
    alpha = 0.9,
    scale = "width"
  ) +
  geom_jitter(
    data = data_for_regions_states,
    aes(y = state, x = risk_severe_H5),
    color = "grey40",
    alpha = 0.3,
    height = 0.25,
    size = 1
  ) +
  geom_point(
    data = state_risks,
    aes(y = state, x = risk_severe_H5),
    color = "grey30",
    alpha = 0.8,
    size = 2
  ) +
  scale_fill_manual(values = region_colors) +
  labs(
    y = "State",
    x = "Proportion at high risk \nfor severe group 1 disease",
    title = "County-level Risk by State for Group 1",
    fill = "U.S. Census Region"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 7),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.position = "bottom"
  )
G1_results

########## E ###########

# Reorder states by weighted average (H7)
state_risksH7 <- state_risks %>%
  mutate(state = fct_reorder(state, risk_severe_H7))

# Apply same ordering to counties
data_for_regions_statesH7 <- data_for_regions_states %>%
  mutate(state = factor(state, levels = levels(state_risksH7$state)))

# --- Region colors (same as Group 1) ---
region_colors <- c(
  "Northeast" = "#99DDFFFF",
  "West"      = "#DCB0F2FF",
  "South"     = "#F89C74FF",
  "Midwest"   = "#C9DB74FF"
)

# --- Plot ---
G2_results <- ggplot() +
  # Violin plots colored by region
  geom_violin(
    data = data_for_regions_statesH7,
    aes(y = state, x = risk_severe_H7, fill = census_region_main),
    color = NA,
    alpha = 0.9,
    scale = "width"
  ) +
  
  # County-level dots
  geom_jitter(
    data = data_for_regions_states,
    aes(y = state, x = risk_severe_H7),
    color = "grey40",
    alpha = 0.3,
    height = 0.25,
    size = 1
  ) +
  
  # State-level weighted averages
  geom_point(
    data = state_risks,
    aes(y = state, x = risk_severe_H7),
    color = "grey30",
    alpha = 0.8,
    size = 2
  ) +
  
  scale_fill_manual(values = region_colors) +
  
  labs(
    y = "State",
    x = "Proportion at high risk \nfor severe group 2 disease",
    title = "County-level Risk by State for Group 2",
    fill = "U.S. Census Region"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 7),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.position = "bottom"
  )
G2_results

# Assign plots
a <- plot_age2
a

b <- plot_risk2
b

c <- plot_unimprinted
c

d <- G1_results
d
e <- G2_results
e

# Base theme

base_theme <- theme(
  plot.title = element_blank(),
  
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text = element_text(size = 12),
  
  legend.position = "bottom",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.6, "cm"),
  
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0)
)

a <- a + base_theme
b <- b + base_theme
c <- c + base_theme
d <- d + base_theme
e <- e + base_theme

# Remove y-axis labels from top row (shared label used instead)

a <- a + labs(y = NULL)
b <- b + labs(y = NULL)
c <- c + labs(y = NULL)

# Extract legends (BEFORE removing them)

legend_groups <- get_legend(a)  # orange/purple group legend
legend_age    <- get_legend(c)  # median age legend
legend_bottom <- get_legend(d)  # state/census legend

# Remove legends from plots

a <- a + theme(legend.position = "none")
b <- b + theme(legend.position = "none")
c <- c + theme(legend.position = "none")

d <- d + theme(legend.position = "none")
e <- e + theme(legend.position = "none")

# Shared y-axis label

shared_y_axis <- ggdraw() +
  draw_label(
    "Proportion at high risk\nfor severe influenza",
    angle = 90,
    size = 14
  )

panel_label_size <- 14

# TOP PANELS (A, B, C)

top_panels <- plot_grid(
  a, b, c,
  ncol = 3,
  
  labels = c("A", "B", "C"),
  label_size = panel_label_size,
  label_fontface = "bold",
  
  # IMPORTANT FIX: keep labels visible (no clipping)
  label_x = -0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  
  align = "hv",
  axis = "tblr"
)

top_main <- plot_grid(
  shared_y_axis,
  top_panels,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

# Legend row under A/B/C (group + age split)

top_legend_row <- plot_grid(
  
  legend_groups,
  
  NULL,
  
  # wrap legend_age in a padded grob
  
  plot_grid(
    
    NULL,
    
    legend_age,
    
    ncol = 2,
    
    rel_widths = c(0.15, 1)
    
  ),
  
  ncol = 3,
  
  rel_widths = c(1, 0.2, 1)
  
)

top_row <- plot_grid(
  top_main,
  top_legend_row,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

# BOTTOM PANELS (D, E)

bottom_main <- plot_grid(
  d, e,
  ncol = 2,
  
  labels = c("D", "E"),
  label_size = panel_label_size,
  label_fontface = "bold",
  
  # IMPORTANT FIX: keep labels visible
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

# FINAL FIGURE

final_panel <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  
  # keeps D/E tall (50-state axis preserved)
  rel_heights = c(0.78, 1.9)
)

# Display

final_panel

# Save

ggsave(
  filename = here("reconstructions-us-counties-main", "figs", "fig_4_11.png"),
  plot = final_panel, width = 12, height = 10, dpi = 300, limitsize = FALSE)
