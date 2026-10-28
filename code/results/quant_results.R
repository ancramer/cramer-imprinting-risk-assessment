#### load in packages ##########################################################
library(tidyverse)
library(here)

##### load in data #############################################################
county_risk <- read.csv(here("data", "processed-data", "map_data_numbers.csv"))

# load county-level vulnerability data
map_data <- read.csv(here("data", 
                          "processed-data", 
                          "map_data_numbers.csv"))

# fix fips to 5-character strings
map_data <- map_data %>%
  mutate(fips = str_pad(as.character(fips), width = 5, pad = "0"))


# load county population totals
county_pop_totals <- read.csv(here("data", "processed-data", "county_pop_totals.csv")) %>%
  mutate(fips = sprintf("%05d", as.integer(county_fips.x))) %>%
  select(fips, county_total)


# join population to vulnerability data
map_data_pop <- map_data %>%
  left_join(county_pop_totals, by = "fips") %>%
  rename(county_pop = county_total)

# calculate weighted average percent vulnerable, group by state
# create state fips column from first 2 digists of county fips
map_data_pop <- map_data_pop %>%
  mutate(state_fips = substr(fips, 1, 2))

#### find risk for a specific county ###########################################
county_risk |> 
  filter(fips == "12119") |> 
  glimpse()

# find number of counties that are at a risk of 0.30 or greater
county_risk |> 
  filter(risk_severe_H5 > 0.35) |> 
  nrow()/3222

county_risk |> 
  filter(risk_severe_H7 > 0.35) |> 
  nrow()/3222

# highest-risk counties (Group 1)
top_group1_counties <- map_data_pop %>%
  select(fips, county_pop, risk_unprot_H5) %>%
  arrange(desc(risk_unprot_H5)) %>%
  slice_head(n = 10)

top_group1_counties

# lowest-risk counties (Group 1)
bottom_group1_counties <- map_data_pop %>%
  select(fips, county_pop, risk_unprot_H5) %>%
  arrange(risk_unprot_H5) %>%
  slice_head(n = 10)

bottom_group1_counties

# highest-risk counties (Group 2)
top_group2_counties <- map_data_pop %>%
  select(fips, county_pop, risk_unprot_H7) %>%
  arrange(desc(risk_unprot_H7)) %>%
  slice_head(n = 10)

top_group2_counties

# lowest-risk counties (Group 2)
bottom_group2_counties <- map_data_pop %>%
  select(fips, county_pop, risk_unprot_H7) %>%
  arrange(risk_unprot_H7) %>%
  slice_head(n = 20)

bottom_group2_counties


#### Vulnerability maps results ################################################

#### calculate weighted average percent vulnerable, group by state
# create state fips column from first 2 digists of county fips
map_data_pop <- map_data_pop %>%
  mutate(state_fips = substr(fips, 1, 2))

state_weighted_vul <- map_data_pop %>%
  group_by(state_fips) %>%
  summarise(
    group1_weighted = sum(risk_unprot_H5 * county_pop, na.rm = TRUE) / sum(county_pop, na.rm = TRUE),
    group2_weighted = sum(risk_unprot_H7 * county_pop, na.rm = TRUE) / sum(county_pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group1_weighted = group1_weighted,
    group2_weighted = group2_weighted
  )

state_lookup <- data.frame(
  state_fips = c("01","02","04","05","06","08","09","10","11","12","13","15","16","17","18","19",
                 "20","21","22","23","24","25","26","27","28","29","30","31","32","33","34",
                 "35","36","37","38","39","40","41","42","44","45","46","47","48","49","50",
                 "51","53","54","55","56"),
  state_abbr = c("AL","AK","AZ","AR","CA","CO","CT","DE","DC","FL","GA","HI","ID","IL","IN","IA",
                 "KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM",
                 "NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA",
                 "WV","WI","WY")
)

state_weighted_vul <- state_weighted_vul %>%
  left_join(state_lookup, by = "state_fips") %>%
  select(state_abbr, everything())

# sort states by Group 1 weighted vulnerability (high to low)
state_sorted_group1 <- state_weighted_vul %>%
  arrange(desc(group1_weighted))

state_sorted_group1

# sort states by Group 2 weighted vulnerability (high to low)
state_sorted_group2 <- state_weighted_vul %>%
  arrange(desc(group2_weighted))

state_sorted_group2



#### Median age results ########################################################

# read in data (this is the same as in figs_for_supp.R)
here()
age2 <- read.csv(here("data", "processed-data", "ACSST5Y2023.S0101-Data.csv"))

# take off first row
age2 <- age2[-1, ]

# clean and prepare median age
age2_median <- age2 %>%
  select(GEO_ID, NAME, median_age = S0101_C01_032E) %>%
  mutate(
    median_age = as.numeric(median_age),
    fips = stringr::str_pad(substr(GEO_ID, 10, 14), width = 5, pad = "0")
  )

# get properly formatted county data
us_counties <- us_map(regions="counties") # county shapefiles for plotting
us_counties

# join median age to counties
plot_data <- left_join(us_counties, age2_median, by = "fips")

# sort data to find youngest and oldest median counties
# youngest counties
plot_data %>%
  filter(!is.na(median_age)) %>%
  arrange(median_age) %>%
  slice(1:10) %>%
  select(NAME, abbr, median_age)

# oldest counties
plot_data %>%
  filter(!is.na(median_age)) %>%
  arrange(desc(median_age)) %>%
  slice(1:10) %>%
  select(NAME, abbr, median_age)

# calcualted oldest and youngest states
state_median_age <- plot_data %>%
  filter(!is.na(median_age)) %>%
  group_by(abbr) %>%
  summarise(
    state_median_age = median(median_age),
    .groups = "drop"
  )

# pull youngest state
youngest_state <- state_median_age %>%
  slice_min(state_median_age, n = 1)

youngest_state

# pull oldest state
oldest_state <- state_median_age %>%
  slice_max(state_median_age, n = 1)

oldest_state


