#Load packages
library(tidyverse)
library(zippeR)
library(zctaCrosswalk)
library(readxl)
library(here)

#read in risk data from Remi, which is comorbidities by zip code
risk <- read_csv(here("data", "raw-data", "county_risk_data", 'COVID High risk population per age group per zip code.csv'))

################################################################################
#### Put counties to zip code risk data, so it's now risk tied to counties

#add column zip of zipcodes, take out old ZCTA5 (same data with bad name)
risk2 <- risk %>%
  mutate(zip = zi_repair(ZCTA5)) %>% 
  select(-ZCTA5)

#isolate just connecticut, deal with this later
risk_ct <- risk2 %>%
  filter(str_starts(zip, "06"))

#save everything non-CT
risk_other <- risk2 %>%
  filter(!str_starts(zip, "06")) #deal with this later


#links zipcode and country into a table, acts as a reference
crosswalk <- zctaCrosswalk::get_zcta_crosswalk() 

#fix name so it will merge
crosswalk <- crosswalk %>%
  rename(zip = zcta) #changed from zip_code to zip so it would merged below

#keep all columns from risk2 (non CT), add info from crosswalk that matches by zip
merged <- risk_other %>%
  left_join(crosswalk, by = "zip")

################################################################################
#### Fixing NAs identified above

#read in file I made that contains information for NAs
fixes <- read_csv(here("data", "raw-data", "county_risk_data", "fixes_for_missing_data.csv"))

#remove unecessary info
fixes <- fixes %>%
  select(-state) %>%
  select(-county_fips_short) %>%
  select(-"...7") %>%
  select(-"...8") %>%
  select(-"...9") %>%
  select(-"...10")

#merge the fixes with original 'merged'
merged2 <- merged %>%
  left_join(fixes, by = "zip") %>%
  mutate(
    county_name = ifelse(is.na(county_name.x), county_name.y, county_name.x),
    county_fips = ifelse(is.na(county_fips.x), county_fips.y, county_fips.x),
    state_fips = ifelse(is.na(state_fips.x), state_fips.y, state_fips.x)
  ) %>%
  select(-county_name.x, -county_fips.x, -state_fips.x,
         -county_name.y, -county_fips.y, -state_fips.y)

########################################################################################
#### Fix Connecticut counties, because they have been remade

#Old connecticut zips are present, but maps use new CT zips
#county_avg2 %>%
#  filter(county_fips=="09009")

#read in new CT county fips with zips
CT <- read_excel(here("data", "raw-data", "county_risk_data", "CT-ZIP-to-PlanningRegion.xlsx"))

#rename ZCTA_5 to zip in CT data
CT <- CT %>%
  mutate(zip = zi_repair(ZCTA_5)) %>% 
  select(-ZCTA_5)

#isolated just connecticut earlier

#merge CT risk by zip to new county reference list
risk_ct_merged <- merge(risk_ct, CT, by = "zip", all.x = TRUE)

#rename columns to match the risk2 columns
risk_ct_merged <- risk_ct_merged %>%
  rename(county_name = PlanningRegion) %>%
  rename(county_fips = PlanningRegion_GeoID) %>%
  mutate(county_fips = paste0("0", county_fips)) %>% #put a 0 in front of all county fips because they got messed up by excel formatting
  mutate(state_fips = "09") #add in state fips

#recombine CT with all other states
merged2 <- bind_rows(merged2, risk_ct_merged)

####################################################################################
#### Weighted by zip code n score
# Accounts for zip codes that occur in multiple counties, assigned a score "n" for how many times they appear

#load population count data
pop <- read.csv(here("data", "raw-data", "county_risk_data", "ZCTA Population per age group.csv"))

# Generated below, note this takes like 48 hours+ to run, so it's been saved into this dup_results.csv
dup_results <- read.csv(here("data", "raw-data", "county_risk_data", "dup_results.csv"))
## Connecticut already contains new county names, might need to be edited a little but they're in good shape

####################################################################################
#### Zip code duplicates: Weighing population by their proportion areas
#rename zip codes in dup_results
dup_results <- dup_results %>%
  rename(zip = ZIP_CODE) %>%
  mutate(
    zip = str_pad(
      as.character(zip),
      width = 5,
      side = "left",
      pad = "0"
    )
  )

#get zip codes from GeographicAreaName to make it easier bc it's not working
pop$zip <- sub("ZCTA5 ", "", pop$GeographicAreaName)

# Add in state fips
dup_results <- dup_results %>%
  mutate(state_fips = sprintf("%02d", county_fips %/% 1000))

#merge population data with dup proportion data
merged_pop_dups <- left_join(pop, dup_results, by = "zip")

#remove NA rows, which are the zips that are normal
merged_pop_dups <- merged_pop_dups %>%
  filter(!is.na(proportion_of_zip))

#multiply the age groups by their proportions
merged_pop_dups <- merged_pop_dups %>%
  mutate(across(c(X0_4, X5_9, X10_14, X15_19, X20_24, X25_29, X30_34, 
                  X35_39, X40_44, X45_49, X50_54, X55_59, X60_64, X65_69, 
                  X70_74, X75_79, X80_84, X85.), 
                ~ .x * proportion_of_zip))

# THIS SHOULD BE READY TO COMBINE WITH REGULAR COUNTIES BELOW.

####################################################################################
#### Weighted by population per zip code
#### This makes zip codes that are entirely in one county be weighted towards the county average 
#### by their respective population within the zip, and dups are weighted above.

#merge by zipcode
#left join by zip code to add in county/fips info fo all zip that are in one county only
#remove age risk columns bc we're trying to get population. We will deal with risk later
merged_pop <- left_join(pop, merged2 %>% select(zip, state_fips, county_fips, county_name), by = "zip")

#Fix county names for dups since we sdid singles about, change name of column and add word "County"
merged_pop_dups <- merged_pop_dups %>%
  rename(county_name = COUNTY) %>%
  mutate(county_name = paste0(county_name, " County"))

# Now we will combine the dups with pop scaled by their % area to the set of single zip codes.

# Convert types & add missing columns to make things combine well
merged_pop_dups <- merged_pop_dups %>%
  mutate(
    zip = as.character(zip),
    county_fips = as.character(county_fips),
    state_fips = as.character(state_fips)
  )

merged_pop <- merged_pop %>%
  mutate(
    zip = as.character(zip),
    county_fips = as.character(county_fips),
    state_fips = as.character(state_fips),
    proportion_of_zip = 1  # fill in 1 bc these zips are 100% in one county
  )

# Combine datasets, prioritize merged_pop_dups
combined_pop <- bind_rows(merged_pop_dups, merged_pop) %>%
  arrange(zip, county_fips, desc(proportion_of_zip)) %>%
  distinct(zip, county_fips, .keep_all = TRUE)


#NAs are puerto rico (which I thought I removed earlier but must have gotten back in somehow)
combined_pop <- combined_pop %>%
  filter(!grepl("^00", as.character(zip)))


# Now we will add in the risks
combined_pop_risks <- left_join(combined_pop, merged2, by = "zip")

# add 0's to the start of county_fips that dropped them again
combined_pop_risks <- combined_pop_risks %>%
  mutate(county_fips.x = sprintf("%05s", as.character(county_fips.x)))
unique(combined_pop_risks$proportion_of_zip)


# Step 1: scale Total by proportion_of_zip for each ZIP
step1 <- combined_pop_risks %>%
  #filter(!is.na(county_fips.x) & county_fips.x != "000NA") %>%  # remove unmatched ZIPs
  mutate(Total = Total * proportion_of_zip)  # create new column with scaled population

# Step 2: sum scaled Total per county
county_totals <- step1 %>%
  filter(!is.na(county_fips.x) & county_fips.x != "000NA") %>%  # remove unmatched ZIPs
  # for ZIPs in multiple counties, scale Total by proportion_of_zip
  group_by(zip, county_fips.x) %>%
  summarise(scaled_total = mean(Total) * mean(proportion_of_zip), .groups = "drop") %>%
  # sum over counties
  group_by(county_fips.x) %>%
  summarise(county_total = sum(scaled_total), .groups = "drop")
# save this to use elsewhere: Total population per county
write.csv(county_totals, here("data", "processed-data", "county_pop_totals.csv"), row.names = FALSE)

# Step 2: Add proportion of total county population per ZIP
combined_weighted <- combined_pop_risks %>%
  left_join(county_totals, by = "county_fips.x") %>%
  mutate(proportion_pop = Total / county_total)

# Step 3: Define column groups
pop_age_cols <- c("X0_4", "X5_9", "X10_14", "X15_19", "X20_24", "X25_29", "X30_34",
                  "X35_39", "X40_44", "X45_49", "X50_54", "X55_59", "X60_64", "X65_69",
                  "X70_74", "X75_79", "X80_84", "X85.")
risk_age_cols <- c("0_0.5", "0.5_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
                   "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75+")

# Step 4: Collapse to 1 row per county
county_avg_both <- combined_weighted %>%
  group_by(county_fips.x) %>%
  summarise(across(all_of(c(pop_age_cols, risk_age_cols)),
                   ~sum(.x * proportion_pop, na.rm = TRUE) / sum(proportion_pop, na.rm = TRUE),
                   .names = "avg_{.col}"),
            county_name = first(county_name.x),
            state_fips = first(state_fips.x),
            .groups = "drop")

# avg_X is weighted average number of people in the age group in the county
# avg_ is the weighted risk for the age group for the county

#remove 1st row which was an average of weird NA counties that have barely any people living in them and no risk data
county_avg_both2 <- county_avg_both[-1, ]

########################################################################################
#### Add in CDC high risk age groups
# 2 and under and 65 and older go to risk=1

county_avg_both_with_ages <- county_avg_both2 %>%
  select(-"avg_X0_4", -"avg_X5_9", -"avg_X10_14", -"avg_X15_19", -"avg_X20_24", -"avg_X25_29", -"avg_X30_34",
         -"avg_X35_39", -"avg_X40_44", -"avg_X45_49", -"avg_X50_54", -"avg_X55_59", -"avg_X60_64", -"avg_X65_69",
         -"avg_X70_74", -"avg_X75_79", -"avg_X80_84", -"avg_X85.") %>%
  mutate(
    avg_0_2 = 1,
    avg_2_4 = avg_0.5_4,
    avg_65_69 = 1,
    avg_70_74 = 1,
    `avg_75+` = 1
  ) %>%
  select(-avg_0_0.5, -avg_0.5_4) %>% #drop old columns
  relocate(avg_0_2, .after = 1) %>%  #move avg_0_2 to be second
  relocate(avg_2_4, .after = avg_0_2)  #move avg_2_4 to be third


county_avg_both_with_only_comorb <- county_avg_both2 %>%
  select(-"avg_X0_4", -"avg_X5_9", -"avg_X10_14", -"avg_X15_19", -"avg_X20_24", -"avg_X25_29", -"avg_X30_34",
         -"avg_X35_39", -"avg_X40_44", -"avg_X45_49", -"avg_X50_54", -"avg_X55_59", -"avg_X60_64", -"avg_X65_69",
         -"avg_X70_74", -"avg_X75_79", -"avg_X80_84", -"avg_X85.") %>%
  mutate(
    avg_0_2 = avg_0.5_4,
    avg_2_4 = avg_0.5_4,
    # avg_65_69 = 1,
    #avg_70_74 = 1,
    #`avg_75+` = 1
  ) %>%
  select(-avg_0_0.5, -avg_0.5_4) %>% #drop old columns
  relocate(avg_0_2, .after = 1) %>%  #move avg_0_2 to be second
  relocate(avg_2_4, .after = avg_0_2)  #move avg_2_4 to be third

#### Plot for comorbs only ####
#plot process starts with making it long
county_avg_both_with_only_comorb_long <- county_avg_both_with_only_comorb %>%
  pivot_longer(
    cols = starts_with("avg_"),
    names_to = "age_group",
    values_to = "risk"
  )

# Correct the order of the age groups bc right now it's messed up
age_levels <- c(
  "avg_0_2", "avg_2_4", "avg_5_9", "avg_10_14", "avg_15_19",
  "avg_20_24", "avg_25_29", "avg_30_34", "avg_35_39", "avg_40_44",
  "avg_45_49", "avg_50_54", "avg_55_59", "avg_60_64",
  "avg_65_69", "avg_70_74", "avg_75+"
)

# Apply that order to the long dataset
county_avg_both_with_only_comorb_long <- county_avg_both_with_only_comorb_long %>%
  mutate(age_group = factor(age_group, levels = age_levels))

# remove avg label
county_avg_both_with_only_comorb_long <- county_avg_both_with_only_comorb_long %>%
  mutate(
    # set factor with explicit order
    age_group = factor(age_group, levels = age_levels, ordered = TRUE),
    # remove "avg_" for display but preserve ordering
    age_group = fct_relabel(age_group, ~ str_remove(.x, "^avg_"))
  )

c <- ggplot(county_avg_both_with_only_comorb_long, 
            aes(x = age_group, y = risk)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.3, color = "#238443", size = 1.5) + 
  geom_violin(fill = "#f7fcb9", alpha = 1, color = "black") +
  labs(
    title = "Risk among counties by age group (comorbidities)",
    x = "Age Group",
    y = "Proportion of population at high risk (with at least 1 comorbidity)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = function(x) gsub("avg_", "", gsub("_", "-", x)))

c
# Save dataset to use elsewhere (in scatterplots)
write.csv(
  county_avg_both_with_only_comorb_long,
  here("data", "processed-data", "county_avg_both_with_only_comorb_long.csv"),
  row.names = FALSE
)


#Save plot
#ggsave(here("risk_by_county", "figures", "comorb_dist3.png"), plot = c)

#### Plot of age effects only ####

county_avg_with_only_ages <- county_avg_both_with_ages %>%
  mutate(
    avg_2_4 = 0,
    avg_5_9 = 0,
    avg_10_14 = 0,
    avg_15_19 = 0,
    avg_20_24 = 0,
    avg_25_29 = 0,
    avg_30_34 = 0,
    avg_35_39 = 0,
    avg_40_44 = 0,
    avg_45_49 = 0,
    avg_50_54 = 0,
    avg_55_59 = 0,
    avg_60_64 = 0)

#making it long
county_avg_with_only_ages <- county_avg_with_only_ages %>%
  pivot_longer(
    cols = starts_with("avg_"),
    names_to = "age_group",
    values_to = "risk"
  )

# remove avg label
county_avg_with_only_ages <- county_avg_with_only_ages %>%
  mutate(
    # set factor with explicit order
    age_group = factor(age_group, levels = age_levels, ordered = TRUE),
    # remove "avg_" for display but preserve ordering
    age_group = fct_relabel(age_group, ~ str_remove(.x, "^avg_"))
  )

#plot
a <- ggplot(county_avg_with_only_ages, aes(x = age_group, y = risk)) +
  geom_violin(fill = "firebrick", alpha = 0.6, color = "black") +
  labs(
    title = "Risk by age group (age)",
    x = "Age Group",
    y = "Proportion at high risk (age based risk)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1))
a

#Save plot
#ggsave(here("risk_by_county", "figures", "age_risk_dist2.png"), plot = a)

#### Plot for both ####
# Pivot to long format
county_avg_both_with_ages <- county_avg_both_with_ages %>%
  pivot_longer(
    cols = starts_with("avg_"),
    names_to = "age_group",
    values_to = "risk"
  )

age_levels <- c(
  "avg_0_2", "avg_2_4", "avg_5_9", "avg_10_14", "avg_15_19",
  "avg_20_24", "avg_25_29", "avg_30_34", "avg_35_39", "avg_40_44",
  "avg_45_49", "avg_50_54", "avg_55_59", "avg_60_64",
  "avg_65_69", "avg_70_74", "avg_75+"
)

# Apply that order to the long dataset
county_avg_both_with_ages <- county_avg_both_with_ages %>%
  mutate(age_group = factor(age_group, levels = age_levels))


q <- ggplot(county_avg_both_with_ages, aes(x = age_group, y = risk)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.3, color = "#a6317d", size = 1.5) +
  geom_violin(fill = "#fcfdbf", alpha = 1, color = "black") +
  labs(
    # title = "Risk among counties by age group (Comorbs + age effects)",
    x = "Age group",
    y = "Proportion of population \nwith a risk factor"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = function(x) gsub("_", "-", gsub("avg_", "", x))) +
  theme_minimal(base_size = 14) +  # base size affects all text elements
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
q

#Save plot
ggsave(here("code", "figs", "comorb_and_age_risk_by_age_group.png"), plot = q)
#save to send to Bren or to put into maps. Final version, contains risk weighted by age group, and
# Population per county is weighted by how much geographical space each zip takes 
# up within county borders.
write.csv(county_avg_both_with_ages, here("data", "processed-data", "county_risk_weighted_with_dups.csv"), row.names = FALSE)
# avg_X is weighted average number of people in the age group in the county
# avg_ is the weighted risk for the age group for the county
# made somewhere around line 378

