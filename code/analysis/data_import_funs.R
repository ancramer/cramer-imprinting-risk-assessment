#########
# Removed functions not used from Gostic et al. repo, to make this US only

#show_available_countries
#show_available_regions
#parse_country_names
#parse_region_names
#get_WHO_region
##########


test_rowsums_subtype <- function(H1N1, H2N2, H3N2){
 stopifnot(all(H1N1+H3N2+H2N2 == 1 | is.na(H1N1+H3N2+H2N2)))
}

test_rowsums_group <- function(group1, group2){
  stopifnot(all(group1+group2 == 1 | is.na(group1+group2)))
}

check_years <- function(years, max_year){
  stopifnot(1997:max_year %in% years)
}


###############
# This stays the same: since it always was US specific 
##############
get_template_data <- function(){
  ## From 1918-1976, only 1 flu A subtype circulated
  ## From 1977-1996, country-specific data is not available
  ## For this time period, output a template of country and region-independent data
  ## Columns:
  ##  * year: 1918:1996
  ##  * `A/H1N1`: fraction of flu A positive specimens of subytpe H1N1
  ##  * `A/H2N2`: fraction of flu A positive specimens of subytpe H2N2
  ##  * `A/H3N2`: fraction of flu A positive specimens of subytpe H3N2
  ##  * `B`: placeholder -- all NA
  ##  * group1: fraction of flu A positive specimens of group1 (H1N1+H2N2)
  ##  * group2: fraction of flu A positive specimens of group2 (H3N2)
  ## Load epidemiologica data from 1977-1996
  ## From Thompson et al. JAMA, 2003 (JAMA. 2003;289(2):179-186. doi:10.1001/jama.289.2.179)
  Thompson_data = read_csv('data/processed-data/Thompson_data.csv', show_col_types = FALSE) %>%
    mutate(`A/H1N1` = n_H1N1/n_A,
           `A/H2N2` = 0,
           `A/H3N2` = n_H3N2/n_A,
            B = NA) %>%
    select(year, `A/H1N1`, `A/H2N2`, `A/H3N2`, B) %>%
    mutate(data_from = 'Thompsonetal_JAMA_2003')
  ## Set up template and fill in years of single-subtype circulation
  ## Then return
  template = bind_rows(
    ## H1N1 era
    tibble(year = 1918:1956, 
           `A/H1N1` = 1,
           `A/H2N2` = 0,
           `A/H3N2` = 0,
           B = NA,
           data_from = 'Historical_assumption'),
    ## H2N2 era
    tibble(year = 1957:1967, 
           `A/H1N1` = 0,
           `A/H2N2` = 1,
           `A/H3N2` = 0,
           B = NA,
           data_from = 'Historical_assumption'),
    ## H3N2 era
    tibble(year = 1968:1976, 
           `A/H1N1` = 0,
           `A/H2N2` = 0,
           `A/H3N2` = 1,
           B = NA,
           data_from = 'Historical_assumption'),
    ## cocirculation 1977-1996
    Thompson_data
    ) %>%
    mutate(group1 = `A/H1N1` + `A/H2N2`,
           group2 = `A/H3N2`,
           A = `A/H1N1`+`A/H2N2`+`A/H3N2`) %>%
    select(year, starts_with('A'), starts_with('B'), starts_with('group'), data_from)
  ## Test
  test_rowsums_group(template$group1, template$group2)
  test_rowsums_subtype(template$`A/H1N1`, template$`A/H2N2`, template$`A/H3N2`)
  return(template)
}

############
# This replaces get_country_inputs_1997_to_present function 
#############
get_US_inputs_1997_to_present <- function(max_year) { 
  # Define file path
  file_path <- 'data/raw-data/americas/Americas_1997_2024.csv'
  
  # Load data
  us_data <- read_csv(file_path, show_col_types = FALSE) %>%
    group_by(ISO_YEAR) %>%  # Aggregate by country and year
    summarise(
      n_H1N1 = sum(AH1, na.rm = TRUE) + sum(AH1N12009, na.rm = TRUE),
      n_H3N2 = sum(AH3, na.rm = TRUE),
      n_A = n_H1N1 + n_H3N2,
      n_BYam = sum(BYAM, na.rm = TRUE),
      n_BVic = sum(BVIC, na.rm = TRUE),
      n_B = n_BYam + n_BVic
    ) %>%
    ungroup()
  
  # Print a message for verification
  cat("Data successfully loaded for the United States.\n")
  
  # Ensure years are valid
  check_years(us_data$ISO_YEAR, max_year)
  
  return(us_data)
}


###############
# This replaces get_country_data function 
###############
get_US_data <- function(max_year, min_samples_per_year = 30) {
  ## Load the template for historical data (1918-1996)
  template <- get_template_data()
  
  ## Get US-specific data from 1997-present
  us_data <- get_US_inputs_1997_to_present(max_year) %>%
    #filter(n_A >= min_samples_per_year) %>%
    mutate(data_from = 'US-specific data')
  
  ## Calculate the proportions of each subtype
  formatted_data <- us_data %>%
    mutate(`A/H1N1` = n_H1N1 / n_A,
           `A/H2N2` = 0,
           `A/H3N2` = n_H3N2 / n_A,
           B = NA) %>%
    mutate(group1 = `A/H1N1` + `A/H2N2`,
           group2 = `A/H3N2`,
           A = `A/H1N1` + `A/H2N2` + `A/H3N2`) %>%
    rename(year = ISO_YEAR) %>%
    select(year, starts_with('A'), starts_with('B'), starts_with('group'), data_from)
  
  ## Combine with the historical template data
  full_outputs <- bind_rows(template, formatted_data)
  
  ## Run validation checks
  check_years(years = full_outputs$year, max_year = max_year)
  test_rowsums_group(full_outputs$group1, group2 = full_outputs$group2)
  test_rowsums_subtype(full_outputs$`A/H1N1`, full_outputs$`A/H2N2`, full_outputs$`A/H3N2`)
  
  ## Convert to matrix with column names as years
  output_matrix <- full_outputs %>%
    as.matrix() %>%
    t()
  
  colnames(output_matrix) <- output_matrix['year', ]
  
  return(output_matrix)
}



