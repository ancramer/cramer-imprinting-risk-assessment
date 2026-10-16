
############
#changed to call the new intensity master file for added years
############

get_p_infection_year = function(birth_year,
                                observation_year, ## Year of data collection, which matters if observation_year is shortly after birth_year
                                baseline_annual_p_infection = 0.28, ## Baseline annual probability of infection
                                max_year){ ## max calendar year for which to output estimates
  ## Function to calculate probs of first exposure in year x, given birth in year y
  ## INPUTS
  ##    - year in which an individual was born (birth.year)
  ##    - year in which the individual became infected with bird flu (infection year)
  ## OUTPUTS
  ##    - vector of 13 probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 13th year of life (age 12)
  stopifnot(observation_year <= max_year)
  cat('Reading annual intensities from ../processed_data/intensity_master_new.csv. See Gostic et al. 2016 for details.- File has been updated')
  intensity_df = read_csv('data/processed-data/intensity_master_new.csv', show_col_types = F)  ##changed to take the new intensitymaster file
  # Weighted attack rate = annual prob infection weighted by circulation intensity
  weighted.attack.rate = baseline_annual_p_infection*(intensity_df$intensity)
  names(weighted.attack.rate) = intensity_df$year
  ################# Calculations ---------------
  possible_imprinting_years = birth_year:min(birth_year+12, observation_year) #Calendar years of first infection (ages 0-12)
  nn = length(possible_imprinting_years) # How many possible years of first infection? (should be 13)
  valid_attack_rates = weighted.attack.rate[as.character(possible_imprinting_years)] #Get weighted attack rates corresponding to possible years of first infection
  attack_rate_complements = matrix(rep(1-valid_attack_rates, nn), nn, nn, byrow = T)
  ## Create matrices of 0s and 1s, which will be used below to vectorize and speed calculations
  infection_year = not_infection_year = matrix(0, nn, nn)
  diag(infection_year) = 1   #Fill in diagonal of one with 1s for years of first infection
  not_infection_year[lower.tri(not_infection_year)] = 1  #Fill in sub-diagonal for all the years since birth in which the individual escaped infection. 
  # Exact probability of escaping infection in the previous (x-1) years, and becoming infected in year x
  prod.mat = (valid_attack_rates*infection_year)+(attack_rate_complements*not_infection_year)
  #Fill in upper triangle with 1s to make multiplication possible
  prod.mat[upper.tri(prod.mat)] = 1
  #Take product across rows
  p_ij = apply(prod.mat, 1, prod)
  p_ij # Output probability of first infection in year i given birth year
}


#########
#This function calculates imprinting probabilities for the US only 
###########

get_imprinting_probabilities <- function(observation_years) {
  ## This function calculates imprinting probabilities for the US only
  
  # Get max observation year
  max_year <- max(observation_years)
  stopifnot(max_year <= as.numeric(format(Sys.Date(), '%Y')))
  
  # Define birth years and infection years
  birth_years <- 1918:max_year
  nn_birth_years <- length(birth_years)
  
  # Initialize matrices for probabilities
  H1N1_probs <- matrix(NA, nrow = length(observation_years), ncol = length(birth_years),
                       dimnames = list(observation_years, birth_years))
  H2N2_probs <- naive_probs <- H3N2_probs <- H1N1_probs
  
  ## Get US-specific epidemiological data
  us_epi_data <- get_US_data(max_year)
  
  ## Extract the fractions of circulating influenza subtypes by birth year
  H1.frac <- as.numeric(us_epi_data["A/H1N1", as.character(birth_years)])
  H2.frac <- as.numeric(us_epi_data["A/H2N2", as.character(birth_years)])
  H3.frac <- as.numeric(us_epi_data["A/H3N2", as.character(birth_years)])
  names(H1.frac) <- names(H2.frac) <- names(H3.frac) <- as.character(birth_years)
  
  ## Loop across observation years and birth years
  for (jj in seq_along(observation_years)) {
    observation_year <- observation_years[jj]
    
    for (ii in seq_along(birth_years)) {
      birth_year <- birth_years[ii]
      
      if (birth_year > observation_year) next  # Skip future birth years
      
      n_infection_years <- min(12, observation_year - birth_year)  # Up to age 12
      
      # Get probability of first infection in each year
      inf.probs <- get_p_infection_year(birth_year, observation_year,
                                        baseline_annual_p_infection = 0.28, max_year)
      
      # Normalize probabilities if all 13 years have passed
      if (length(inf.probs) == 13) inf.probs <- inf.probs / sum(inf.probs)
      
      # Fill in imprinting probabilities
      valid_infection_years <- as.character(seq(birth_year, birth_year + n_infection_years))
      
      if (all(valid_infection_years %in% names(H1.frac))) {
        H1N1_probs[jj, ii] <- sum(inf.probs * H1.frac[valid_infection_years])
        H2N2_probs[jj, ii] <- sum(inf.probs * H2.frac[valid_infection_years])
        H3N2_probs[jj, ii] <- sum(inf.probs * H3.frac[valid_infection_years])
      } else {
        H1N1_probs[jj, ii] <- NA
        H2N2_probs[jj, ii] <- NA
        H3N2_probs[jj, ii] <- NA
      }
      
      naive_probs[jj, ii] <- round(1 - sum(inf.probs, na.rm = TRUE), digits = 8)
    }
  }
  
  ## Ensure the birth years are displayed in descending order
  descending_order <- rev(birth_years)
  H1N1_probs <- H1N1_probs[, as.character(descending_order)]
  H2N2_probs <- H2N2_probs[, as.character(descending_order)]
  H3N2_probs <- H3N2_probs[, as.character(descending_order)]
  naive_probs <- naive_probs[, as.character(descending_order)]
  
  ## Validate that all probabilities sum to 1
  total <- H1N1_probs + H2N2_probs + H3N2_probs + naive_probs
  if (any(!(round(total, 4) %in% c(0, 1)))) {
    warning("Weights do not sum to 1")
  }
  
  ## Normalize so that sum of probabilities is exactly 1
  total[total == 0] <- 1  # Avoid division by zero
  H1N1_probs <- H1N1_probs / total
  H2N2_probs <- H2N2_probs / total
  H3N2_probs <- H3N2_probs / total
  naive_probs <- naive_probs / total
  
  ## Return results
  return(list(
    H1N1_probs = H1N1_probs,
    H2N2_probs = H2N2_probs,
    H3N2_probs = H3N2_probs,
    naive_probs = naive_probs
  ))
}


