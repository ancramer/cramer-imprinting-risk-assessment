## Generate imprinting probabilities
rm(list = ls())
getwd()
source("code/analysis/calculation_funs.R")
source("code/analysis/data_import_funs.R")

library(tidyverse)

#################
# Imprinting probs, US only, through 2024
get_imprinting_probabilities(observation_years = 2024)  

