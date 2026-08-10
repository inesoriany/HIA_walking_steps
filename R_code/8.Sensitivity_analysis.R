#################################################
########     SENSITIVITY ANALYSIS         #######
#################################################

# Modified parameters - files outputted :
  # Walking speed 5.3 km/h
  # Former DRF


################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  purrr,        # Loop
  stringr,      # Text extraction
  srvyr,        # Survey
  survey,
  ggplot2,      # Data visualization
  progress      # Progression bar
)


################################################################################################################################
################################################################################################################################
#                                                       ALTERNATIVE DRF                                                        #
################################################################################################################################
################################################################################################################################

# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Disability weights distribution table
dw_distrib_table <- import(here("data_clean", "Diseases", "dw_distrib_table.xlsx"))


# Import functions
source(here("R_code", "Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("R_code", "Parameters.R"))

# Diseases considered
alt_dis_vec <- c("mort", "bc", "cc", "cvd", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly")


################################################################################################################################
#                                                    3. DATA PREPARATION                                                       #
################################################################################################################################
