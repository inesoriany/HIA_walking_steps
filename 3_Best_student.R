#################################################
#############     GOOD STUDENTS     #############
#################################################



################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  srvyr,        # Survey
  tidyr,        # Table - Data organization, extraction
  tidyverse,    # Data manipulation and visualization
  ggplot2       # Plotting
)


################################################################################################################################
#                                                    2. IMPORT DATA                                                            #
################################################################################################################################
# Walkers dataset
emp_walkers <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Walking trips dataset
emp_walk_trip <- import(here("data_clean", "EMP_dis_walk_trips.xlsx"))


# RR by step, simulated dose-response relationships
rr_central_table <- import(here("data_clean", "Diseases", "DRF", "rr_central_interpolated.rds"))


# Disability weights
dw_table <- import(here("data", "dw_table.xlsx"))



# Import functions
source(here("0_Functions.R"))


################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "bc", "cvd", "cancer", "diab2", "dem", "dep")

# Bound
bound_vec <- c("mid", "low", "up")


################################################################################################################################
#                                                   4. AGE DISTRIBUTION                                                        #
################################################################################################################################
# Survey design ponderated by day
jour <- emp_walkers %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)


# Walking pyramid : Age distribution of walking volume for France
distrib_walk_EMP2019 <- jour %>% 
  group_by(age_grp10) %>% 
  summarise(mean_ind = survey_mean(step_commute, na.rm = TRUE))
