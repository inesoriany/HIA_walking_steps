###############################################
############ CREATING DATASETS ################
###############################################

# This code puts together datasets at the individual and trip level, combining information on:
# - walking, driving exposure
# - disease incidence



# Files outputted :
  # EMP_walkers.xlsx
  # EMP_drivers.xlsx


################################################################################################################################
################################################################################################################################
#                                                      WALKERS DATABASE                                                        #
################################################################################################################################
################################################################################################################################


################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################

pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data manipulation
  forcats,      # Factor conversion
  tidyr,        # Table - Data organization, extraction
  tidyverse,    # Data management, ggplot included
  epikit,       # Age categories creation
  janitor,      # De-duplication
  survey        # Survey management
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# EMP 2019 : walking and driving exposure
emp_walk_ind <- import(here("data", "emp_dataset_walk_individual.xlsx")) 
emp_walk_trip <- import(here("data", "emp_dataset_walk_trip.xlsx")) 
emp_car_trip <- import(here("data", "emp_dataset_car_trip.xlsx")) 

# Diseases incidence data
diseases <- import(here("data", "GBD_diseases.xlsx"), sheet = "Central")

# INSEE mortality data
insee <- import(here("data", "INSEE_2019.RDS"))



# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                     3. SET PARAMETERS                                                        #
################################################################################################################################
# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")
morbi_vec = c("cvd", "bc", "cancer", "diab2", "dem", "dep")



################################################################################################################################
#                                                4. DISEASES INCIDENCE DATASET                                                 #
################################################################################################################################
diseases <- diseases %>% 
  rename(age_grp = age_grp.x)


# 10 year age categories
diseases_10 <- diseases %>%
  mutate(age_min = as.numeric(sub("-.*", "", age_grp)),
         age_grp10 = paste0(floor(age_min / 10) * 10, "-", floor(age_min / 10) * 10 + 9)) %>%
  group_by(age_grp10, sex) %>%
  summarise(pop_age_grp10 = sum(pop_age_grp, na.rm = TRUE),
            across(all_of(paste0(morbi_vec, "_incidence")), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop") 



################################################################################################################################
#                                           5. CREATION OF SUBSET OF EMP DATASET                                               #
################################################################################################################################
# --------------------------------------
# INDIVIDUAL
# --------------------------------------

# Total walking distance
walkers <- emp_walk_ind %>% 
  mutate(nbkm_intermodal_walk = intermodal_walk_time * walk_speed / 60,
         nbkm_tot_walking = nbkm_main_walk + nbkm_intermodal_walk,
         nbkm_tot_walking_jour = nbkm_tot_walking * pond_jour / (pond_indc * 7))

# Create walkers dataset combing diseases incidence and walking exposure for each individual
walkers <- walk_dataset(walkers, diseases, diseases_10, insee, morbi_vec, 
                           walk_dist_var = "nbkm_tot_walking", 
                           walk_dist_jour_var = "nbkm_tot_walking_jour", 
                           step_length = step_length, 
                           walk_speed = walk_speed)



# --------------------------------------
# TRIPS
# --------------------------------------

# Total walking distance per trip
walking_trip <- emp_walk_trip  %>% 
  mutate(nbkm_intermodal_walk = intermodal_walk_time * walk_speed / 60,
         nbkm_tot_walking = nbkm_main_walk + nbkm_intermodal_walk,
         nbkm_tot_walking_jour = nbkm_tot_walking * pond_jour / (pond_indc * 7))

# Create walkers dataset combing diseases incidence and walking exposure for each individual
walking_trip <- walk_dataset(walking_trip, diseases, diseases_10, insee, morbi_vec, 
                           walk_dist_var = "nbkm_tot_walking", 
                           walk_dist_jour_var = "nbkm_tot_walking_jour", 
                           step_length = step_length, 
                           walk_speed = walk_speed)

# Area type depend on density
walking_trip <- walking_trip %>%
  mutate(
    area_type = case_when(
      densitecom_ori == 1    ~ "urban",
      densitecom_ori == 2    ~ "periurban",
      TRUE                   ~ "rural"),
    area_type = factor(area_type, levels = c("rural", "periurban", "urban"))
  )



################################################################################################################################
#                                                   5. DISEASE WALKERS EMP SUBSET                                                      #
################################################################################################################################

# Individual
walkers_long <- walkers %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )

# Trips
walking_trip_long <- walking_trip %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )


################################################################################################################################
#                                                       6. EXPORT DATA                                                         #
################################################################################################################################

# Tables of walkers
  export(walkers, here("data_clean", "EMP_walkers.xlsx"))
  export(walkers_long, here("data_clean", "EMP_dis_walkers.xlsx"))

# Tables of walking trips
  export(walking_trip, here("data_clean", "EMP_walking_trips.xlsx"))
  export(walking_trip_long, here("data_clean", "EMP_dis_walking_trips.xlsx"))




################################################################################################################################
################################################################################################################################
#                                                      DRIVERS DATABASE                                                        #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                  4. CREATION OF SUBSET OF EMP SUBSET WITH ONLY VARIABLES NEEDED                              #
################################################################################################################################

# Selecting only variables of interests for drivers / removing non-relevant variables 
emp_drivers <- emp_subset %>% 
  select(-nbkm_walking,
         -day_time)


# Associate drive speed
emp_drivers <- emp_drivers %>%
  mutate(drive_speed = case_when(
    area_type == c("urban", "periurban") ~ urban_car_speed,
    tuu2017_res == 8     ~ paris_car_speed,
    TRUE                 ~ rural_car_speed
    )
  )


# Day time spent walking if those car distances were walked (min)
emp_drivers <- emp_drivers %>% 
  mutate(day_time_shift = mdisttot_fin1*60 / walk_speed)

# Day time spent walking if those car distances were walked (min)
emp_drivers <- emp_drivers %>% 
  mutate(week_time_shift = 7*mdisttot_fin1*60 / walk_speed)


# Daily steps if those car distances were walked
emp_drivers <- emp_drivers %>% 
  mutate(day_step_commute_shift = mdisttot_fin1/step_length)


################################################################################################################################
#                                                   5. DISEASE EMP SUBSET                                                      #
################################################################################################################################

emp_drivers_long <- emp_drivers %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )



################################################################################################################################
#                                                    6. EXPORT EMP SUBSET                                                      #
################################################################################################################################

# Tables of drivers
  export(emp_drivers, here("data_clean", "EMP_drivers.xlsx"))
  export(emp_drivers_long, here("data_clean", "EMP_dis_drivers.xlsx"))

