###############################################
############ CREATING DATASETS ################
###############################################

# This code puts together datasets at the individual and trip level, combining information on:
# - walking, driving exposure
# - disease incidence


# Files outputted :
  # EMP_walkers.xlsx : dataset at the individual level with information on walking exposure and diseases incidence for each individual
  # EMP_dis_walkers.xlsx : Long format with one row per disease
  # EMP_walking_trips.xlsx : dataset at the trip level with information on walking exposure and diseases incidence for each trip
  # EMP_dis_walking_trips.xlsx : Long format with one row per disease
  # EMP_car_trips.xlsx : dataset at the trip level with information on driving exposure and diseases incidence for each trip
  # EMP_dis_car_trips.xlsx : Long format with one row per disease




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


# Disease incidence
 # Central values
  incid_table_mid <- import(here("data", "GBD_diseases.xlsx"), sheet = "Central")
  
  # IC95 lower
  incid_table_low <- import(here("data", "GBD_diseases.xlsx"), sheet = "Lower")
  
  # IC95 upper
  incid_table_up <- import(here("data", "GBD_diseases.xlsx"), sheet = "Upper") 



# INSEE mortality data
insee <- import(here("data", "INSEE_2019.RDS"))


# RR by step, simulated dose-response relationships
rr_central_table <- import(here("data_clean", "Diseases", "DRF", "rr_central_interpolated.rds"))
rr_distrib_table <- import(here("data_clean", "Diseases", "DRF", "rr_sim_interpolated.rds"))


# Disability weights
dw_table <- import(here("data", "dw_table.xlsx"))


# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                     3. SET PARAMETERS                                                        #
################################################################################################################################
# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "cancer", "diab2", "dem", "dep")



################################################################################################################################
################################################################################################################################
#                                                      HIA INPUT DATASET                                                       #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                                4. DISEASES INCIDENCE DATASET                                                 #
################################################################################################################################

# Data preparation : One unique table, aggregated to 10-year age categories
incid_mid_10 <- incid_table_mid %>%
  mutate(age_grp10 = age10_categ(age_grp)) %>%
  group_by(age_grp10, sex) %>%
  summarise(pop_age_grp10 = sum(pop_age_grp, na.rm = TRUE),
            across(ends_with("_incidence"), sum, na.rm = TRUE),
            .groups = "drop")  %>% 
  rename_with(~ str_replace(.x, "_incidence$", "_incidence_mid"), ends_with("_incidence"))

incid_low_10 <- incid_table_low %>%
  mutate(age_grp10 = age10_categ(age_grp)) %>%
  group_by(age_grp10, sex) %>%
  summarise(pop_age_grp10 = sum(pop_age_grp, na.rm = TRUE),
            across(ends_with("_incidence"), sum, na.rm = TRUE),
            .groups = "drop") %>%
  rename_with(~ str_replace(.x, "_incidence$", "_incidence_low"), ends_with("_incidence"))

incid_up_10 <- incid_table_up %>%
  mutate(age_grp10 = age10_categ(age_grp)) %>%
  group_by(age_grp10, sex) %>%
  summarise(pop_age_grp10 = sum(pop_age_grp, na.rm = TRUE),
            across(ends_with("_incidence"), sum, na.rm = TRUE),
            .groups = "drop") %>%
  rename_with(~ str_replace(.x, "_incidence$", "_incidence_up"), ends_with("_incidence"))

incidence_table <- incid_mid_10 %>%
  left_join(incid_low_10 %>% select(-pop_age_grp10), by = c("age_grp10", "sex")) %>%
  left_join(incid_up_10 %>% select(-pop_age_grp10), by = c("age_grp10", "sex")) %>%
  mutate(across(where(is.character) & !matches("age_grp10|sex"), as.numeric))


# Pivot longer
incidence_long_table <- incidence_table %>%
  pivot_longer(
    cols = matches(".*_incidence_(mid|low|up)$"),
    names_to = c("disease", ".value"),
    names_sep = "_incidence_")



##############################################################
#               DISEASE INCIDENCE DISTRIBUTIONS              #
##############################################################
# Generate incidence normal distributions
set.seed(123)
incidence_distrib_table <- incidence_long_table %>%
  rowwise() %>%
  mutate(incidence_distrib = list(
    generate_RR_distrib(mid, low, up, 1000)))


# Sort the incidence distributions in ascending order
incidence_distrib_table <- incidence_distrib_table %>% 
  rowwise() %>% 
  mutate(incidence_distrib = list(sort(unlist(incidence_distrib)))) %>%
  ungroup()

# Separate the incidence_distrib column into multiple columns
incidence_distrib_table <- incidence_distrib_table  %>% 
  unnest_wider(incidence_distrib, names_sep = "_") %>% 
  pivot_longer(
    cols = starts_with("incidence_distrib_"), 
    names_to = "simulation_id",                     # column name for the simulation ID
    values_to = "simulated_incidence")  %>%         # column name for the simulated incidence values
  mutate(simulation_id = as.numeric(str_remove(simulation_id, "incidence_distrib_")))      # simulation ID as a numeric value


# Calculate incidence rate
incidence_distrib_table <- incidence_distrib_table %>%
  mutate(rate = simulated_incidence / pop_age_grp10)



################################################################################################################################
#                                                  5. REDUCTION RISKS DATASET                                                  #
################################################################################################################################

##############################################################
#                REDUCTION RISKS CENTRAL VALUE               #
##############################################################
# Baseline RR
rr_central_baseline <- rr_central_table %>%
  filter(step == baseline_step) %>%
  select(disease, rr2000 = mid)


# Compute reduction risk
reduc_central_table <- rr_central_table %>%
  left_join(rr_central_baseline, by = c("disease")) %>%
  mutate(reduction_risk = (rr2000 - mid) / rr2000)



##############################################################
#                REDUCTION RISKS DISTRIBUTIONS               #
##############################################################
# Baseline RR
rr_baseline <- rr_distrib_table %>%
  filter(step == baseline_step) %>%
  select(disease, simulation_id, rr2000 = rr_interpolated)


# Compute reduction risk
reduc_distrib_table <- rr_distrib_table %>%
  left_join(rr_baseline, by = c("disease", "simulation_id")) %>%
  mutate(reduction_risk = (rr2000 - rr_interpolated) / rr2000)



################################################################################################################################
#                                                6. DISABILITY WEIGHTS DATASET                                                 #
################################################################################################################################

##############################################################
#              DISABILITY WEIGHTS DISTRIBUTIONS              #
##############################################################
# Generate DW normal distributions
set.seed(123)
dw_distrib_table <- dw_table %>%
  rowwise() %>%
  mutate(dw_distrib = list(
    generate_RR_distrib(dw_mid, dw_low, dw_up, 1000)))



# Sort the DW distributions in ascending order
dw_distrib_table <- dw_distrib_table %>% 
  rowwise() %>% 
  mutate(dw_distrib = list(sort(unlist(dw_distrib)))) %>%
  ungroup()


# Separate the dw_distrib column into multiple columns
dw_distrib_table <- dw_distrib_table  %>% 
  unnest_wider(dw_distrib, names_sep = "_") %>% 
  pivot_longer(
    cols = starts_with("dw_distrib_"), 
    names_to = "simulation_id",                     # column name for the simulation ID
    values_to = "simulated_dw"                      # column name for the simulated RR values
  ) %>% 
  mutate(simulation_id = as.numeric(str_remove(simulation_id, "dw_distrib_")))      # simulation ID as a numeric value



################################################################################################################################
#                                                       7. EXPORT DATA                                                         #
################################################################################################################################

# Incidence distribution table
  export(incidence_distrib_table, here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))


# Reduction risks central table
  export(reduc_central_table, here("data_clean", "Diseases", "DRF", "reduction_risk_central_table.xlsx"))
# Reduction risks distribution table
  export(reduc_distrib_table, here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))


# Disability weights distribution table
  export(dw_distrib_table, here("data_clean", "Diseases", "dw_distrib_table.xlsx"))





################################################################################################################################
################################################################################################################################
#                                                      WALKERS DATABASE                                                        #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                           5. CREATION OF SUBSET OF EMP DATASET                                               #
################################################################################################################################
# --------------------------------------
# INDIVIDUAL
# --------------------------------------

# Total walking distance
walk_individual <- emp_walk_ind %>% 
  mutate(nbkm_intermodal_walk = intermodal_walk_time * walk_speed / 60,
         nbkm_tot_walking = nbkm_main_walk + nbkm_intermodal_walk,
         nbkm_tot_walking_jour = nbkm_tot_walking * pond_jour / (pond_indc * 7))


# Create walkers dataset combing walking exposure for each individual
walkers <- walk_dataset(walk_individual, incid_mid_10, insee, morbi_vec, 
                           walk_dist_var = "nbkm_tot_walking", 
                           walk_dist_jour_var = "nbkm_tot_walking_jour", 
                           step_length = step_length, 
                           walk_speed = walk_speed)



# --------------------------------------
# TRIPS
# --------------------------------------

# Total walking distance per trip
walk_trip <- emp_walk_trip  %>% 
  mutate(nbkm_intermodal_walk = intermodal_walk_time * walk_speed / 60,
         nbkm_tot_walking = nbkm_main_walk + nbkm_intermodal_walk,
         nbkm_tot_walking_jour = nbkm_tot_walking * pond_jour / (pond_indc * 7))


# Area type depend on density
walk_trip <- walk_trip %>%
  mutate(
    area_type = case_when(
      densitecom_ori == 1    ~ "urban",
      densitecom_ori == 2    ~ "periurban",
      TRUE                   ~ "rural"),
    area_type = factor(area_type, levels = c("rural", "periurban", "urban")))


# Create walkers dataset combing diseases incidence and walking exposure for each individual
walking_trip <- walk_dataset(walk_trip, incid_mid_10, insee, morbi_vec, 
                           walk_dist_var = "nbkm_tot_walking", 
                           walk_dist_jour_var = "nbkm_tot_walking_jour", 
                           step_length = step_length, 
                           walk_speed = walk_speed)




################################################################################################################################
#                                                   5. DISEASE WALKERS EMP SUBSET                                                      #
################################################################################################################################
# --------------------------------------
# INDIVIDUAL
# --------------------------------------
walkers_long <- walkers %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )


# --------------------------------------
# TRIPS
# --------------------------------------
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

# Total walking distance had those car trips been walked per trip
car_trip <- emp_car_trip  %>% 
  mutate(nbkm_car_jour = nbkm_car * pond_jour / (pond_indc * 7))

# Create drives dataset combing diseases incidence and walking exposure for each individual
car_trip <- walk_dataset(car_trip, incid_mid_10, insee, morbi_vec, 
                           walk_dist_var = "nbkm_car", 
                           walk_dist_jour_var = "nbkm_car_jour", 
                           step_length = step_length, 
                           walk_speed = walk_speed)


# Area type depend on density
car_trip <- car_trip %>%
  mutate(
    area_type = case_when(
      densitecom_ori == 1    ~ "urban",
      densitecom_ori == 2    ~ "periurban",
      TRUE                   ~ "rural"),
    area_type = factor(area_type, levels = c("rural", "periurban", "urban")))


# Associate drive speed
car_trip <- car_trip %>%
  mutate(drive_speed = case_when(
    area_type == c("urban", "periurban") ~ urban_car_speed,
    tuu2017_ori == 8     ~ paris_car_speed,
    TRUE                 ~ rural_car_speed))



################################################################################################################################
#                                                   5. DISEASE EMP SUBSET                                                      #
################################################################################################################################

car_trip_long <- car_trip %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )



################################################################################################################################
#                                                    6. EXPORT EMP SUBSET                                                      #
################################################################################################################################

# Tables of drivers
  export(car_trip, here("data_clean", "EMP_car_trips.xlsx"))
  export(car_trip_long, here("data_clean", "EMP_dis_car_trips.xlsx"))
