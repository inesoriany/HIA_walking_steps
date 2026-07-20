###############################################
############ EXTRACT EMP DATABASE #############
###############################################

# This code imports EMP data files from csv and puts together datasets at the individual and trip level, with information on:
# - individual characteristics
# - total number of kilometers by walk, car



################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data manipulation
  tidyverse    # Data management
  )


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Individual data
ind <- import(here("data", "emp_2019_donnees_individuelles_anonymisees_novembre2024", "k_individu_public_V3.csv" ))

# Individual characteristics
ind_kish <- import(here("data", "emp_2019_donnees_individuelles_anonymisees_novembre2024", "tcm_ind_kish_public_V3.csv" ))

# Household characteristics 
household <- import(here("data", "emp_2019_donnees_individuelles_anonymisees_novembre2024", "tcm_men_public_V3.csv"))

# Trip data
trip <- import(here("data", "emp_2019_donnees_individuelles_anonymisees_novembre2024", "5. k_deploc_public_V4.csv" ))




################################################################################################################################
#                                                     3. SELECT DATA                                                           #
################################################################################################################################

# Individual ID and ponderations
ind <- ind %>% 
  rename_with(tolower) %>% 
  select(ident_ind, 
         pond_indc) %>% 
  mutate(ident_ind = as.character(ident_ind))
  

# Individual characteristics
ind_kish <- ind_kish %>% 
  rename_with(tolower) %>% 
  select(ident_men,
         ident_ind,
         sexe, 
         age) %>% 
  mutate(ident_men = as.character(ident_men),
         ident_ind = as.character(ident_ind))


# Household characteristics
household <- household  %>% 
  rename_with(tolower) %>% 
  select(ident_men,
        densitecom_res)  %>% 
  mutate(ident_men = as.character(ident_men))


# Trips
trip <- trip %>% 
  rename_with(tolower) %>% 
  select(ident_ind, 
         ident_dep,
         pond_jour,
         mdisttot_fin,                               # Trip length
         mtempsmap,                                  # Walking duration
         mtp,                                        # Main mean of transportation
         densitecom_ori,                             # Departure commune density  
         densitecom_des,                             # Arrival commune density
         tuu2017_ori) %>%
  mutate(ident_ind = as.character(ident_ind),
         ident_dep = as.character(ident_dep))




##############################################################
#                          WALKING                           #
##############################################################
# --------------------------------------
# INDIVIDUAL
# --------------------------------------
# Walking levels
walk_ind <- trip %>% 
  mutate(intermodal_walk_time = if_else(is.na(mtempsmap), 0, mtempsmap)) %>% 
  mutate(main_walk = mtp == 1.1,
         nbkm_main_walk = if_else(main_walk, mdisttot_fin, 0)) %>%
  group_by(ident_ind, pond_jour) %>%
  summarise(intermodal_walk_time = sum(intermodal_walk_time, na.rm = TRUE),
            nbkm_main_walk = sum(nbkm_main_walk, na.rm = TRUE),
            .groups = "drop") %>% 
  select(ident_ind,
         pond_jour,
         intermodal_walk_time,
         nbkm_main_walk) %>% 
  
# Add individual characteristics
  left_join(ind, by = "ident_ind") %>%
  left_join(ind_kish, by = "ident_ind")  %>% 

# Add household characteristics
  left_join(household, by = "ident_men")




# --------------------------------------
# TRIPS
# --------------------------------------
walk_trip <- ind %>%                                  
  left_join(trip, by ="ident_ind", relationship = "many-to-many") %>% 
  mutate(intermodal_walk_time = if_else(is.na(mtempsmap), 0, mtempsmap)) %>% 
  mutate(main_walk = mtp == 1.1,
         nbkm_main_walk = if_else(main_walk, mdisttot_fin, 0)) %>% 
  
# Add individual characteristics
  left_join(ind_kish, by = "ident_ind")



##############################################################
#                          DRIVING                           #
##############################################################
# --------------------------------------
# TRIPS
# --------------------------------------
car_trip <- trip %>% 
  mutate(main_car = mtp %in% c(3.1, 3.2, 3.3, 3.4)) %>% 
  filter(main_car) %>% 
  mutate(nbkm_car = if_else(main_car, mdisttot_fin, 0)) %>% 
  group_by(ident_ind, ident_dep, pond_jour) %>% 
  select(ident_ind,
         ident_dep,
         pond_jour,
         nbkm_car,
         densitecom_ori,
         densitecom_des,
         tuu2017_ori) %>% 
  
# Add individual characteristics 
  left_join(ind, by = "ident_ind") %>% 
  left_join(ind_kish, by = "ident_ind")



################################################################################################################################
#                                                  4. EXPORT EMP DATASETS                                                      #
################################################################################################################################
# Walking at the individual level
  export(walk_ind, here("data", "emp_dataset_walk_individual.xlsx"))

# Walking trips
  export(walk_trip, here("data", "emp_dataset_walk_trip.xlsx"))

# Car trips
  export(car_trip, here("data", "emp_dataset_car_trip.xlsx"))
