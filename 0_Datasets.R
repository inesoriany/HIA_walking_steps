###############################################
############ CREATING DATASETS ################
###############################################

# Files needed :
  # emp_dataset_km_bike_and_car_and_walk_individual.csv
  # out_merged.csv
  # INSEE_2019.RDS


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
  epikit,       # Age categories creation
  janitor,      # De-duplication
  survey        # Survey management
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# EMP 2019 : distances for bike, cars, walking
emp <- import(here("data", "emp_dataset_km_bike_and_car_and_walk_individual.csv")) 

# Diseases and mortality data
diseases <- import(here("data", "out_merged.csv"))

# INSEE data
insee <- import(here("data", "INSEE_2019.RDS"))


################################################################################################################################
#                                                     3. SET PARAMETERS                                                        #
################################################################################################################################
# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "diab2", "dem")


################################################################################################################################
#                                  4. CREATION OF SUBSET OF EMP SUBSET WITH ONLY VARIABLES NEEDED                              #
################################################################################################################################

# Creation of subset uniting variables of EMP and calculations
emp_subset <- emp %>% 
  select(
    ident_ind,
    sexe,
    age,
    quartile_rev,
    tuu2017_res,
    pond_indc,
    pond_jour,
    nbkm_walking_lower,
    nbkm_walking_upper,
    mdisttot_fin1
  )


# Re-write nbkm_walking, easier to be used in functions
emp_subset <- emp_subset %>%
  rename(nbkm_walking = nbkm_walking_lower) %>% 

# Re-write sexe as female and male and convert as factors
  mutate(sexe = as.character(sexe)) %>%                                 # Conversion in character for function to work well
  mutate(sexe = fct_recode(sexe, "Male" = "1", "Female" = "2"))         # Replacing


# Daily steps
emp_subset <- emp_subset %>% 
  mutate(step = nbkm_walking/step_length)

# Day time spent walking (min)
emp_subset <- emp_subset %>% 
  mutate(day_time = nbkm_walking*60/walk_speed) %>% 
  mutate(day_time_upper = nbkm_walking_upper*60/walk_speed)


# Create age categories
emp_subset <- emp_subset %>% 
  mutate(
    age_grp.x = age_categories(
      age,
      lower = 0,
      upper = 110,
      by = 5)
    ) %>% 
  mutate(age_grp.x = as.character(age_grp.x))


# Add population counts per sex
diseases <- diseases %>% 
  rename(sexe = sex) %>% 
  mutate(sexe = fct_recode(sexe, "Male" = "male", "Female" = "female"))         # Replacing

emp_subset <- emp_subset %>% 
  left_join(
      diseases %>% select(pop_age_grp, sexe, age_grp.x),    # Matching columns
      by = c("sexe", "age_grp.x")                           # Fill the variables depending on sex and age group
    ) %>% 
  rename(pop_age_sex = pop_age_grp)


# Add diseases incidences
emp_subset <- emp_subset %>% 
  left_join(
    diseases %>% select(cc_incidence, dem_incidence, bc_incidence, cvd_incidence, diab2_incidence, sexe, age_grp.x),    # Matching columns
    by = c("sexe", "age_grp.x")                                                                                         # Fill the variables depending on sex and age group
  ) 

# Add mortality rates
emp_subset <- emp_subset %>% 
  left_join(
    insee %>% select(MR, sexe, age),       # Matching columns
    by = c("sexe", "age")                  # Fill the variables depending on sexe and age
  ) %>% 
  rename(mort_rate = MR)


# Calculate death incidence
emp_subset <-  emp_subset %>% 
  mutate(mort_incidence = mort_rate * pop_age_sex)


# Calculate incidence rates
for (dis in dis_vec){
  emp_subset <- emp_subset %>%
    mutate(
      !!paste0(dis, "_rate") := if_else(
        !is.na(pop_age_sex),
        .data[[paste0(dis, "_incidence")]] / pop_age_sex,
        NA_real_
      )
    )
}


# Add life-expectancy for each sex
for (i in 1:nrow(emp_subset)) {
  emp_subset[i, "life_exp"] <- ifelse (
    emp_subset$sexe[i]== "Female", 
    85.99324,                            # Life expectancy for women = 85.99324
    79.59503                             # Life expectancy for men = 79.59503
    )  
  }

# Add the years of life remaining, potentially affected by diseases or premature death
for (i in 1:nrow(emp_subset)) {
  emp_subset[i,"years_remaining"] <- ifelse(
    emp_subset[i, "life_exp"]-emp_subset[i, "age"] >=0,
    emp_subset[i, "life_exp"]-emp_subset[i, "age"],
    0                                                   # No negatives for individuals above life expectancy
    )
}

# Only keep ages 20-89 years 
emp_subset <-  emp_subset %>% 
  filter(age >= 20 & age <90)



# Area type
emp_subset <- emp_subset %>%
  mutate(
    area_type = case_when(
      tuu2017_res %in% 2:4 ~ "semi_urban",
      tuu2017_res %in% 5:8 ~ "urban",
      TRUE                ~ "rural"
    ),
    area_type = factor(area_type, levels = c("rural", "semi_urban", "urban"))
  )



################################################################################################################################
#                                                    5. EXPORT EMP SUBSET                                                      #
################################################################################################################################

export(emp_subset, here("data_clean", "EMP_walkers.xlsx"))






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
         -nbkm_walking_upper,
         -day_time,
         -day_time_upper)


# Associate drive speed
emp_drivers <- emp_drivers %>%
  mutate(drive_speed = case_when(
    tuu2017_res %in% 2:7 ~ urban_car_speed,
    tuu2017_res == 8     ~ paris_car_speed,
    TRUE                 ~ rural_car_speed
    )
  )


# Day time spent walking if those car distances were walked (min)
emp_drivers <- emp_drivers %>% 
  mutate(day_time_shift = mdisttot_fin1*60 / walk_speed)


# Daily steps if those car distances were walked
emp_drivers <- emp_drivers %>% 
  mutate(day_step_shift = mdisttot_fin1/step_length)



################################################################################################################################
#                                                    5. EXPORT EMP SUBSET                                                      #
################################################################################################################################

export(emp_drivers, here("data_clean", "EMP_drivers.xlsx"))











