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
dis_vec = c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")
morbi_vec = c("cvd", "bc", "cancer", "diab2", "dem", "dep")

################################################################################################################################
#                                                4. DISEASES DATA PREPARATION                                                  #
################################################################################################################################
# For simplification of coding, rename dep_prevalence in dep_incidence but these are prevalences
diseases_5 <-  diseases %>% 
  rename (dep_incidence = dep_prevalence)

# Age categories : 10 years
diseases_5 <- diseases_5 %>%
  mutate(
    age_grp10 = case_when(
      age_grp.x %in% c("20-24") ~ "20-24",
      age_grp.x %in% c("25-29","30-34") ~ "25-34",
      age_grp.x %in% c("35-39","40-44") ~ "35-44",
      age_grp.x %in% c("45-49","50-54") ~ "45-54",
      age_grp.x %in% c("55-59","60-64") ~ "55-64",
      age_grp.x %in% c("65-69","70-74") ~ "65-75",
      age_grp.x %in% c("75-79","80-84","85-89") ~ "75-89",
      TRUE ~ NA_character_
    )
  )


diseases_10 <- diseases_5 %>%
  group_by(sex, age_grp10) %>%
  summarise(
    pop_age_sex = sum(pop_age_grp, na.rm = TRUE),
    across(
      .cols = paste0(morbi_vec, "_incidence"),
      .fns  = ~ sum(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  )



################################################################################################################################
#                                 5. CREATION OF SUBSET OF EMP DATASET WITH ONLY VARIABLES NEEDED                              #
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
    mdisttot_fin1
  )


# Re-write nbkm_walking, easier to be used in functions
emp_subset <- emp_subset %>%
  rename (nbkm_walking = nbkm_walking_lower) %>% 

# Re-write sexe as female and male and convert as factors
  mutate (sexe = as.character(sexe)) %>%                                 # Conversion in character for function to work well
  mutate (sexe = fct_recode(sexe, "Male" = "1", "Female" = "2")) %>%     # Replacing
  rename (sex = sexe)


# Daily steps
emp_subset <- emp_subset %>% 
  mutate(step_commute = nbkm_walking/step_length)

# Day time spent walking (min)
emp_subset <- emp_subset %>% 
  mutate(day_time = nbkm_walking*60/walk_speed)

# Week time spent walking (min)
emp_subset <- emp_subset %>% 
  mutate(week_time = 7*nbkm_walking*60/walk_speed)


# Create age categories
emp_subset <- emp_subset %>% 
  mutate(
    age_grp10 = case_when(
      age >= 20 & age <= 24 ~ "20-24",
      age >= 25 & age <= 34 ~ "25-34",
      age >= 35 & age <= 44 ~ "35-44",
      age >= 45 & age <= 54 ~ "45-54",
      age >= 55 & age <= 64 ~ "55-64",
      age >= 65 & age <= 75 ~ "65-75",
      age >= 75 & age <= 89 ~ "75-89",
      TRUE ~ NA_character_
    )
  )



# Add population counts per sex and age group
emp_subset <- emp_subset %>% 
  left_join(
      diseases_10 %>% 
        select(pop_age_sex, sex, age_grp10),     # Matching columns
        by = c("sex", "age_grp10")               # Fill the variables depending on sex and age group
    )


# Add diseases incidences / prevalences
emp_subset <- emp_subset %>% 
  left_join(
    diseases_10 %>% 
      select(cvd_incidence, bc_incidence, cancer_incidence, diab2_incidence, dem_incidence, dep_incidence, sex, age_grp10),    # Matching columns
      by = c("sex", "age_grp10")                                                                                               # Fill the variables depending on sex and age group
  ) 


# Add mortality rates
insee <- insee %>% 
  rename(sex = sexe)

emp_subset <- emp_subset %>% 
  left_join(
    insee %>% select(MR, sex, age),       # Matching columns
    by = c("sex", "age")                  # Fill the variables depending on sexe and age
  ) %>% 
  rename(mort_rate = MR)


# Calculate death incidence
emp_subset <-  emp_subset %>% 
  mutate(mort_incidence = mort_rate * pop_age_sex)



# Calculate disease incidence rates
for (dis in morbi_vec){
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
emp_subset <- emp_subset %>%
  mutate(
    life_exp = if_else(
      sex == "Female",
      85.99324,
      79.59503
    )
  )


# Add the years of life remaining, potentially affected by diseases or premature death
emp_subset <- emp_subset %>%
  mutate(
    years_remaining = pmax(life_exp - age, 0)
  )


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
         -day_time)


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

# Day time spent walking if those car distances were walked (min)
emp_drivers <- emp_drivers %>% 
  mutate(week_time_shift = 7*mdisttot_fin1*60 / walk_speed)


# Daily steps if those car distances were walked
emp_drivers <- emp_drivers %>% 
  mutate(day_step_commute_shift = mdisttot_fin1/step_length)



################################################################################################################################
#                                                    5. EXPORT EMP SUBSET                                                      #
################################################################################################################################

export(emp_drivers, here("data_clean", "EMP_drivers.xlsx"))











