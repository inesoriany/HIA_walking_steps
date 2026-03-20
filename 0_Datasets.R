###############################################
############ CREATING DATASETS ################
###############################################

# Files needed :
  # emp_dataset_km_bike_and_car_and_walk_individual.csv : Walking levels (EMP)
  # out_merged.csv : incidence data
  # INSEE_2019.RDS : population data


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

# EMP 2019 : distances for bike, cars, walking
emp <- import(here("data", "emp_dataset_km_bike_and_car_and_walk_individual.csv")) 

# Diseases and mortality data
diseases <- import(here("data", "GBD_diseases.xlsx"), sheet = "Central")

# INSEE data
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
#                                 5. CREATION OF SUBSET OF EMP DATASET WITH ONLY VARIABLES NEEDED                              #
################################################################################################################################

# Creation of subset uniting variables of EMP and calculations
emp_subset <- emp %>% 
  mutate(ID = row_number()) %>%                   # ID number for every individual
  select(
    ID,
    ident_ind,
    sexe,
    age,
    quartile_rev,
    DENSITECOM_RES,
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
  mutate(age_grp = age_grp(age),
         age_grp10 = age_grp_10(age)) %>%
  filter(age >= 20 & age <90)                        # Only keep ages 20-89 years 



# Add population counts per sex and age group
emp_subset <- emp_subset %>% 
  left_join(
    diseases %>% 
      select(pop_age_grp, sex, age_grp),     # Matching columns
    by = c("sex", "age_grp")               # Fill the variables depending on sex and age group
  )

emp_subset <- emp_subset %>% 
  left_join(
      diseases_10 %>% 
        select(pop_age_grp10, sex, age_grp10),     # Matching columns
        by = c("sex", "age_grp10")               # Fill the variables depending on sex and age group
    )


# Add diseases incidences / prevalences
emp_subset <- emp_subset %>% 
  left_join(
    diseases %>% 
      select(cvd_incidence, bc_incidence, cancer_incidence, diab2_incidence, dem_incidence, dep_incidence, sex, age_grp),    # Matching columns
      by = c("sex", "age_grp")                                                                                               # Fill the variables depending on sex and age group
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
  mutate(mort_incidence = mort_rate * pop_age_grp)



# Calculate disease incidence rates
for (dis in morbi_vec){
  emp_subset <- emp_subset %>%
    mutate(
      !!paste0(dis, "_rate") := if_else(
        !is.na(pop_age_grp),
        .data[[paste0(dis, "_incidence")]] / pop_age_grp,
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



# Area type depend on density
emp_subset <- emp_subset %>%
  mutate(
    area_type = case_when(
      DENSITECOM_RES == 1    ~ "urban",
      DENSITECOM_RES == 2    ~ "periurban",
      TRUE                   ~ "rural"
    ),
    area_type = factor(area_type, levels = c("rural", "periurban", "urban"))
  )



################################################################################################################################
#                                                   5. DISEASE EMP SUBSET                                                      #
################################################################################################################################

emp_long <- emp_subset %>%
  pivot_longer(
    cols = matches("(_rate|_incidence)$"),   # all rate and incidence columns
    names_to = c("disease", ".value"),       # disease name + output column name
    names_pattern = "(.*)_(rate|incidence)"  # regex: capture disease + type
  )


################################################################################################################################
#                                                       6. EXPORT DATA                                                         #
################################################################################################################################

# Tables of walkers
  export(emp_subset, here("data_clean", "EMP_walkers.xlsx"))
  export(emp_long, here("data_clean", "EMP_dis_walkers.xlsx"))




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

