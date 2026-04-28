#################################################
#############     BEST STUDENTS     #############
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

# Walking trips dataset
emp_walk_trip <- import(here("data_clean", "EMP_dis_walking_trips.xlsx"))


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
#                                                4. MAIN WALKING TRIPS DATASET                                                 #
################################################################################################################################
# Only exclusively walking trips (no intermodal walk)
  emp_main_walk_trip <- emp_walk_trip %>% 
    filter(intermodal_walk_time == 0)

# Add population by age group and area type
emp_main_walk_trip <- emp_main_walk_trip %>%
  left_join(emp_main_walk_trip %>%
      select(ident_ind, age_grp10, area_type, pond_indc) %>%
      distinct() %>%
      group_by(age_grp10, area_type) %>%
      summarise(pop_age_area = sum(pond_indc, na.rm = TRUE),
      .groups = "drop"), by = c("age_grp10", "area_type"))


# Sum steps by individual, area type, and disease
steps_by_individual_area_disease <- emp_main_walk_trip %>% 
  group_by(ident_ind, area_type, disease) %>% 
  summarise(
    total_steps = sum(step_commute, na.rm = TRUE),
    n_trips = n_distinct(ident_dep),
    pond_indc = first(pond_indc),
    .groups = "drop"
  )

# Merge aggregated data back to original trips dataset
emp_main_walk_trip <- emp_main_walk_trip %>% 
  left_join(
    steps_by_individual_area_disease %>% 
      select(ident_ind, area_type, disease, total_steps, n_trips),
    by = c("ident_ind", "area_type", "disease")
  )


################################################################################################################################
#                                                   5. AGE DISTRIBUTION                                                        #
################################################################################################################################
# Survey design ponderated by day
main_jour <- emp_main_walk_trip %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



# Walking pyramid : Age distribution of walking volume for each territory (in steps)
distrib_main_walk_EMP2019 <- main_jour %>% 
  group_by(age_grp10, area_type) %>% 
  summarise(mean_ind = survey_mean(step_commute, na.rm = TRUE))


# Distribution coefficient by area type
distrib_main_walk_EMP2019 <- distrib_main_walk_EMP2019 %>% 
  group_by(area_type) %>% 
  mutate(rho = mean_ind / mean_ind[age_grp10 == "20-29"]) %>% 
  ungroup()

emp_main_walk_trip <- emp_main_walk_trip %>% 
  left_join(distrib_main_walk_EMP2019 %>% 
  select(age_grp10, rho), by = "age_grp10")


# Calculate target to reach for each individual depending on their area type and age group
emp_main_walk_trip <- emp_main_walk_trip %>% 
  mutate(target_step = rho * case_when(
    area_type == "urban" ~ urban_target,
    area_type == "periurban" ~ periurban_target,
    area_type == "rural" ~ rural_target,
    TRUE ~ NA_real_
  ) / sum(rho*pop_age_area))


# Filter individuals below targets ? and do HIA (to know the added value)
  ################## See how to integrate that later ----
  # Check if individuals meet the targets
  below_targets <- emp_main_walk_trip %>% 
    mutate(
      target = case_when(
        area_type == "urban" ~ urban_target,
        area_type == "periurban" ~ periurban_target,
        area_type == "rural" ~ rural_target,
        TRUE ~ NA_real_
      ),
      meets_target = total_steps >= target
    ) %>% 
    filter(!meets_target)

  cat("\n=== PERSONNES EN DESSOUS DES CIBLES ===\n")
  cat("Total (pondéré) :", sum(below_targets$pond_indc, na.rm = TRUE), "personnes\n")
  cat("Total (non pondéré) :", nrow(below_targets), "personnes\n\n")

  by_area <- below_targets %>% 
    group_by(area_type) %>% 
    summarise(
      n_below_weighted = sum(pond_indc, na.rm = TRUE),
      n_below_unweighted = n(),
      mean_steps = mean(total_steps, na.rm = TRUE),
      .groups = "drop"
    )
  print(by_area)

# data preparation for HIA