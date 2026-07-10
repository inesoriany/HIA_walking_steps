#################################################
########          BEST PRACTICE          ########
########        Commune de départ        ########
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
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Import drivers dataset
emp_walk_trip <- import(here("data_clean", "EMP_dis_walking_trips.xlsx"))


# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Risk reduction distribution table
reduction_risk_distrib_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))

# Disability weights distribution table
dw_distrib_table <- import(here("data_clean", "Diseases", "dw_distrib_table.xlsx"))



# Import functions
source(here("0_Functions.R"))




################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec <- c("mort", "cvd", "cancer", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")




################################################################################################################################
#                                                   4. AGE DISTRIBUTION                                                        #
################################################################################################################################
# Survey design ponderated by day
main_jour <- emp_walk_trip %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



# Walking pyramid : Age distribution of walking volume for each territory (in steps)
distrib_main_walk_EMP2019 <- main_jour %>% 
  group_by(age_grp10, area_type) %>% 
  summarise(mean_ind = survey_mean(step_main, na.rm = TRUE))


# Distribution coefficient by area type
distrib_main_walk_EMP2019 <- distrib_main_walk_EMP2019 %>% 
  group_by(area_type) %>% 
  mutate(rho = mean_ind / mean_ind[age_grp10 == "20-29"]) %>% 
  ungroup()

emp_target_walk <- emp_walk_trip %>% 
  left_join(distrib_main_walk_EMP2019 %>% 
  select(age_grp10, area_type, rho), by = c("age_grp10", "area_type"))



# Calculate population per age group and area
emp_target_walk <- emp_target_walk  %>% 
  group_by(age_grp10, area_type, disease)  %>% 
  mutate(pop_age_area = sum(pond_indc, na.rm = TRUE))  %>% 
  ungroup()


# Calculate target to reach for each individual depending on their area type and age group
#(steps per person by age group and area type)
emp_target_walk <- emp_target_walk %>% 
  group_by(area_type) %>%
  mutate(
    target_area = case_when(
      area_type == "urban" ~ urban_target,
      area_type == "periurban" ~ periurban_target,
      area_type == "rural" ~ rural_target),

    
    # target steps per person for each individual, adjusted by the distribution coefficient and population structure
    target_step = rho * target_area * sum(pop_age_area) / sum(pop_age_area * rho)) %>%
  
  ungroup()


###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                       HIA - BEST PRACTICE                                                                               #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                   5. DATA PREPARATION                                                        #
################################################################################################################################

# Initialization
emp_target_walk <- emp_target_walk %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000, round(target_step / 100) * 100 + 2000))  %>% 
  filter(mtp ==1.1)             # Exclusively walking


# EMP Dataset per disease
PRACT_list <- list() 

for (dis in dis_vec) {
  PRACT_list[[dis]] <- emp_target_walk %>% 
    filter(disease == dis)
}



###############################################################################################################################
#                                              6. HEALTH IMPACT ASSESSMENT                                                    #
###############################################################################################################################

# Total of prevented burden of each disease per area type for each simulation 
set.seed(123)
PRACT_burden_total <- HIA_burden_total(PRACT_list,
                                       incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                                       dis_vec,
                                       vsl,
                                       group = "area_type",
                                       N = 1000)


# Export : Table of HIA outcomes per simulation
export(PRACT_burden_total, here("output", "RDS", "Best practice", "HIA_best_practice_1000replicate.rds"))





################################################################################################################################
#                                   7. HIA OUTCOMES : Total of prevented cases, DALY, costs                                    #
################################################################################################################################

# Import data
PRACT_burden_total <- import(here("output", "RDS", "Best practice", "HIA_best_practice_1000replicate.rds"))


##############################################################
#                        PER DISEASE                         #
##############################################################

# --------------------------------------
# MONTE-CARLO
# --------------------------------------

# Prevented diseases per area type
set.seed(123)
PRACT_burden_per_area <- HIA_burden_IC(PRACT_burden_total, dis_vec, outcome_vec, calc_replicate_IC)


# Total for morbidity
PRACT_burden_morbidity <- PRACT_burden_per_area %>%
  filter(disease != "mort") %>% 
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "Morbidity") %>%
  select(disease, everything()) 
  
  
# Total for all diseases
PRACT_burden_global <- PRACT_burden_per_area %>%
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "All") %>%
  select(disease, everything()) 
  

# Gather results
PRACT_burden <- bind_rows(PRACT_burden_per_area, PRACT_burden_morbidity, PRACT_burden_global)




###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                          HIA - 2019                                                                                     #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                   5. DATA PREPARATION                                                        #
################################################################################################################################

# Initialization
walk_2019 <- emp_walk_trip  %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000, round(step_main/ 100) * 100 + 2000))  %>% 
  filter(mtp == 1.1)


# EMP Dataset per disease and bound
walk_trips_2019_list <- list()

for (dis in dis_vec) {
  walk_trips_2019_list[[dis]] <- walk_2019 %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                6. HEALTH IMPACT ASSESSMENT                                                   #
################################################################################################################################
# Total of prevented burden of each disease per area type for each simulation 
set.seed(123)
burden_2019_total <- HIA_burden_total(walk_trips_2019_list,
                                      incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                                      dis_vec,
                                      vsl,
                                      group = "area_type",
                                      N = 1000)

# Export HIA outcomes of 1000 replications
export(burden_2019_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_area_2019_1000replicate.rds"))




################################################################################################################################
#                         7. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

# Import data
burden_2019_total <- import(here("output", "RDS", "2019", "Monte Carlo", "HIA_area_2019_1000replicate.rds"))


##############################################################
#                       PER DISEASES                         #
##############################################################

# --------------------------------------
# MONTE-CARLO
# --------------------------------------

# Prevented diseases per area type
set.seed(123)
burden_per_area <- HIA_burden_IC(burden_2019_total, dis_vec, outcome_vec, calc_replicate_IC)


# Total for morbidity
burden_morbidity <- burden_per_area %>%
  filter(disease != "mort") %>% 
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "Morbidity") %>%
  select(disease, everything()) 
  
  
# Total for all diseases
burden_global <- burden_per_area %>%
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "All") %>%
  select(disease, everything()) 
  

# Gather results
burden_2019 <- bind_rows(burden_per_area, burden_morbidity, burden_global)




###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                               COMPARISON 2019 vs BEST PRACTICE                                                                          #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                     8. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Cases prevented by disease and area type
plot_PRACT_cases_prev <-
  ggplot() +

  # 2019 baseline
  geom_col(data = burden_2019 %>%
      filter(!disease %in% c("All", "Morbidity")) %>%
      mutate(disease = factor(disease, levels = c("mort", "cancer", "cvd", "diab2", "dem", "dep")),
             area_type = factor(area_type, levels = c("urban", "periurban", "rural"))),
      aes(x = disease,
          y = tot_cases,
          fill = area_type,
          alpha = "2019 baseline"),
          position = position_dodge(width = 0.8),
          width = 0.7) +

  geom_errorbar(data = burden_2019 %>%
      filter(!disease %in% c("All", "Morbidity")) %>%
      mutate(disease = factor(disease, levels = c("mort", "cancer", "cvd", "diab2", "dem", "dep")),
             area_type = factor(area_type, levels = c("urban", "periurban", "rural"))),
    aes(x = disease,
        ymin = tot_cases_low,
        ymax = tot_cases_up,
        group = area_type,
        alpha = "2019 baseline"),
        position = position_dodge(width = 0.8),
        width = 0.2) +


  # Best practice
  geom_col(data = PRACT_burden %>%
      filter(!disease %in% c("All", "Morbidity")) %>%
      mutate(disease = factor(disease, levels = c("mort", "cancer", "cvd", "diab2", "dem", "dep")),
             area_type = factor(area_type, levels = c("urban", "periurban", "rural"))),
    aes(x = disease,
        y = tot_cases,
        fill = area_type,
        alpha = "Best practice scenario"),
        position = position_dodge(width = 0.8),
        width = 0.7) +

  geom_errorbar(data = PRACT_burden %>%
      filter(!disease %in% c("All", "Morbidity")) %>%
      mutate(disease = factor(disease, levels = c("mort", "cancer", "cvd", "diab2", "dem", "dep")),
             area_type = factor(area_type, levels = c("urban", "periurban", "rural"))),
    aes(x = disease,
        ymin = tot_cases_low,
        ymax = tot_cases_up,
        group = area_type,
        alpha = "Best practice scenario"),
        position = position_dodge(width = 0.8),
        width = 0.2) +

  scale_fill_manual(
    values = c(urban = "#1b9e77",
               periurban = "#d95f02",
               rural = "#7570b3"),
    name = "Area type") +

  scale_alpha_manual(
    name = "Scenario",
    values = c("2019 baseline" = 1,
               "Best practice scenario" = 0.4)) +

  scale_x_discrete(labels = names_disease) +

  labs(x = "Disease",
       y = "Cases prevented") +

  theme_minimal()


plot_PRACT_cases_prev 

################################################################################################################################
#                                                       8. DESCRIPTION                                                         #
################################################################################################################################

# Description : how many people below the tragets and age/sex distribution below targets by area type



################################################################################################################################
#                                                       9. EXPORT DATA                                                         #
################################################################################################################################

# Plot
ggsave(here("output", "Plots", "Best practice", "best_practice_cases_prev.png"), plot = plot_PRACT_cases_prev)

# Tables
export(PRACT_burden, here("output", "Tables", "Best practice", "HIA_best_practice_1000replicate.xlsx"))
export(burden_2019, here("output", "Tables", "Best practice", "HIA_area_2019_1000replicate.xlsx"))
