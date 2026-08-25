#################################################
########          BEST PRACTICE          ########
########       Commune de résidence      ########
#################################################



################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  purrr,        # Loop
  srvyr,        # Survey
  tidyr,        # Table - Data organization, extraction
  tidyverse,    # Data manipulation and visualization
  ggplot2       # Plotting
)




################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Drivers dataset
emp_walk <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Depression duration distribution table
dep_distrib_table <- import(here("data_clean", "Diseases", "dep_duration_distrib_table.xlsx"))

# Risk reduction distribution table
reduction_risk_distrib_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))

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
dis_vec <- c("mort", "cvd", "cancer", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")



################################################################################################################################
#                                                     4. STEP TARGETS                                                          #
################################################################################################################################
# Mean intermodal walking distance associated to different modes of transport
intermodal_area <- emp_walk_trips %>%
  filter(!is.na(pond_jour)) %>% 
  as_survey_design(ids = ident_ind, 
                   weights = pond_jour) %>% 
  group_by(area_type, mode) %>%
  summarise(intermodal_mean_km = survey_mean(nbkm_intermodal_walk, na.rm = TRUE))  %>% 
  mutate(intermodal_mean_step = intermodal_mean_km / step_length)


# Modal share of each mode of transport by area type (ADEME, September 2025 - "Contribution de la marche et du vélo à la décarbonation et à l'amélioration de la qualité de l'air en France")
intermodal_area <- intermodal_area %>%
  mutate(
    modal_share = case_when(
      area_type == "urban" & mode == "car" ~ 0.27,
      area_type == "urban" & mode == "public_transport" ~ 0.16,
      area_type == "urban" & mode == "bike" ~ 0.13,
      area_type == "urban" & mode == "walk" ~ 0.44,
      area_type == "urban" & mode == "other" ~ 0,
      area_type == "periurban" & mode == "car" ~ 0.64,
      area_type == "periurban" & mode == "public_transport" ~ 0.04,
      area_type == "periurban" & mode == "bike" ~ 0.05,
      area_type == "periurban" & mode == "walk" ~ 0.27,
      area_type == "periurban" & mode == "other" ~ 0,
      area_type == "rural" & mode == "car" ~ 0.74,
      area_type == "rural" & mode == "public_transport" ~ 0.03,
      area_type == "rural" & mode == "bike" ~ 0.04,
      area_type == "rural" & mode == "walk" ~ 0.18,
      area_type == "rural" & mode == "other" ~ 0.01,
      TRUE ~ NA_real_
    )
  )


# Target steps by area_type : weighted mean of intermodal walking distance by mode of transport
target_steps_area <- intermodal_area %>%
  group_by(area_type) %>%
  summarise(target_steps = case_when(
    first(area_type) == "urban"     ~ urban_target + sum(intermodal_mean_step * modal_share, na.rm = TRUE),
    first(area_type) == "periurban" ~ periurban_target + sum(intermodal_mean_step * modal_share, na.rm = TRUE),
    first(area_type) == "rural"     ~ rural_target + sum(intermodal_mean_step * modal_share, na.rm = TRUE),
  ))


# Associate target steps to individuals
urban_target <- target_steps_area %>%
  filter(area_type == "urban") %>%
  pull(target_steps)

periurban_target <- target_steps_area %>%
  filter(area_type == "periurban") %>%
  pull(target_steps)

rural_target <- target_steps_area %>%
  filter(area_type == "rural") %>%
  pull(target_steps)

emp_walk <- emp_walk  %>% 
  mutate(target_area = case_when(
    area_type == "urban"     ~ urban_target,
    area_type == "periurban" ~ periurban_target,
    area_type == "rural"     ~ rural_target
  ))



################################################################################################################################
#                                                   5. AGE DISTRIBUTION                                                        #
################################################################################################################################
# Survey design ponderated by day
main_jour <- emp_walk %>% 
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

emp_target_walk <- emp_walk %>% 
  left_join(distrib_main_walk_EMP2019 %>% 
  select(age_grp10, area_type, rho), by = c("age_grp10", "area_type"))



# Calculate population per age group and area
emp_target_walk <- emp_target_walk  %>% 
  group_by(age_grp10, area_type, disease)  %>% 
  mutate(pop_age_area = sum(pond_indc, na.rm = TRUE))  %>% 
  ungroup()


# Calculate target steps per person for each individual, adjusted by the distribution coefficient and population structure
emp_target_walk <- emp_target_walk %>% 
  mutate(target_ind = rho * target_area * sum(pop_age_area) / sum(pop_age_area * rho))



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
  mutate(step = pmin(12000, round(target_ind / 100) * 100 + baseline_step))  %>% 
  mutate(step_2019 = pmin(12000, round(step_commute/ 100) * 100 + baseline_step))


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
                                       calc_HIA_replicate,
                                       incidence_distrib_table, dep_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                                       dis_vec,
                                       prop_relapse, duration_recovery, vsl,
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
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "Morbidity",
         area_type = "All") %>%
  select(disease, area_type, everything()) 


# Total per disease without area_type stratification
PRACT_burden_by_disease <- PRACT_burden_per_area %>%
  group_by(disease) %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(area_type = "All") %>%
  select(disease, area_type, everything())


# Total for all diseases by area_type
PRACT_burden_by_area <- PRACT_burden_per_area %>%
  group_by(area_type) %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "All") %>%
  select(disease, area_type, everything())


# Total for all diseases across all area types
PRACT_burden_global <- PRACT_burden_per_area %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "All",
         area_type = "All") %>%
  select(disease, area_type, everything())


# Gather results
PRACT_burden <- bind_rows(PRACT_burden_per_area,
                          PRACT_burden_by_disease,
                          PRACT_burden_by_area,
                          PRACT_burden_global,
                          PRACT_burden_morbidity)




###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                          HIA - 2019                                                                                     #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                   5. DATA PREPARATION                                                        #
################################################################################################################################

# Initialization
emp_2019 <- emp_walk  %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000, round(step_commute/ 100) * 100 + 2000))


# EMP Dataset per disease and bound
emp_2019_list <- list()

for (dis in dis_vec) {
  emp_2019_list[[dis]] <- emp_2019 %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                6. HEALTH IMPACT ASSESSMENT                                                   #
################################################################################################################################
# Total of prevented burden of each disease per area type for each simulation 
set.seed(123)
burden_2019_total <- HIA_burden_total(emp_2019_list,
                                      calc_HIA_replicate,
                                      incidence_distrib_table, dep_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                                      dis_vec,
                                      prop_relapse, duration_recovery, vsl,
                                      group = "area_type",
                                      N = 1000)

# Export HIA outcomes of 1000 replications
export(burden_2019_total, here("output", "RDS", "2019", "HIA_area_2019_1000replicate.rds"))




################################################################################################################################
#                         7. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

# Import data
burden_2019_total <- import(here("output", "RDS", "2019", "HIA_area_2019_1000replicate.rds"))


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
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "Morbidity",
         area_type = "All") %>%
  select(disease, area_type, everything())

# Total per disease without area_type stratification
burden_by_disease <- burden_per_area %>%
  group_by(disease) %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(area_type = "All") %>%
  select(disease, area_type, everything())

# Total for all diseases by area_type
burden_by_area <- burden_per_area %>%
  group_by(area_type) %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "All") %>%
  select(disease, area_type, everything())

# Total for all diseases across all area types
burden_global <- burden_per_area %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(disease = "All",
         area_type = "All") %>%
  select(disease, area_type, everything())

# Gather results
burden_2019 <- bind_rows(burden_per_area,
                         burden_by_disease,
                         burden_by_area,
                         burden_global,
                         burden_morbidity)




###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                               COMPARISON 2019 vs BEST PRACTICE                                                                          #
###########################################################################################################################################################################
###########################################################################################################################################################################

# Import data
burden_2019 <- import(here("output", "Tables", "Best practice", "HIA_area_2019_1000replicate.xlsx"))
PRACT_burden <- import(here("output", "Tables", "Best practice", "HIA_best_practice_1000replicate.xlsx"))


# Incremental benefits
PRACT_burden_add <- PRACT_burden %>% 
    mutate(tot_cases = tot_cases - burden_2019[["tot_cases"]],
           tot_cases_low = tot_cases_low - burden_2019[["tot_cases_low"]], 
           tot_cases_up = tot_cases_up - burden_2019[["tot_cases_up"]],
           tot_daly = tot_daly - burden_2019[["tot_daly"]],
           tot_daly_low = tot_daly_low - burden_2019[["tot_daly_low"]],
           tot_daly_up = tot_daly_up - burden_2019[["tot_daly_up"]],
           tot_medic_costs = tot_medic_costs - burden_2019[["tot_medic_costs"]],
           tot_medic_costs_low = tot_medic_costs_low - burden_2019[["tot_medic_costs_low"]],
           tot_medic_costs_up = tot_medic_costs_up - burden_2019[["tot_medic_costs_up"]],
           tot_soc_costs = tot_soc_costs - burden_2019[["tot_soc_costs"]],
           tot_soc_costs_low = tot_soc_costs_low - burden_2019[["tot_soc_costs_low"]],
           tot_soc_costs_up = tot_soc_costs_up - burden_2019[["tot_soc_costs_up"]])  %>% 
    select(disease, area_type, tot_cases, tot_cases_low, tot_cases_up, tot_daly, tot_daly_low, tot_daly_up,
           tot_medic_costs, tot_medic_costs_low, tot_medic_costs_up, tot_soc_costs, tot_soc_costs_low, tot_soc_costs_up)


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
      filter(!area_type %in% c("All")) %>%
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
      filter(!area_type %in% c("All")) %>%
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
    values = colors_area,
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

# Description : how many people below the targets and age/sex distribution below targets by area type

# Proportion of French adult below targets by area type
emp_below_target <- PRACT_list[["mort"]]  %>% 
  mutate(below_target = step_2019 < step)  %>% 
  group_by(area_type, sex, age_grp10)  %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc)

# Proportion of French adult below targets by area type
prop_below_target <- emp_below_target  %>% 
  group_by(area_type)  %>% 
  summarise(perc = 100 * survey_mean(below_target, na.rm = TRUE, vartype = "ci"))
  


# Proportion of French adult below targets by area type depending on age
prop_below_target_age <- emp_below_target  %>% 
  group_by(area_type, age_grp10)  %>% 
  summarise(perc = 100 * survey_mean(below_target, na.rm = TRUE, vartype = "ci"))

plot_prop_below_target_age <- ggplot(prop_below_target_age, aes(x = age_grp10, y = perc,
                                            ymin = perc_low, ymax = perc_upp, fill = area_type)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_area) +
  ylab ("Proportion of walkers below the best practice targets in the past day (%)") +
  xlab("Age group") +
  theme_minimal()
plot_prop_below_target_age


# Proportion of French adult below targets by area type depending on sex
prop_below_target_sex <- emp_below_target  %>% 
  group_by(area_type, sex)  %>% 
  summarise(perc = 100 * survey_mean(below_target, na.rm = TRUE, vartype = "ci"))

plot_prop_below_target_sex <- ggplot(prop_below_target_sex, aes(x = area_type, y = perc,
                                            ymin = perc_low, ymax = perc_upp, fill = sex)) +
  geom_col(width = 0.7, position = position_dodge2(0.4))+
  geom_errorbar(position = position_dodge(0.7), width = 0.25) +
  scale_fill_manual(values = colors_sex) +
  ylab ("Proportion of walkers below the best practice targets in the past day (%)") +
  xlab("Area type") +
  theme_minimal()
plot_prop_below_target_sex



################################################################################################################################
#                                                       9. EXPORT DATA                                                         #
################################################################################################################################

# Plot
ggsave(here("output", "Plots", "Best practice", "best_practice_cases_prev.png"), plot = plot_PRACT_cases_prev)
ggsave(here("output", "Plots", "Best practice", "prop_below_target_age.png"), plot = plot_prop_below_target_age)
ggsave(here("output", "Plots", "Best practice", "prop_below_target_sex.png"), plot = plot_prop_below_target_sex)

# Tables
export(PRACT_burden, here("output", "Tables", "Best practice", "HIA_best_practice_1000replicate.xlsx"))
export(burden_2019, here("output", "Tables", "Best practice", "HIA_area_2019_1000replicate.xlsx"))
export(PRACT_burden_add, here("output", "Tables", "Best practice", "HIA_best_practice_added_1000replicate.xlsx"))
export(prop_below_target, here("output", "Tables", "Best practice", "prop_below_target.xlsx"))
