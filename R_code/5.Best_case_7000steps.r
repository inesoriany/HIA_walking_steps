#################################################
##########      BEST CASE SCENARIO      #########
##########          Monte Carlo         #########
#################################################




###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                       HIA - 7000 steps                                                                                  #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  purrr,        # Loop
  stringr,      # Text extraction
  srvyr,        # Survey
  survey,
  ggplot2,      # Data visualization
  progress      # Progression bar
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Risk reduction distribution table
reduction_risk_distrib_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))

# Duration of depression distribution table
dep_duration_distrib_table <- import(here("data_clean", "Diseases", "dep_duration_distrib_table.xlsx"))

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
#                                                    3. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
BEST_emp_long <- emp_long %>% 
  # Recommendation of 7000 steps and baseline at 2000
  mutate(step = 7000 + baseline_step)



# EMP Dataset per disease
BEST_replicate_list <- list() 

for (dis in dis_vec) {
  BEST_replicate_list[[dis]] <- BEST_emp_long %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                      ALL DISEASES                          #
##############################################################
# Total of prevented burden of each disease for each simulation 
set.seed(123) 
BEST_burden_total <- HIA_burden_total(BEST_replicate_list, calc_HIA_replicate, incidence_distrib_table, dep_duration_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 prop_relapse, duration_recovery,vsl, 
                                 group = NULL,
                                 N = 1000)

  
# Export : Table of HIA outcomes per simulation
export(BEST_burden_total, here("output", "RDS", "7000 steps", "HIA_7000steps_1000replicate.rds"))


##############################################################
#                        PER SEX                             #
##############################################################
# Total of prevented burden of each disease per sex for each simulation
set.seed(123)
BEST_burden_sex_total <- HIA_burden_total(BEST_replicate_list, calc_HIA_replicate, incidence_distrib_table, dep_duration_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 prop_relapse, duration_recovery,vsl, 
                                 group ="sex", 
                                 N = 1000)


# Export : Table of HIA outcomes per simulation
export(BEST_burden_sex_total, here("output", "RDS", "7000 steps", "HIA_sex_7000_steps_1000replicate.rds"))



##############################################################
#                        PER AGE                             #
##############################################################
# Total of prevented burden of each disease per age for each simulation
set.seed(123)
BEST_burden_age_total <- HIA_burden_total(BEST_replicate_list, calc_HIA_replicate, incidence_distrib_table, dep_duration_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                     dis_vec, 
                                     prop_relapse, duration_recovery, vsl, 
                                     group = "age_grp10", 
                                     N = 1000)

# Export : Table of HIA outcomes per simulation
export(BEST_burden_age_total, here("output", "RDS", "7000 steps", "HIA_age_7000steps_1000replicate.rds"))




################################################################################################################################
#                         5. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
##############################################################

# Import data
BEST_burden_total <- import(here("output", "RDS", "7000 steps", "HIA_7000steps_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  BEST_burden_per_disease <- HIA_burden_IC(BEST_burden_total, dis_vec, outcome_vec, calc_replicate_IC) 


  # Total for morbidity
  BEST_burden_morbidity <- BEST_burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  
  # Total for all diseases
  BEST_burden_global <- BEST_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  BEST_burden <- bind_rows(BEST_burden_per_disease, BEST_burden_morbidity, BEST_burden_global)
  

  
  
# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
  # Per disease
  BEST_Rubin_burden_per_disease <- HIA_burden_IC(BEST_burden_total, dis_vec, outcome_vec, calc_IC_Rubin) 

  # Total for morbidity
  BEST_Rubin_burden_morbidity <- BEST_Rubin_burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  # Total for all diseases
  BEST_Rubin_burden_global <- BEST_Rubin_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
 BEST_Rubin_burden <- bind_rows(BEST_Rubin_burden_per_disease, BEST_Rubin_burden_morbidity, BEST_Rubin_burden_global)




# Import data
burden_2019 <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))
BEST_burden <- import(here("output", "Tables", "7000 steps", "HIA_7000steps_1000replicate.xlsx"))


# Incremental benefits
BEST_burden_add <- BEST_burden %>% 
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
    select(disease, tot_cases, tot_cases_low, tot_cases_up, tot_daly, tot_daly_low, tot_daly_up,
           tot_medic_costs, tot_medic_costs_low, tot_medic_costs_up, tot_soc_costs, tot_soc_costs_low, tot_soc_costs_up)





##############################################################
#                         PER SEX                            #
##############################################################

# Import data
BEST_burden_sex_total <- import(here("output", "RDS", "7000 steps", "HIA_sex_7000steps_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median (Monte Carlo)
set.seed(123)
BEST_burden_per_sex <- HIA_burden_IC(BEST_burden_sex_total, dis_vec, outcome_vec, calc_replicate_IC)


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
BEST_Rubin_burden_per_sex <- HIA_burden_IC(BEST_burden_sex_total, dis_vec, outcome_vec, calc_IC_Rubin)




##############################################################
#                         PER AGE                            #
##############################################################

# Import data
BEST_burden_age_total <- import(here("output", "RDS", "7000 steps", "HIA_age_7000steps_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median (Monte Carlo)
set.seed(123)
BEST_burden_per_age <- HIA_burden_IC(BEST_burden_age_total, dis_vec, outcome_vec, calc_replicate_IC)


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
BEST_Rubin_burden_per_age <- HIA_burden_IC(BEST_burden_age_total, dis_vec, outcome_vec, calc_IC_Rubin)




################################################################################################################################
#                                                     9. VISUALIZATION                                                         #
################################################################################################################################

# Import 2019 data
burden_sex_2019 <- import(here("output", "Tables", "2019", "HIA_per_sex.xlsx"))


# Plot : Cases prevented had the 2019 population walked 7000 steps according to sex, compared to 2019 levels
plot_BEST_cases_prev <- 
  ggplot() +
  geom_bar(data = burden_sex_2019 %>%  
      mutate(disease = factor(disease, levels = c("mort", "cvd", "cancer", "diab2", "dem", "dep"))),
    mapping = aes(x = disease, y = tot_cases, fill = sex, alpha = "2019 baseline"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(data = burden_sex_2019,
    mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "2019 baseline"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_fill_manual(values = colors_sex) +

  geom_bar(data = BEST_burden_per_sex, 
    mapping = aes(x = disease, y = tot_cases, fill = sex, alpha = "7000 steps recommendation"),
    width = 0.7,
    position = position_dodge2(0.7),
    stat = "identity") +
  geom_errorbar(
    data = BEST_burden_per_sex,
    mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "7000 steps recommendation"),
    position = position_dodge(0.7),
    width = 0.25) +
  scale_alpha_manual(
    name = "Scenario",
    values = c("2019 baseline" = 1, "7000 steps recommendation" = 0.4)) +
  scale_x_discrete(labels = names_disease) + 
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal()

plot_BEST_cases_prev 




################################################################################################################################
#                                                      11. EXPORT DATA                                                         #
################################################################################################################################
# Plot
  ggsave(here("output", "Plots", "7000 steps", "7000steps_cases_prev.png"), plot = plot_BEST_cases_prev)


# HIA of best case scenario
  export(BEST_burden, here("output", "Tables", "7000 steps", "HIA_7000steps_1000replicate.xlsx"))
  export(BEST_Rubin_burden, here("output", "Tables", "7000 steps", "HIA_Rubin_7000steps_1000replicate.xlsx"))
  export(BEST_burden_per_sex, here("output", "Tables", "7000 steps", "HIA_sex_7000steps_1000replicate.xlsx"))
  export(BEST_Rubin_burden_per_sex, here("output", "Tables", "7000 steps", "HIA_sex_Rubin_7000steps_1000replicate.xlsx"))
  export(BEST_burden_per_age, here("output", "Tables", "7000 steps", "HIA_age_7000steps_1000replicate.xlsx"))
  export(BEST_Rubin_burden_per_age, here("output", "Tables", "7000 steps", "HIA_age_Rubin_7000steps_1000replicate.xlsx"))
  export(BEST_burden_add, here("output", "Tables", "7000 steps", "HIA_added_7000steps_1000replicate.xlsx"))
