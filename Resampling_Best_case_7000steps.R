#################################################
################## RESAMPLING ###################
#################################################



###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                      HIA - 7000 STEPS                                                                                   #
###########################################################################################################################################################################
###########################################################################################################################################################################


################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
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
emp_step <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_distrib_table <- import(here("data_clean", "Diseases", "DRF", "rr_sim_interpolated.rds"))


# Disability weights
dw_distrib_table <- import(here("data_clean", "Diseases", "dw_sim.rds"))

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
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_step <- emp_step %>% 
  # Recommendation of 7000 steps and baseline at 2000
  mutate(step = 7000 + baseline_step)



# EMP Dataset per disease
RECO_replicate_list <- list() 

for (dis in dis_vec) {
  RECO_replicate_list[[dis]] <- emp_step %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################
# Calculate for each individual the number of prevented cases, DALY and costs with the 1000 simulated parameters
RECO_HIA_replicate_list <- calc_HIA_replicate(data_list = RECO_replicate_list,
                                              rr_distrib_table = rr_distrib_table,
                                              dw_distrib_table = dw_distrib_table,
                                              dis_vec = dis_vec,
                                              vsl,
                                              baseline_step)



################################################################################################################################
#                                               5. TOTAL BURDEN PER SIMULATION                                                 #
################################################################################################################################

##############################################################
#                       ALL DISEASES                         #
##############################################################
# Total of prevented burden of each disease for each simulation 
RECO_burden_replicate <- burden_replicate_prevented(data_list = RECO_HIA_replicate_list,
                                                    dis_vec,
                                                    group = NULL,
                                                    N = 1000)


# Export : Table of HIA outcomes per simulation
export(RECO_burden_replicate, here("output", "RDS", "7000 steps", "Resampling", "HIA_1000replicate_7000steps.rds"))



##############################################################
#                         PER SEX                            #
##############################################################
# Total of prevented burden of each disease for each simulation 
RECO_burden_replicate_sex <- burden_replicate_prevented(data_list = RECO_HIA_replicate_list,
                                                    dis_vec,
                                                    group = "sex",
                                                    N = 1000)


# Export : Table of HIA outcomes per simulation
export(RECO_burden_replicate_sex , here("output", "RDS", "7000 steps", "Resampling", "HIA_per_sex_1000replicate_7000steps.rds"))




################################################################################################################################
#                         6. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       ALL DISEASES                         #
##############################################################

# Import data
RECO_burden_replicate <- import(here("output", "RDS", "7000 steps", "Resampling", "HIA_1000replicate_7000steps.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
# Per disease
set.seed(123)
RECO_burden_per_disease <- HIA_burden_IC(RECO_burden_replicate, dis_vec, outcome_vec, calc_replicate_IC) 

# Total for morbidity
RECO_burden_morbidity <- RECO_burden_per_disease %>%
  filter(disease != "mort") %>% 
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "Morbidity") %>%
  select(disease, everything()) 


# Total for all diseases
RECO_burden_global <- RECO_burden_per_disease %>%
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "All") %>%
  select(disease, everything()) 

# Gather results
RECO_burden <- bind_rows(RECO_burden_per_disease, RECO_burden_morbidity, RECO_burden_global)




# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
# Per disease
RECO_Rubin_burden_per_disease <- HIA_burden_IC(RECO_burden_replicate, dis_vec, outcome_vec, calc_IC_Rubin)

# Total for morbidity
RECO_Rubin_burden_morbidity <- RECO_Rubin_burden_per_disease %>%
  filter(disease != "mort") %>% 
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "Morbidity") %>%
  select(disease, everything()) 

# Total for all diseases
RECO_Rubin_burden_global <- RECO_Rubin_burden_per_disease %>%
  summarise(across(where(is.numeric), 
                   ~ sum(.x, na.rm = TRUE) )) %>%
  mutate(disease = "All") %>%
  select(disease, everything()) 

# Gather results
RECO_Rubin_burden <- bind_rows(RECO_Rubin_burden_per_disease,RECO_Rubin_burden_morbidity, RECO_Rubin_burden_global)



##############################################################
#                  ALL DISEASES (Additional)                 #
##############################################################
# Import 2019 data
burden_2019 <- import(here("output", "Tables", "2019", "Resampling", "HIA_per_disease.xlsx"))


# Data preparation
burden_2019_row <- burden_2019 %>% 
  rename_with(.fn = ~ paste0(.x, "_2019"),
              .cols = -c (disease))

RECO_burden_row <- RECO_burden %>% 
  rename_with(.fn = ~ paste0(.x, "_RECO"),
              .cols = -c (disease))


# Additional prevented cases for each disease
add_RECO_burden <- burden_2019_row %>%
  left_join(RECO_burden_row, by = "disease", suffix = c("_2019", "_RECO")) %>%
  mutate(across(
    matches("_RECO$"),
    ~ . - get(sub("_RECO$", "_2019", cur_column())),
    .names = "{.col}_add"
  ))



##############################################################
#                         PER SEX                            #
##############################################################
# Import data
RECO_burden_replicate_sex <- import(here("output", "RDS", "7000 steps", "Resampling", "HIA_per_sex_1000replicate_7000steps.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
# Per disease
set.seed(123)
RECO_burden_per_sex <- HIA_burden_IC(RECO_burden_replicate_sex , dis_vec, outcome_vec, calc_replicate_IC) 


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
# Per sex
RECO_Rubin_burden_per_sex <- HIA_burden_IC(RECO_burden_replicate_sex, dis_vec, outcome_vec, calc_IC_Rubin)




################################################################################################################################
#                                                     11. VISUALIZATION                                                        #
################################################################################################################################
# Import 2019 data
burden_sex_2019 <- import(here("output", "Tables", "2019", "Resampling", "HIA_per_sex.xlsx"))


# Plot : Cases prevented had the 2019 population walked 7000 steps according to sex, compared to 2019 levels
plot_RECO_cases_prev <- 
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
  
  
  geom_bar(data = RECO_burden_per_sex, 
           mapping = aes(x = disease, y = tot_cases, fill = sex, alpha = "7000 steps recommendation"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "7000 steps recommendation" = 0.4)) +
  
  geom_errorbar(data = RECO_burden_per_sex,
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "7000 steps recommendation"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal() 

plot_RECO_cases_prev 




################################################################################################################################
#                                                      12. DESCRIPTION                                                         #
################################################################################################################################
# DALYs comparison between 2019 and best case scenario
DALY_prev_fraction <- RECO_burden_row %>% 
  inner_join(burden_2019_row,
    by = c("disease")) %>% 
  mutate(
    daly_mid_fraction = tot_daly_RECO / tot_daly_2019,
    daly_low_fraction = tot_daly_RECO / tot_daly_low_2019,
    daly_up_fraction  = tot_daly_RECO  / tot_daly_up_2019)  %>% 
  select(disease, starts_with("daly"))


################################################################################################################################
#                                                      13. EXPORT DATA                                                         #
################################################################################################################################

# Plot
  ggsave(here("output", "Plots", "7000 steps", "cases_prev_7000steps.png"), plot = plot_RECO_cases_prev)

# Tables of HIA outcomes per simulation
  export(RECO_burden_replicate, here("output", "RDS", "7000 steps", "Resampling", "HIA_1000replicate_7000steps.rds"))


# Tables of HIA outcomes
  export(RECO_burden, here("output", "Tables", "7000 steps", "Resampling", "HIA_per_disease_7000steps.xlsx"))
  export(RECO_Rubin_burden, here("output", "Tables", "7000 steps", "Resampling", "HIA_per_disease_Rubin_7000steps.xlsx"))
