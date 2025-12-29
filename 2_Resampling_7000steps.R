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
dis_vec <- c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")


# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")

# Age categories
age_vec <- c("20-24", "25-34", "35-44", "45-54", "55-64", "65-75", "75-89")

# Sex categories
sex_vec <- c("Female", "Male")



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
#                      ALL DISEASES                          #
##############################################################
# Total of prevented burden of each disease for each simulation 
RECO_burden_replicate <- burden_replicate_prevented(data_list = RECO_HIA_replicate_list,
                                                    dis_vec,
                                                    group = NULL,
                                                    N = 1000)


# Export : Table of HIA outcomes per simulation
export(RECO_burden_replicate, here("output", "RDS", "7000 steps", "Resampling", "HIA_1000replicate_7000steps.rds"))



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
RECO_burden_per_disease <- HIA_burden_IC(RECO_burden_replicate, dis_vec, NULL, NULL, outcome_vec, calc_replicate_IC) %>% 
  select(-c(age_grp10, sex))

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
RECO_Rubin_burden_per_disease <- HIA_burden_IC(RECO_burden_replicate, dis_vec, NULL, NULL, outcome_vec, calc_IC_Rubin) %>% 
  select(-c(age_grp10, sex))

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



################################################################################################################################
#                                                      11. EXPORT DATA                                                         #
################################################################################################################################

# Tables of HIA outcomes per simulation
  export(RECO_burden_replicate, here("output", "RDS", "7000 steps", "Resampling", "HIA_1000replicate_7000steps.rds"))


# Tables of HIA outcomes
  export(RECO_burden, here("output", "Tables", "7000 steps", "Resampling", "HIA_per_disease_7000steps.xlsx"))
  export(RECO_Rubin_burden, here("output", "Tables", "7000 steps", "Resampling", "HIA_per_disease_Rubin_7000steps.xlsx"))






