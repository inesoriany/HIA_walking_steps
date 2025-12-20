#################################################
############ HEALTH IMPACT ASSESSMENT ###########
#################################################



###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                          HIA - 2019                                                                                     #
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
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_distrib_table <- import(here("data_clean", "Diseases", "DRF", "rr_sim_interpolated.rds"))


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
dis_vec = c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")


# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")

# Age categories
age_vec <- c("20-24", "25-34", "35-44", "45-54", "55-64", "65-75", "75-89")

# Sex categories
sex_vec <- c("Female", "Male")




################################################################################################################################
#                                                    3. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_long <- emp_long %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000 steps
  mutate(step = pmin(12000, round(step_commute / 100) * 100 + 2000))


# EMP Dataset per disease
replicate_list <- list() 

for (dis in dis_vec) {
  replicate_list[[dis]] <- emp_long %>% 
    filter(disease == dis)
}



##############################################################
#              DISABILITY WEIGHTS DISTRIBUTIONS              #
##############################################################

# Generate DW normal distributions
set.seed(123)
dw_distrib_table <- dw_table %>%
  rowwise() %>%
  mutate(dw_distrib = list(
    generate_RR_distrib(dw_mid, dw_low, dw_up, 1000)
  ))



# Sort the DW distributions in ascending order
dw_distrib_table <- dw_distrib_table %>% 
  rowwise() %>% 
  mutate(dw_distrib = list(sort(unlist(dw_distrib)))) %>%
  ungroup()


# Separate the dw_distrib column into multiple columns
dw_distrib_table <- dw_distrib_table  %>% 
  unnest_wider(dw_distrib, names_sep = "_") %>% 
  pivot_longer(
    cols = starts_with("dw_distrib_"), 
    names_to = "simulation_id",                     # column name for the simulation ID
    values_to = "simulated_dw"                      # column name for the simulated RR values
  ) %>% 
  mutate(simulation_id = as.numeric(str_remove(simulation_id, "dw_distrib_")))      # simulation ID as a numeric value




################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################
# Calculate for each individual the number of prevented cases, DALY and costs with the 1000 simulated parameters
HIA_replicate_list <- calc_HIA_replicate(data_list = replicate_list,
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
burden_replicate <- burden_replicate_prevented(data_list = HIA_replicate_list,
                                               dis_vec,
                                               group = NULL)
  
  
# Export : Table of HIA outcomes per simulation
export(burden_replicate, here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))


##############################################################
#                        PER AGE                             #
##############################################################
# Total of prevented burden per age for each simulation
burden_replicate_age <- burden_replicate_prevented(data_list = HIA_replicate_list,
                                                   dis_vec,
                                                   group = "age_grp10")


# Export : Table of HIA outcomes per simulation
export(burden_replicate_age, here("output", "RDS", "2019", "Resampling", "HIA_per_age_1000replicate.rds"))





################################################################################################################################
#                         6. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       ALL DISEASES                         #
##############################################################

# Import data
burden_replicate <- import(here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  burden_per_disease <- HIA_burden_IC(burden_replicate, dis_vec, NULL, NULL, outcome_vec, calc_replicate_IC) %>% 
    select(-c(age_grp10, sex))

  # Total for morbidity
  burden_morbidity <- burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  
  # Total for all diseases
  burden_global <- burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  burden <- bind_rows(burden_per_disease, burden_morbidity, burden_global)
  

  
  
# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
  # Per disease
  Rubin_burden_per_disease <- HIA_burden_IC(burden_replicate, dis_vec, NULL, NULL, outcome_vec, calc_IC_Rubin) %>% 
    select(-c(age_grp10, sex))

  # Total for morbidity
  Rubin_burden_morbidity <- Rubin_burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  # Total for all diseases
  Rubin_burden_global <- Rubin_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  Rubin_burden <- bind_rows(Rubin_burden_per_disease,Rubin_burden_morbidity, Rubin_burden_global)



##############################################################
#                         PER AGE                            #
##############################################################

# Import data
burden_replicate_age <- import(here("output", "RDS", "2019", "Resampling", "HIA_per_age_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median (Monte Carlo)
set.seed(123)
burden_per_age <- HIA_burden_IC(burden_replicate_age, dis_vec, age_vec, NULL, outcome_vec, calc_replicate_IC) %>% 
  select(-c(sex))
  


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
Rubin_burden_per_age <- HIA_burden_IC(burden_replicate_age, dis_vec, age_vec, NULL, outcome_vec, calc_IC_Rubin) %>% 
  select(-c(sex))



##############################################################
#                         PER SEX ?                          #
##############################################################




################################################################################################################################
#                                              8. REDUCTION IN MORTALITY RISK                                                  #
################################################################################################################################




################################################################################################################################
#                                                     9. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Median DALY prevented by walking in 2019 according to age group
plot_daly_prevented <- burden_per_age %>% filter(disease != "bc") %>% 
  ggplot(aes(x = age_grp10, y = tot_daly, fill = disease)) +
  geom_bar(width = 0.7, position = "stack", stat = "identity")  +
  scale_fill_manual(values = colors_disease, labels = names_disease) +
  xlab("Age group") +
  ylab("Median DALY") +
  theme_minimal()

plot_daly_prevented



################################################################################################################################
#                                               10. ECONOMIC UNIT VALUE (€)                                                    #
################################################################################################################################



################################################################################################################################
#                                                      11. EXPORT DATA                                                         #
################################################################################################################################
# Tables of disability weights distribution
export(dw_distrib_table, here("data_clean", "Diseases", "dw_sim.rds"))

# Tables of HIA outcomes per simulation
export(burden_replicate, here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))
export(burden_replicate_age, here("output", "RDS", "2019", "Resampling", "HIA_per_age_1000replicate.rds"))


# Tables of HIA outcomes
export(burden, here("output", "Tables", "2019", "Resampling", "HIA_per_disease.xlsx"))
export(Rubin_burden, here("output", "Tables", "2019", "Resampling", "HIA_per_disease_Rubin.xlsx"))
export(burden_per_age, here("output", "Tables", "2019", "Resampling", "HIA_per_age.xlsx"))
export(Rubin_burden_per_age, here("output", "Tables", "2019", "Resampling", "HIA_per_age_Rubin.xlsx"))

       
# Plot
ggsave(here("output", "Plots", "2019", "DALY_prevented.png"), plot = plot_daly_prevented)






