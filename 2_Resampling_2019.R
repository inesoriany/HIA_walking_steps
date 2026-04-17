#################################################
################## RESAMPLING ###################
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
emp_walk <- import(here("data_clean", "EMP_walkers.xlsx"))
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
  mutate(step = pmin(12000, round(step_commute / 100) * 100 + baseline_step))


# EMP Dataset per disease
replicate_list <- list() 

for (dis in dis_vec) {
  replicate_list[[dis]] <- emp_long %>% 
    filter(disease == dis)
}


##############################################################
#               DISEASE INCIDENCE DISTRIBUTIONS              #
##############################################################

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
                                               group = NULL,
                                               N = 1000)
  
  
# Export : Table of HIA outcomes per simulation
export(burden_replicate, here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))


##############################################################
#                        PER AGE                             #
##############################################################
# Total of prevented burden per age for each simulation
burden_replicate_age <- burden_replicate_prevented(data_list = HIA_replicate_list,
                                                   dis_vec,
                                                   group = "age_grp10",
                                                   N = 1000)


# Export : Table of HIA outcomes per simulation
export(burden_replicate_age, here("output", "RDS", "2019", "Resampling", "HIA_per_age_1000replicate.rds"))





################################################################################################################################
#                         6. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
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



################################################################################################################################
#                                              8. REDUCTION IN MORTALITY RISK                                                  #
################################################################################################################################
##############################################################
#               ‰ REDUCTION IN MORTALITY RISK                #
##############################################################
# Mortality reduction risk per RR simulated
reduc_mortality_risk <- HIA_replicate_list[["mort"]] %>% 
  as_survey_design(ids = ident_ind, weights = pond_indc) %>% 
  group_by(simulation_id) %>% 
  summarise(
    mean_mort_reduction_risk = survey_mean(reduction_risk, na.rm = TRUE)
  ) %>% 
  ungroup()


# Export reduction in mortality risk for 1000 replications
export(reduc_mortality_risk, here("output", "RDS", "Log linear", "2019",  "reduc_mortality_risk_1000_rep.RDS"))



# Load reduction in mortality risk for 1000 replications
reduc_mortality_risk <- import(here("output", "RDS", "Log linear", "2019",  "reduc_mortality_risk_1000_rep.RDS"))

# IC95 and mean
N = 1000
IC <-  calc_replicate_IC(reduc_mortality_risk, "mean_mort_reduction_risk")
reduc_mortality_risk_IC <- data.frame(
  reduc_mortality_risk = paste0(round(IC["50%"], 3), " (", round(IC["2.5%"], 3), " - ", round(IC["97.5%"],3),  ")"),
  N_replications = N)



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
# Survey design ponderated by day
jour <- emp_walk %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)

# Total walked distance in 2019
km_total_2019 <- as.numeric(svytotal(~nbkm_walking, jour)) *365.25/7                              # Total km per year
km_total_2019_IC <- confint(svytotal(~nbkm_walking, jour) *365.25/7 )                             # Confidence interval

# Total steps in 2019
step_total_2019 <- as.numeric(svytotal(~step_commute, jour)) *365.25/7                            # Total steps per year
step_total_2019_IC <- confint(svytotal(~step_commute, jour) *365.25/7 )                           # Confidence interval


# Setting parameters 
km_low_2019 <- km_total_2019_IC[1, 1]
km_sup_2019 <- km_total_2019_IC[1, 2]

step_low_2019 <- step_total_2019_IC[1, 1]
step_sup_2019 <- step_total_2019_IC[1, 2]

euro <- burden_global[["tot_medic_costs"]] [[1]] 
euro_low <- burden_global[["tot_medic_costs_low"]] [[1]]
euro_sup <- burden_global[["tot_medic_costs_sup"]] [[1]] 

soc_euro <- burden_global[["tot_soc_costs"]] [[1]] 
soc_euro_low <- burden_global[["tot_soc_costs_low"]] [[1]] 
soc_euro_sup <- burden_global[["tot_soc_costs_sup"]] [[1]] 



##############################################################
#                        VALUE OF 1km                        #
##############################################################
# --------------------------------------
# MEDICAL COSTS
# --------------------------------------
  # Calculate economic value of 1 km walked (medical costs) (in €)
  set.seed(123)
  unit_2019 <- unit_value(km_total_2019, km_low_2019, km_sup_2019, euro, euro_low, euro_sup, N=1000)
  unit_value_2019 <- as.data.frame(t(quantile(unit_2019, probs = c(0.025, 0.5, 0.975)))) %>% 
    rename(euro_2.5 = "2.5%",
           euro_50 = "50%",
           euro_97.5 = "97.5%") %>% 
    mutate(km = 1)


# --------------------------------------
# SOCIAL COSTS
# --------------------------------------
  # Calculate economic value of 1 km walked (intangible costs) (in €)
  set.seed(123)
  unit_soc_2019 <- unit_value(km_total_2019, km_low_2019, km_sup_2019, soc_euro, soc_euro_low, soc_euro_sup, N=1000)
  unit_soc_value_2019 <- as.data.frame(t(quantile(unit_soc_2019, probs = c(0.025, 0.5, 0.975)))) %>% 
    rename(soc_euro_2.5 = "2.5%",
           soc_euro_50 = "50%",
           soc_euro_97.5 = "97.5%") %>% 
    mutate(km = 1)



##############################################################
#                           SAVE 1€                          #
##############################################################
# --------------------------------------
# MEDICAL COSTS
# --------------------------------------

# Calculate distance walked to save 1€ of medical costs (km)
  set.seed(123)
  euro_km_2019 <- euro_km_unit(km_total_2019, km_low_2019, km_sup_2019, euro, euro_low, euro_sup, N = 1000)
  euro_km_unit_2019 <- as.data.frame(t(quantile(euro_km_2019, probs = c(0.025, 0.5, 0.975))))

# Calculate number of steps to save 1€ of medical costs (km)
  set.seed(123)
  euro_step_2019 <- euro_step_unit(step_total_2019, step_low_2019, step_sup_2019, euro, euro_low, euro_sup, N = 1000)
  euro_step_unit_2019 <- as.data.frame(t(quantile(euro_step_2019, probs = c(0.025, 0.5, 0.975))))
  
# Save 1€ of medical costs
euro_unit_2019 <- euro_step_unit_2019 %>% 
  rename(step_2.5  = "2.5%",
           step_50   = "50%",
           step_97.5 = "97.5%") %>% 
  bind_cols(euro_km_unit_2019 %>% 
              rename(km_2.5  = "2.5%",
                     km_50   = "50%",
                     km_97.5 = "97.5%")) %>% 
  mutate(min_2.5  = km_2.5  * 60 / walk_speed,              # Duration walked to save 1€ of medical costs (min)
         min_50   = km_50   * 60 / walk_speed,
         min_97.5 = km_97.5 * 60 / walk_speed,
         medic_costs = 1)
    


# --------------------------------------
# SOCIAL COSTS
# --------------------------------------
# Calculate distance walked to save 1€ of intangible costs (km)
  set.seed(123)
  soc_euro_km_2019 <- euro_km_unit(km_total_2019, km_low_2019, km_sup_2019, soc_euro, soc_euro_low, soc_euro_sup, N = 1000)
  soc_euro_km_unit_2019 <- as.data.frame(t(quantile(soc_euro_km_2019, probs = c(0.025, 0.5, 0.975))))
  
# Calculate number of steps walked to save 1€ of intangible costs (km)
  set.seed(123)
  soc_euro_step_2019 <- euro_step_unit(step_total_2019, step_low_2019, step_sup_2019, soc_euro, soc_euro_low, soc_euro_sup, N = 1000)
  soc_euro_step_unit_2019 <- as.data.frame(t(quantile(soc_euro_step_2019, probs = c(0.025, 0.5, 0.975))))

# Save 1€ of medical costs
soc_euro_unit_2019 <- soc_euro_step_unit_2019 %>% 
  rename(step_2.5  = "2.5%",
           step_50   = "50%",
           step_97.5 = "97.5%") %>% 
  bind_cols(soc_euro_km_unit_2019 %>% 
              rename(km_2.5  = "2.5%",
                     km_50   = "50%",
                     km_97.5 = "97.5%")) %>% 
  mutate(min_2.5  = km_2.5  * 60 / walk_speed,                # Duration walked to save 1€ of intangible costs (min)
         min_50   = km_50   * 60 / walk_speed,
         min_97.5 = km_97.5 * 60 / walk_speed,
         soc_costs = 1)
 



################################################################################################################################
#                                                      11. EXPORT DATA                                                         #
################################################################################################################################
# Plot
  ggsave(here("output", "Plots", "2019", "DALY_prevented.png"), plot = plot_daly_prevented)


  
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

  
# Table of reduction of mortality risk
  export(reduc_mortality_risk_IC, here("output", "Tables", "2019", "Resampling", "reduc_mortality_risk.xlsx"))
  
  
# Tables economic unit value
  # Economic value of 1 km walked
  export(unit_value_2019, here("output", "Tables", "2019", "Resampling", "1km_value_1000replicate.xlsx"))
  # Social economic value of 1 km walked
  export(unit_soc_value_2019, here("output", "Tables", "2019", "Resampling", "1km_soc_value_1000replicate.xlsx"))
  
  # Number of steps, distance and duration to save 1€ of medical costs in 2019
  export(euro_unit_2019, here("output", "Tables", "2019", "Resampling", "1€_step_km_min_1000replicate.xlsx"))
  # Number of steps, distance and duration to save 1€ of intangible costs in 2019
  export(soc_euro_unit_2019, here("output", "Tables", "2019", "Resampling", "soc_1€_step_km_min_1000replicate.xlsx"))
  


