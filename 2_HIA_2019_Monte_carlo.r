#################################################
##############      HIA 2019        #############
##############     Monte Carlo      #############
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



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                      ALL DISEASES                          #
##############################################################
# Total of prevented burden of each disease for each simulation 
set.seed(123)
burden_total <- HIA_burden_total(replicate_list, incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 vsl, 
                                 group = NULL,
                                 N = 1000)
 

# Export : Table of HIA outcomes per simulation
export(burden_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_1000replicate.rds"))


##############################################################
#                        PER SEX                             #
##############################################################
# Total of prevented burden of each disease per sex for each simulation
set.seed(123)
burden_sex_total <- HIA_burden_total(replicate_list, incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 vsl, 
                                 group ="sex", 
                                 N = 1000)


# Export : Table of HIA outcomes per simulation
export(burden_sex_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_per_sex_1000replicate.rds"))



##############################################################
#                        PER AGE                             #
##############################################################
# Total of prevented burden of each disease per age for each simulation
set.seed(123)
burden_age_total <- HIA_burden_total(replicate_list, incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                     dis_vec, 
                                     vsl, 
                                     group = "age_grp10", 
                                     N = 1000)

# Export : Table of HIA outcomes per simulation
export(burden_age_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_per_age_1000replicate.rds"))




################################################################################################################################
#                         5. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
##############################################################

# Import data
burden_total <- import(here("output", "RDS", "2019", "Monte Carlo", "HIA_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  burden_per_disease <- HIA_burden_IC(burden_total, dis_vec, outcome_vec, calc_replicate_IC) 


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
  Rubin_burden_per_disease <- HIA_burden_IC(burden_total, dis_vec, outcome_vec, calc_IC_Rubin) 

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
#                         PER SEX                            #
##############################################################

# Import data
burden_sex_total <- import(here("output", "RDS", "2019", "Monte Carlo", "HIA_per_sex_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median (Monte Carlo)
set.seed(123)
burden_per_sex <- HIA_burden_IC(burden_sex_total, dis_vec, outcome_vec, calc_replicate_IC)


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
Rubin_burden_per_sex <- HIA_burden_IC(burden_sex_total, dis_vec, outcome_vec, calc_IC_Rubin)



##############################################################
#                         PER AGE                            #
##############################################################

# Import data
burden_age_total <- import(here("output", "RDS", "2019", "Monte Carlo", "HIA_per_age_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median (Monte Carlo)
set.seed(123)
burden_per_age <- HIA_burden_IC(burden_age_total, dis_vec, outcome_vec, calc_replicate_IC)


# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
Rubin_burden_per_age <- HIA_burden_IC(burden_age_total, dis_vec, outcome_vec, calc_IC_Rubin)



################################################################################################################################
#                                              8. REDUCTION IN MORTALITY RISK                                                  #
################################################################################################################################

##############################################################
#               ‰ REDUCTION IN MORTALITY RISK                #
##############################################################
set.seed(123)

N <- 10
reduc_mortality_risk <- data.frame()

for (i in 1:N) {

  HIA_mortality_replicate <- calc_HIA_replicate(replicate_list,
                                                incidence_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                                                "mort",
                                                vsl)

  reduc_mortality_risk_replicate <- HIA_mortality_replicate[["mort"]] %>% 
    as_survey_design(ids = ident_ind, weights = pond_indc) %>% 
    summarise(mean_mort_reduction_risk = survey_mean(reduction_risk, na.rm = TRUE)) %>% 
    mutate(simulation_id = i)

  reduc_mortality_risk <- bind_rows(reduc_mortality_risk, reduc_mortality_risk_replicate)

}

# Export reduction in mortality risk for 1000 replications
export(reduc_mortality_risk, here("output", "RDS", "2019", "Monte Carlo",  "reduc_mortality_risk_1000_rep.RDS"))



# Load reduction in mortality risk for 1000 replications
reduc_mortality_risk <- import(here("output", "RDS", "2019", "Monte Carlo",  "reduc_mortality_risk_1000_rep.RDS"))

# IC95 and median mean
N <- 10
IC <-  calc_replicate_IC(reduc_mortality_risk, "mean_mort_reduction_risk")
reduc_mortality_risk_IC <- data.frame(
  reduc_mortality_risk = paste0(round(IC["50%"], 3), " (", round(IC["2.5%"], 3), " - ", round(IC["97.5%"],3),  ")"),
  N_replications = N)




################################################################################################################################
#                                                     9. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Cases prevented by walking in 2019 according to sex
plot_cases_prev <- burden_per_sex %>% 
  mutate(disease = factor(disease, levels = c("mort", "cvd", "cancer", "diab2", "dem", "dep"))) %>%
  ggplot(aes(x = disease, y = tot_cases, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x.top = element_blank(),      # delete labels X at the top
    axis.ticks.x.top = element_blank()      # delete ticks X at the top
  )

plot_cases_prev



# Plot : Median DALY prevented by walking in 2019 according to age group
plot_daly_prevented <- burden_per_age %>% 
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
km_total_2019 <- as.numeric(svytotal(~nbkm_tot_walking, jour)) *365.25/7                              # Total km per year
km_total_2019_IC <- confint(svytotal(~nbkm_tot_walking, jour) *365.25/7 )                             # Confidence interval

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
  ggsave(here("output", "Plots", "2019", "Monte Carlo", "cases_prevented.png"), plot = plot_cases_prev)
  ggsave(here("output", "Plots", "2019", "Monte Carlo", "DALY_prevented.png"), plot = plot_daly_prevented)

  
# Tables of HIA outcomes per simulation
  export(burden_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_1000replicate.rds"))
  export(burden_sex_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_per_sex_1000replicate.rds"))
  export(burden_age_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_per_age_1000replicate.rds"))

  
# Tables of HIA outcomes
  export(burden, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_disease.xlsx"))
  export(Rubin_burden, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_disease_Rubin.xlsx"))
  export(burden_per_sex, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_sex.xlsx"))
  export(burden_per_sex, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_sex_Rubin.xlsx"))
  export(burden_per_age, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_age.xlsx"))
  export(Rubin_burden_per_age, here("output", "Tables", "2019", "Monte Carlo", "HIA_per_age_Rubin.xlsx"))

  
# Table of reduction of mortality risk
  export(reduc_mortality_risk_IC, here("output", "Tables", "2019", "Monte Carlo", "reduc_mortality_risk.xlsx"))
  
  
# Tables economic unit value
  # Economic value of 1 km walked
  export(unit_value_2019, here("output", "Tables", "2019", "Monte Carlo", "1km_value_1000replicate.xlsx"))
  # Social economic value of 1 km walked
  export(unit_soc_value_2019, here("output", "Tables", "2019", "Monte Carlo", "1km_soc_value_1000replicate.xlsx"))
  
  # Number of steps, distance and duration to save 1€ of medical costs in 2019
  export(euro_unit_2019, here("output", "Tables", "2019", "Monte Carlo", "1€_step_km_min_1000replicate.xlsx"))
  # Number of steps, distance and duration to save 1€ of intangible costs in 2019
  export(soc_euro_unit_2019, here("output", "Tables", "2019", "Monte Carlo", "soc_1€_step_km_min_1000replicate.xlsx"))





