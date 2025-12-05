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
  srvyr,        # Survey
  survey,
  ggplot2,      # Data visualization
  cli           # Progression bar
)




################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################
# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_distrib_table <- import(here("data_clean", "DRF", "rr_sim_interpolated.rds"))


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
#                                                 4. SIMULATED REDUCTION RR                                                    #
################################################################################################################################

## ALL DISEASES (EXCEPT BREAST CANCER)
  # Baseline step 2000
  rr_baseline <- rr_distrib_table %>%
    filter(step == 2000) %>%
    select(disease, simulation_id, rr2000 = rr_interpolated)
  
  
  # Reduction risk 
  rr_distrib_table <- rr_distrib_table %>%
    left_join(rr_baseline, by = c("disease", "simulation_id")) %>%
    mutate(reduction_risk = (rr2000 - rr_interpolated) / rr2000)
  



################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_long <- emp_long %>% 
  # Round the number of steps to the nearest hundred
  mutate(step = pmin(12000, round(step_commute / 100) * 100 + 2000))


# EMP Dataset per disease
for (dis in dis_vec) {
  assign(
    paste0(dis, "_replicate"),
    emp_long %>% filter(disease == dis)
  )
}


################################################################################################################################
#                                                 5. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################
  
for (dis in dis_vec) {

  replicate_name <- paste0(dis, "_replicate")
  dis_replicate <- get(replicate_name)

  # --- breast cancer ---
  if (dis == "bc") {

    params <- dis_setting(dis)
    set.seed(123)

    rr_women_sim <- generate_RR_distrib(params$rr_women, params$rr_women_low, params$rr_women_up, N = 1000)
    rr_men_sim   <- generate_RR_distrib(params$rr_men, params$rr_men_low, params$rr_men_up, N = 1000)

    # Dupliquer chaque individu 1000 fois + assigner les 1000 RR
    dis_replicate <- dis_replicate %>%
      slice(rep(1:n(), each = 1000)) %>%
      group_by(ID) %>%
      mutate(
        simulation_id = 1:1000,
        rr_sim = if (sex[1] == "Female") rr_women_sim else rr_men_sim
      ) %>%
      ungroup()

    # Calcul réduction de risque
    dis_replicate <- log_reduction_risk(
      data = dis_replicate,
      dis = dis,
      rr_women = dis_replicate$rr_sim,
      rr_men   = dis_replicate$rr_sim,
      ref_women = params$ref_women,
      ref_men   = params$ref_men,
      week_base = week_base
    )

  } else {

    # --- other diseases ---
    dis_rr <- rr_distrib_table %>%
      filter(disease == dis) %>%
      select(simulation_id, step, reduction_risk)

    dis_replicate <- dis_replicate %>%
      left_join(dis_rr, by = "step", relationship = "many-to-many")
  }

  assign(replicate_name, dis_replicate)
}

  
  
  
  

##############################################################
#                      CASES PREVENTED                       #
##############################################################

# Calculate prevented cases
for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
    
  dis_replicate <- get(replicate_name) %>% 
    mutate(reduc_incidence = reduction_risk * rate)
     
  assign(replicate_name, dis_replicate)
}
  
  
    
  


##############################################################
#                    DISABILITY WEIGHTS                      #
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



# Randomly associate DW to each individual for each disease
set.seed(123)
for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
  
  dis_dw <- dw_distrib_table %>%
    filter(disease == dis) %>%
    pull(simulated_dw)
  
  dis_replicate <- get(replicate_name) %>%
    mutate(
      dw = sample(dis_dw, size = n(), replace = TRUE) 
    )
  
  assign(replicate_name, dis_replicate)
}





##############################################################     
#                           DALY                             #
##############################################################

for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
  dis_replicate <- get(replicate_name) 

    dis_replicate <- daly(dis_replicate)
    
  assign(replicate_name, dis_replicate)
}



  

##############################################################
#                    ECONOMIC IMPACT                         #
##############################################################

for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
  dis_replicate <- get(replicate_name) 

  dis_replicate <- medic_costs (dis_replicate, dis)
  dis_replicate <- dis_replicate %>% 
    mutate(soc_costs = daly*vsl)
  

  assign(replicate_name, dis_replicate)
}



################################################################################################################################
#                                               6. TOTAL BURDEN PER SIMULATION                                                 #
################################################################################################################################

##############################################################
#                      ALL DISEASES                          #
##############################################################

cli_progress_bar("Burden calculations", total = length(dis_vec) * 1000)

for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
  burden_replicate_name <- paste0(dis, "_burden_replicate")
  
  dis_burden_replicate <- list() 
  
  for(i in 1:1000) {
    cli_progress_update()  
    message("Simulation ", i, " ~ ", dis)
    
    # Survey design
    surv_dis_replicate <- get(replicate_name) %>% 
      filter(simulation_id == i) %>% 
      as_survey_design(ids = ident_ind, weights = pond_indc)
    
    # For 1 simulation, total burden
    dis_burden_replicate[[i]] <- surv_dis_replicate %>% 
      summarise(
        tot_cases = survey_total(reduc_incidence, na.rm = TRUE),
        tot_daly = survey_total(daly, na.rm = TRUE),
        tot_medic_costs = survey_total(medic_costs, na.rm = TRUE),
        tot_soc_costs = survey_total(soc_costs, na.rm = TRUE)
      ) %>%
      mutate(disease = dis, simulation_id = i)
  }
  
  # All simulations in a data.frame
  dis_burden_replicate <- bind_rows(dis_burden_replicate) %>% as.data.frame()
  
  assign(burden_replicate_name, dis_burden_replicate)
}



# Gather the results in a data.frame
burden_replicate <- data.frame()

for (dis in dis_vec) {
  burden_replicate_name <- paste0(dis, "_burden_replicate")
  
  dis_burden_replicate <- get(burden_replicate_name)
  
  burden_replicate <- bind_rows(burden_replicate, dis_burden_replicate)
}


# Export results
export(burden_replicate, here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))



##############################################################
#                         PER AGE                            #
##############################################################
cli_progress_bar("Burden calculations", total = length(dis_vec) * 1000)

for (dis in dis_vec) {
  replicate_name <- paste0(dis, "_replicate")
  burden_replicate_age_name <- paste0(dis, "_burden_replicate_age")
  
  dis_burden_replicate <- list() 
  
  for(i in 1:1000) {
    cli_progress_update()
    message("Simulation ", i, " ~ ", dis)
    
    surv_dis_replicate <- get(replicate_name) %>% 
      filter(simulation_id == i) %>% 
      as_survey_design(ids = ident_ind, weights = pond_indc)
    
    dis_burden_replicate[[i]] <- surv_dis_replicate %>% 
      group_by(age_grp10) %>% 
      summarise(
        tot_cases = survey_total(reduc_incidence, na.rm = TRUE),
        tot_daly = survey_total(daly, na.rm = TRUE),
        tot_medic_costs = survey_total(medic_costs, na.rm = TRUE),
        tot_soc_costs = survey_total(soc_costs, na.rm = TRUE)
      ) %>%
      mutate(disease = dis, simulation_id = i)
  }
  
  # Combiner toutes les simulations
  dis_burden_replicate <- bind_rows(dis_burden_replicate) %>% as.data.frame()
  assign(burden_replicate_age_name, dis_burden_replicate)
  
}


# Gather the results in a data.frame
burden_replicate_age <- data.frame()

for (dis in dis_vec) {
  burden_replicate_age_name <- paste0(dis, "_burden_replicate_age")
  
  dis_burden_replicate_age <- get(burden_replicate_age_name)
  
  burden_replicate_age <- bind_rows(burden_replicate_age, dis_burden_replicate_age)
}


# Export results
export(burden_replicate_age, here("output", "RDS", "2019", "Resampling", "HIA_per_age_1000replicate.rds"))





################################################################################################################################
#                         7. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       ALL DISEASES                         #
##############################################################

# Import data
burden_replicate <- import(here("output", "RDS", "2019", "Resampling", "HIA_1000replicate.rds"))


# IC95 and median (Monte Carlo)
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
  

  
  
# Rubin's rule
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


# IC95 and median (Monte Carlo)
set.seed(123)
burden_per_age <- HIA_burden_IC(burden_replicate_age, dis_vec, age_vec, NULL, outcome_vec, calc_replicate_IC) %>% 
  select(-c(sex))
  


# Rubin's rule
Rubin_burden_per_age <- HIA_burden_IC(burden_replicate_age, dis_vec, age_vec, NULL, outcome_vec, calc_IC_Rubin) %>% 
  select(-c(sex))



##############################################################
#                         PER SEX                            #
##############################################################




################################################################################################################################
#                                              8. REDUCTION IN MORTALITY RISK                                                  #
################################################################################################################################




################################################################################################################################
#                                                     9. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Median DALY prevented by walking in 2019 according to age group
plot_daly_prevented <- ggplot(burden_per_age, aes(x = age_grp10, y = tot_daly, fill = disease)) +
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

# Tables of HIA outcomes
export(burden, here("output", "Tables", "2019", "Resampling", "HIA_per_disease.xlsx"))
export(Rubin_burden, here("output", "Tables", "2019", "Resampling", "HIA_per_disease_Rubin.xlsx"))
export(burden_per_age, here("output", "Tables", "2019", "Resampling", "HIA_per_age.xlsx"))
export(Rubin_burden_per_age, here("output", "Tables", "2019", "Resampling", "HIA_per_age_Rubin.xlsx"))

       
# Plot
ggsave(here("output", "Plots", "2019", "DALY_prevented.png"), plot = plot_daly_prevented)






