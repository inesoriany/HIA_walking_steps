#################################################
############ HEALTH IMPACT ASSESSMENT ###########
#################################################

# Files needed :
# EMP_walkers.xlsx
# rr_central_interpolated.rds
# 0_Functions.R
# 0_Parameters.R

# Files outputted :



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
  ggplot2       # Data visualization
)




################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################
# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_distrib_table <- import(here("data_clean", "DRF", "rr_sim_interpolated.rds"))


# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")
no_bc_vec = c("mort", "cvd", "cancer", "diab2", "dem", "dep")

# Bound 
bound_vec = c("mid", "low", "up")


################################################################################################################################
#                                                  4. SIMULATED RELATIVE RR                                                    #
################################################################################################################################

# Baseline step 2000
rr_baseline <- rr_distrib_table %>%
  filter(step == 2000) %>%
  select(disease, simulation_id, rr2000 = rr_interpolated)

# Relative RR
rr_distrib_table <- rr_distrib_table %>%
  left_join(rr_baseline, by = c("disease", "simulation_id")) %>%
  mutate(relative_rr = rr2000 / rr_interpolated)





################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_long <- emp_long %>% 
  # Round the number of steps to the nearest hundred
  mutate(step = pmin(12000, round(step_commute / 100) * 100 + 2000))


# Associate simulated RR
emp_rr_adjusted <- emp_long %>%
  inner_join(rr_distrib_table, by = c("disease", "step"), relationship = "many-to-many") %>%
  mutate(adjusted_rate = rate * relative_rr) %>%
  ungroup()



# EMP Dataset per disease
for (dis in dis_vec) {
  assign(
    paste0(dis, "_walkers"),
    emp_long %>% filter(disease == dis)
  )
}


################################################################################################################################
#                                                 5. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################
# Associate relative RR to each individual for each disease

  for (dis in no_bc_vec) {
    walkers_name <- paste0(dis, "_walkers")
    cases_name <- paste0(dis, "_cases")
    
    dis_rr <- rr_distrib_table %>%
      filter(disease == dis) %>%
      select(simulation_id, step, relative_rr)
    
    dis_cases <- get(walkers_name) %>%
      left_join(dis_rr, by = "step", relationship = "many-to-many")
    
    assign(cases_name, dis_cases)
  }


# For breast cancer
  params <- dis_setting("bc")
  bc_walkers <- log_reduction_risk(bc_walkers, "bc", 
                                   params$rr_women, params$rr_women_low, params$rr_women_up, 
                                   params$rr_men, params$rr_men_low, params$rr_men_up,
                                   params$ref_women, params$ref_men)



##############################################################
#                      CASES PREVENTED                       #
##############################################################

# Calculate prevented cases
  for (dis in no_bc_vec) {
    cases_name <- paste0(dis, "_cases")
      
    dis_cases <- get(cases_name) %>%
      mutate(reduc_incidence = relative_rr * rate)
      
    assign(cases_name, dis_cases)
  }
  

    
# Mean and 95% CI of prevented cases by individual
  for (dis in no_bc_vec) {
    walkers_name <- paste0(dis, "_walkers")
    cases_name <- paste0(dis, "_cases")
    cases_name_ic95 <- paste0(dis, "_cases_ic95")
    
    dis_cases_ic95 <- get(cases_name) %>%
      group_by(ID) %>%
      summarise(
        cases_mid = mean(reduc_incidence, na.rm = TRUE),
        cases_low = quantile(reduc_incidence, 0.025, na.rm = TRUE),
        cases_up  = quantile(reduc_incidence, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    assign(cases_name_ic95, dis_cases_ic95)
    
    # Add cases prevented in walkers dataset
    walkers_updated <- get(walkers_name) %>% 
      left_join(get(cases_name_ic95), by = "ID")
    
    assign(walkers_name, walkers_updated)
  }
  

  
# For breast cancer
  bc_walkers <- reduc_incidence(bc_walkers)
  

    

##############################################################     
#                           DALY                             #
##############################################################
for (dis in dis_vec) {
  walkers_name <- paste0(dis, "_walkers")
  dis_walkers <- get(walkers_name) 

  for (bound in bound_vec) {
    dis_walkers <- daly(dis_walkers, dis, bound)
  }

  assign(walkers_name, dis_walkers)
}

  
  

##############################################################
#                ECONOMIC IMPACT (MEDICAL)                   #
##############################################################

for (dis in dis_vec) {
  walkers_name <- paste0(dis, "_walkers")
  dis_walkers <- get(walkers_name) 

  for (bound in bound_vec) {
    dis_walkers <- medic_costs (dis_walkers, dis, bound)
  }

  assign(walkers_name, dis_walkers)
}




##############################################################
#       TOTAL OF PREVENTED CASES, DALY, MEDICAL COSTS        #
##############################################################


# Associate all the disease results to walkers dataset
health_walkers_list <- list()
  
  for (dis in dis_vec) {
    walkers_name <- paste0(dis, "_walkers")
    
    walkers_data <- get(walkers_name) %>%
      select(
        ID, disease,
        cases_mid, cases_low, cases_up,
        daly_mid, daly_low, daly_up,
        medic_costs_mid, medic_costs_low, medic_costs_up
      )
    
    # Walkers dataset for 1 disease
    emp_dis <- emp_long %>% filter(disease == dis)
    
    # Join
    emp_dis_joined <- emp_dis %>%
      left_join(walkers_data, by = c("ID", "disease"))
    
    health_walkers_list[[dis]] <- emp_dis_joined
  }
  
health_walkers <- bind_rows(health_walkers_list)
  



# SURVEY DESIGNS
surv_dis <- health_walkers %>% 
  as_survey_design(ids = ID,
                   weights = pond_indc,
                   strata = c(sex, age_grp10, quartile_rev, disease),           # by sex and age group
                   nest = TRUE)

burden <- data.frame()
for (dis in dis_vec) {
  dis_burden <- burden_prevented(surv_dis, dis, NULL)
  burden <- bind_rows(burden, dis_burden)   
}




