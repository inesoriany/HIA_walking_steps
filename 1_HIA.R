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
  ggplot2       # Data visualization
)




################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################
# Walkers dataset
emp_walk <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_central_table <- import(here("data_clean", "DRF", "rr_central_interpolated.rds"))


# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "bc", "cancer", "diab2", "dem", "dep")

# Bound
bound_vec <- c("mid", "low", "up")


################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_walk <- emp_walk %>% 
  # Round the number of steps to the nearest hundred and baselin at 2000
  mutate(step = pmin(12000, round(step_commute / 10) * 10 + 2000))


# EMP Dataswet per disease
for (dis in dis_vec) {
  for (bound in bound_vec) {
    assign(
      paste0(dis, "_walkers_", bound),
      emp_walk %>% filter(disease == dis))
  }
}


################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################

## --------------------------------------------
## ALL DISEASES (EXCEPT BREAST CANCER)
## --------------------------------------------

# Baseline step 2000
rr_baseline <- rr_central_table %>%
  filter(step == 2000) %>%
  select(disease, rr2000_mid = mid, rr2000_low = low, rr2000_up = up)

rr_central_table <- rr_central_table %>%
  left_join(rr_baseline, by = "disease")


# Reduction risk for all diseases except breast cancer
for (bound in bound_vec) {
  rr_central_table <- rr_central_table %>%
    mutate(!!paste0("reduction_risk_", bound) := 
             (.data[[paste0("rr2000_", bound)]] - .data[[bound]]) /.data[[paste0("rr2000_", bound)]]) 
}
  # Rename column: To calculate the upper bound of reduction of the relative risk, use RR lower bound because the decrease will be higher i.e. the person exposed (walking) is less likely to have the disease 
  rr_central_table <- rr_central_table %>% 
    rename(reduction_risk_low = reduction_risk_up,
           reduction_risk_up = reduction_risk_low)





# Associate disease risk reductions to individuals
for (dis in dis_vec) {
  for (bound in bound_vec) {
    walkers_name_bound <- paste0(dis, "_walkers_", bound)
    dis_walkers_bound <- get(walkers_name_bound)
    
    # --------------------------------------
    # BREAST CANCER (special case)
    # --------------------------------------
    if (dis == "bc") {
      
      params <- dis_setting(dis)
      
      dis_walkers_bound <- log_reduction_risk(data = dis_walkers_bound,
                                              dis  = dis,
                                              rr_women = params[[paste0("rr_women_", bound)]],
                                              rr_men   = params[[paste0("rr_men_",   bound)]],
                                              ref_women = params$ref_women,
                                              ref_men   = params$ref_men,
                                              week_base = week_base)
      
      
      # --------------------------------------
      # OTHER DISEASES
      # --------------------------------------
    } else {
      
      dis_rr <- rr_central_table %>%
        filter(disease == dis) %>%
        select(step, reduction_risk =!!sym(paste0("reduction_risk_", bound)))
      
      dis_walkers_bound <- dis_walkers_bound %>%
        left_join(dis_rr, by = "step")
    }
    
    assign(walkers_name_bound, dis_walkers_bound)
  }
}




##############################################################
#                      CASES PREVENTED                       #
##############################################################

# Calculate prevented cases
for (dis in dis_vec) {
  for (bound in bound_vec) {
    walkers_name_bound <- paste0(dis, "_walkers_", bound)
    dis_walkers_bound <- get(walkers_name_bound)
    
    dis_walkers_bound <- reduc_incidence(dis_walkers_bound)
    
    assign(walkers_name_bound, dis_walkers_bound)
  }
}


# Associate in a data.frame
for (bound in bound_vec) {
  
  health_walkers <- data.frame()
  
  for (dis in dis_vec) {
    
    walkers_name_bound <- paste0(dis, "_walkers_", bound)
    dis_walkers_bound <- get(walkers_name_bound)
    
    health_walkers <- bind_rows(health_walkers, dis_walkers_bound)
  }
  
  assign(paste0("health_walkers_", bound), health_walkers)
}






################################################################################################################################
#                                              5. TOTAL OF PREVENTED CASES                                                     #
################################################################################################################################

# Total of prevented cases per bound
for(bound in bound_vec) {
  cases_prev_name <- paste0("cases_prev_", bound)
  cases_prev_bound <- data.frame()
  
  health_walkers_bound <- get(paste0("health_walkers_", bound))
  
  # Survey design
  health_walkers_svy <- health_walkers_bound %>%
    as_survey_design(ids = ident_ind,
                     strata = c(age_grp10, sex),
                     weights = pond_indc,
                     nest = TRUE)
  
  for (dis in dis_vec) {
    dis_cases <- health_walkers_svy %>%
      filter(disease == dis) %>%
      group_by(age_grp10, sex) %>%
      summarise(tot_cases = survey_total(reduc_incidence, na.rm = TRUE)) %>% 
      mutate(disease = dis)
    
    # Results for all diseases
    cases_prev_bound <- bind_rows(cases_prev_bound, dis_cases)
  }
  
  # Results for all diseases per bound
  assign(cases_prev_name, cases_prev_bound)
}



# Gather results with IC
IC_cases_prev <- cases_prev_mid %>% 
  mutate(tot_cases_low = cases_prev_low[,"tot_cases"], 
         tot_cases_up = cases_prev_up[,"tot_cases"])




################################################################################################################################
#                                                     6. VISUALIZATION                                                         #
################################################################################################################################
  











################################################################################################################################
#                                                      7. EXPORT DATA                                                          #
################################################################################################################################
# Tables
export(IC_cases_prev, here("output", "Tables", "2019", "Prev_cases_2019.xlsx"))






