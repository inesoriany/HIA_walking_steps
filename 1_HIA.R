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
emp_walk <- import(here("data_clean", "EMP_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_table <- import(here("data_clean", "DRF", "rr_central_interpolated.rds"))


# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "diab2", "dem")



################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
health_walkers <- emp_walk %>% 
  # Round the number of steps to the nearest ten
  mutate(step = pmin(12000, round(step_commute / 10) * 10 + 2000))


# Associate simulated RR
  # Central
health_walkers_mid <- health_walkers
  for (dis in dis_vec) {
    health_walkers_mid <- health_walkers_mid %>% 
      left_join(
        rr_table %>% 
          filter(disease == dis) %>% 
          select(step, mid) %>% 
          rename(!!paste0(dis, "_rr") := mid),
        by = "step" 
      )
  }

  # Lower bound
health_walkers_low <- health_walkers
  for (dis in dis_vec) {
    health_walkers_low <- health_walkers_low %>% 
      left_join(
        rr_table %>% 
          filter(disease == dis) %>% 
          select(step, low) %>% 
          # To calculate the upper bound of reduction of the relative risk, use RR lower bound because the decrease will be higher,
          # i.e. the person exposed (walking) is less likely to have the disease 
          rename(!!paste0(dis, "_rr") := up),
        by = "step"
      )
  }

  # Upper bound
health_walkers_up <- health_walkers
  for (dis in dis_vec) {
    health_walkers_up <- health_walkers_up %>% 
      left_join(
        rr_table %>% 
          filter(disease == dis) %>% 
          select(step, up) %>% 
          rename(!!paste0(dis, "_rr") := low),
        by = "step"
      )
  }



################################################################################################################################
#                                                 5. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

health_walkers_mid <- calc_HIA (health_walkers_mid, "mid", rr_table, dis_vec)

health_walkers_low <- calc_HIA (health_walkers_low, "low", rr_table, dis_vec)

health_walkers_up <- calc_HIA (health_walkers_up, "up", rr_table, dis_vec)




##############################################################
#                      HIA OUTCOMES                          #     with cases, DALY and medical costs
##############################################################

## Survey design ----
surv_dis <- health_walkers_mid %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_indc,
                   strata = c(sexe, age_grp.x, quartile_rev),           # by sex and age group
                   nest = TRUE)
  # IC
    # Upper bound
    surv_dis_up <- health_walkers_up %>% 
      as_survey_design(ids = ident_ind,
                       weights = pond_indc,
                       strata = c(sexe, age_grp.x, quartile_rev),           
                       nest = TRUE)
    
    # Lower bound
    surv_dis_low <- health_walkers_low %>% 
      as_survey_design(ids = ident_ind,
                       weights = pond_indc,
                       strata = c(sexe, age_grp.x, quartile_rev),      
                       nest = TRUE)


## Total of prevented cases, DALY and saved costs, for each disease in 2019
burden <- data.frame()
for (dis in dis_vec) {
  dis_burden <- burden_prevented(surv_dis, dis, NULL)
  burden <- bind_rows(burden, dis_burden)   
}
    
  # IC
    # Upper bound
    burden_up <- data.frame()
    burden_low <- data.frame()
    for(dis in dis_vec) {
      dis_burden_up <- burden_prevented(surv_dis_up, dis, NULL)
      burden_up <- bind_rows(burden_up, dis_burden_up) 
      
      # Lower bound
      dis_burden_low <- burden_prevented(surv_dis_low, dis, NULL)
      burden_low <- bind_rows(burden_low, dis_burden_low) 
    }
    
    
  # Gather results with IC
  burden_IC <- burden %>% 
    mutate(cases_low = burden_low[,"tot_cases"], cases_up = burden_up[,"tot_cases"],
           daly_low = burden_low[,"tot_daly"], daly_up = burden_up[,"tot_daly"],
           medic_costs_low = burden_low[,"tot_medic_costs"], medic_costs_up = burden_up[,"tot_medic_costs"])


  
  
##############################################################
#                    ECONOMIC IMPACT (2)                     #
##############################################################
  
  
## SOCIAL COSTS (intangible)----
  # Add social costs
  burden_IC <- burden_IC %>% 
    mutate(tot_soc_costs = tot_daly*vsl,
           soc_costs_low = daly_low*vsl,
           soc_costs_up = daly_up*vsl)
  
  # Reorganize columns
  burden_IC <- burden_IC %>% 
    select(disease,
           tot_cases, tot_cases_se, cases_low, cases_up,
           tot_daly, tot_daly_se, daly_low, daly_up,
           tot_medic_costs, tot_medic_costs_se, medic_costs_low, medic_costs_up,
           tot_soc_costs, soc_costs_low, soc_costs_up) %>% 
    mutate(disease = recode_factor(disease, 
                                   bc = "Breast cancer", 
                                   cc="Colon cancer" , 
                                   cvd ="CVD" , 
                                   dem ="Dementia",
                                   diab2 ="T2 Diabetes" , 
                                   dep = "Depression",
                                   mort ="Mortality")) 
  
    






