#################################################
##############     MODAL SHIFT     ##############
#################################################



################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  srvyr,        # Survey
  tidyr,        # Table - Data organization, extraction
  tidyverse,    # Data manipulation and visualization
  ggplot2,      # Plotting
  patchwork,    # Plots combining
  viridis       # Color palette
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Import drivers dataset
emp_drive <- import(here("data_clean", "EMP_dis_drivers.xlsx"))

# RR by step, simulated dose-response relationships
rr_distrib_table <- import(here("data_clean", "Diseases", "DRF", "rr_sim_interpolated.rds"))

# Disability weights distribution
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

# Driven distances shifted to walked (km)
dist_vec <- c(0.5, 1, 1.5, 2)

# Percentages of car trips shifted to walking
perc_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5)




################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################
# Initialization
emp_drive <- emp_drive %>% 
  # Day step and baseline at 2000
  mutate(step = pmin(12000,round(day_step_commute_shift / 100) * 100 + baseline_step))



################################################################################################################################
#                                              5. HEALTH IMPACT ASSESSMENT                                                     #
################################################################################################################################
set.seed(123)
burden_tot <- data.frame()

for(dist in dist_vec){
  print(paste0("Distance = ", dist))
  
  drivers_dist <- emp_drive %>% 
    filter(!is.na(mdisttot_fin1) & mdisttot_fin1 <= dist)
  
  for(dis in dis_vec){
    
    for(perc in perc_vec){
      print(paste0("Disease = ", dis, ", Share = ", perc))
      N <- 100
      burden_run <- list()
      
      for(i in 1:N){
        print(paste0("Run = ", i))
        
        ## ----------------------------------------
        ## 1. Scenario sample
        ## ----------------------------------------
        burden_sample <- drivers_dist %>%
          filter(disease == dis) %>%
          slice_sample(prop = perc) %>%
          mutate(
            run = i,
            percentage = perc,
            distance = dist,
            disease = dis
          )
        
        ## ----------------------------------------
        ## 2. Draw RR
        ## ----------------------------------------
        rr_disease <- if (dis == "bc") "cancer" else dis
        
        sim_id <- rr_distrib_table %>%
          filter(disease == rr_disease) %>%
          pull(simulation_id) %>%
          unique() %>%
          sample(1)
        
        rr_baseline <- rr_distrib_table %>%
          filter(
            disease == rr_disease,
            step == baseline_step,
            simulation_id == sim_id
          ) %>%
          pull(rr_interpolated)
        
        rr_step <- rr_distrib_table %>%
          filter(
            disease == rr_disease,
            simulation_id == sim_id
          ) %>%
          select(step, rr_interpolated) %>%
          mutate(
            reduction_risk = (rr_baseline - rr_interpolated) / rr_baseline
          )
        
        burden_sample <- burden_sample %>%
          left_join(
            rr_step %>% select(step, reduction_risk),
            by = "step"
          ) %>%
          mutate(rr_simulation_id = sim_id)
        
        ## ----------------------------------------
        ## 3. Draw DW
        ## ----------------------------------------
        dw_draw <- dw_distrib_table %>%
          filter(disease == dis) %>%
          pull(simulated_dw) %>%
          sample(1)
        
        burden_sample <- burden_sample %>%
          mutate(dw = dw_draw)
        
        ## ----------------------------------------
        ## 4. Cases, DALY & costs (individual)
        ## ----------------------------------------
        burden_sample <- reduc_incidence(burden_sample)
        burden_sample <- daly(burden_sample)
        burden_sample <- medic_costs(burden_sample, dis)
        
        burden_sample <- burden_sample %>%
          mutate(soc_costs = daly * vsl)
        
        ## ----------------------------------------
        ## 5. Aggregate burden by run
        ## ----------------------------------------
        surv_run <- burden_sample %>%
          as_survey_design(
            ids = ident_ind,
            weights = pond_indc
          )
        
        burden_run[[i]] <- surv_run %>%
          summarise(
            tot_cases = survey_total(cases, na.rm = TRUE),
            tot_daly = survey_total(daly, na.rm = TRUE),
            tot_medic_costs = survey_total(medic_costs, na.rm = TRUE),
            tot_soc_costs = survey_total(soc_costs, na.rm = TRUE)
          ) %>%
          mutate(
            run = i,
            disease = dis,
            distance = dist,
            percentage = perc
          )
      }
      
      ## ----------------------------------------
      ## 6. Combine runs
      ## ----------------------------------------
      burden_tot <- bind_rows(burden_tot, bind_rows(burden_run))
    }
  }
}


# Export HIA outcomes of 100 replications for each scenario of modal shifts
export(burden_tot, here("output", "RDS", "Modal shift", "HIA_modal_shift_100replicate.rds"))




################################################################################################################################
#                                  5. TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                               #
################################################################################################################################

# Import HIA outcomes of 100 replications for each scenario of modal shifts
burden_tot <- import(here("output", "RDS", "Modal shift",  "HIA_modal_shift_100replicate.rds"))



##############################################################
#                          GLOBAL                            #
##############################################################

global_shift <- burden_tot %>% 
  group_by(distance, percentage) %>% 
  summarise(tot_cases = sum(tot_cases),
            tot_cases_se = sum(tot_cases_se),
            tot_daly = sum(tot_daly),
            tot_daly_se = sum(tot_daly_se),
            tot_medic_costs = sum(tot_medic_costs) * 1e-6,                                # in million €
            tot_medic_costs_se = sum(tot_medic_costs_se) * 1e-6,
            tot_soc_costs = sum(tot_soc_costs * 1e-6),
            tot_soc_costs_se = sum(tot_soc_costs_se * 1e-6))



# IC per scenario
set.seed(123)
global_shift_IC <- data.frame()

for (dist in dist_vec) {
  for (perc in perc_vec) {
    scenario <- global_shift %>% 
      filter(distance == dist & percentage == perc)
    cases_IC <- as.data.frame(t(calc_replicate_IC(scenario, "tot_cases"))) %>% 
      rename(cases_low = "2.5%", cases_mid = "50%", cases_sup = "97.5%") %>% 
      mutate(distance = dist, percentage = perc)
    
    daly_IC <- as.data.frame(t(calc_replicate_IC(scenario, "tot_daly"))) %>% 
      rename(daly_low = "2.5%", daly_mid = "50%", daly_sup = "97.5%") %>% 
      mutate(distance = dist, percentage = perc)
    
    medic_costs_IC <- as.data.frame(t(calc_replicate_IC(scenario, "tot_medic_costs"))) %>% 
      rename(medic_costs_low = "2.5%", medic_costs_mid = "50%", medic_costs_sup = "97.5%") %>% 
      mutate(distance = dist, percentage = perc)
    
    soc_costs_IC <- as.data.frame(t(calc_replicate_IC(scenario, "tot_soc_costs"))) %>% 
      rename(soc_costs_low = "2.5%", soc_costs_mid = "50%", soc_costs_sup = "97.5%") %>% 
      mutate(distance = dist, percentage = perc)
    
    scenario_IC <- left_join(cases_IC, daly_IC, by = c("distance", "percentage"))
    scenario_IC <- left_join(scenario_IC, medic_costs_IC, by = c("distance", "percentage"))
    scenario_IC <- left_join(scenario_IC, soc_costs_IC, by = c("distance", "percentage"))
    global_shift_IC <- bind_rows(global_shift_IC, scenario_IC)
  }
}


# Export HIA for each scenario
export(global_shift_IC, here("output", "Tables", "Modal shift", "HIA_modal_shift_100replicate.xlsx"))






