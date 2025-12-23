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
        
        ## ---------------------------
        ## 1. Échantillon individuel
        ## ---------------------------
        burden_sample <- drivers_dist %>%
          filter(disease == dis) %>%
          slice_sample(prop = perc) %>%
          mutate(
            run = i,
            percentage = perc,
            distance = dist,
            disease = dis
          )
        
        ## ---------------------------
        ## 2. Tirage RR
        ## ---------------------------
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
        
        ## ---------------------------
        ## 3. Tirage DW
        ## ---------------------------
        dw_draw <- dw_distrib_table %>%
          filter(disease == dis) %>%
          pull(simulated_dw) %>%
          sample(1)
        
        burden_sample <- burden_sample %>%
          mutate(dw = dw_draw)
        
        ## ---------------------------
        ## 4. Cases, DALY & coûts
        ## ---------------------------
        burden_sample <- reduc_incidence(burden_sample)
        burden_sample <- daly(burden_sample)
        burden_sample <- medic_costs(burden_sample, dis)
        
        burden_sample <- burden_sample %>%
          mutate(soc_costs = daly * vsl)
        
        ## ---------------------------
        ## 5. AGRÉGATION PAR RUN
        ## ---------------------------
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
      
      ## ---------------------------
      ## 6. Combiner les runs
      ## ---------------------------
      burden_tot <- bind_rows(burden_tot, bind_rows(burden_run))
    }
  }
}




# Export HIA outcomes of 100 replications for each scenario of modal shifts
export(burden_tot, here("output", "RDS", "Modal shift", "HIA_modal_shift_100replicate.rds"))




set.seed(123)
burden_tot <- data.frame()

for (dist in dist_vec) {
  
  message("Distance = ", dist)
  
  drivers_dist <- emp_drive %>%
    filter(!is.na(step),
           step <= (round((dist / step_length + baseline_step) / 100) * 100))
  
  for (perc in perc_vec) {
    
    message("  Share = ", perc)
    
    for (i in 1:100) {
      
      message("    Run = ", i)
      
      # 1. Select individuals for the scenario
      shifted_list <- vector("list", length(dis_vec))
      names(shifted_list) <- dis_vec
      
      for (dis in dis_vec) {
        
        dis_drivers <- drivers_dist %>% 
          filter(disease == dis)
        
        shifted_list[[dis]] <- dis_drivers %>% slice_sample(prop = perc)
        
      }
      
      # 2. Calculate HIA
      HIA_shifted_list <- calc_HIA_shift(data_list = shifted_list,
                                         rr_distrib_table = rr_distrib_table,
                                         dw_distrib_table = dw_distrib_table,
                                         dis_vec = dis_vec,
                                         vsl = vsl,
                                         baseline_step = baseline_step)
      
      # 3. Add simulation_id for compatibility
      for (dis in dis_vec) {
        HIA_shifted_list[[dis]] <- HIA_shifted_list[[dis]] %>%
          mutate(simulation_id = i)
      }
      
      # 4. Burden prevented
      burden_i <- burden_replicate_prevented(
        data_list = HIA_shifted_list,
        dis_vec = dis_vec,
        group = NULL,
        N = 1
      ) %>%
        mutate(run = i,
               percentage = perc,
               distance = dist)
      
      burden_tot <- bind_rows(burden_tot, burden_i)
    }
  }
}
