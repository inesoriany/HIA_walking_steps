#################################################
########     SENSITIVITY ANALYSIS         #######
#################################################

# Modified parameters - files outputted :
  # Walking speed 5.3 km/h
  # Former DRF


################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  purrr,        # Loop
  stringr,      # Text extraction
  srvyr,        # Survey
  survey,
  ggplot2,      # Data visualization
  progress      # Progression bar
)


################################################################################################################################
################################################################################################################################
#                                                       ALTERNATIVE DRF                                                        #
################################################################################################################################
################################################################################################################################

# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Depression duration distribution table
dep_distrib_table <- import(here("data_clean", "Diseases", "dep_duration_distrib_table.xlsx"))

# Disability weights distribution table
dw_distrib_table <- import(here("data_clean", "Diseases", "dw_distrib_table.xlsx"))


# Import functions
source(here("R_code", "Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("R_code", "Parameters.R"))

# Diseases considered
alt_dis_vec <- c("mort", "bc", "cc", "cvd", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly")



################################################################################################################################
#                                                3. RELATIVE RISKS DISTRIBUTION                                                #
################################################################################################################################
rr_women <- data.frame()
rr_men <- data.frame()
alt_rr_distrib <- data.frame()
set.seed(123)
for (dis in alt_dis_vec){
  params <- dis_setting(dis)

  # Generate RR normal distributions
  rr_women <-  data.frame(
    disease = dis,
    sex = "Female", 
    ref = params$ref_women) %>%
  mutate(rr_distrib =list(generate_RR_distrib (RR = params$rr_women, params$rr_women_lb, params$rr_women_ub, N=1000)))  %>% 
  mutate(rr_distrib = list(sort(unlist(rr_distrib)))) %>%
  ungroup()
  if (dis %in% c("bc", "cc")) {                                                          # if disease is bc or cc
    rr_men <- data.frame(
      dis = dis,
      sex = "Male", 
      ref = params$ref_men) %>%
    mutate(disease = dis) %>%
    mutate(sex = "Male") %>%
    mutate(rr_distrib =list(generate_RR_distrib (RR = params$rr_men, params$rr_men_lb, params$rr_men_ub, N=1000)))  %>% 
    mutate(rr_distrib = list(sort(unlist(rr_distrib)))) %>%

    ungroup()
  } else {
    rr_men <- rr_women  %>% 
    mutate(sex = "Male")
  }
alt_rr_distrib <- bind_rows(alt_rr_distrib, rr_women, rr_men)
}


# Separate the rr_distrib column into multiple columns
alt_rr_distrib_table <- alt_rr_distrib  %>% 
  unnest_wider(rr_distrib, names_sep = "_") %>% 
  pivot_longer(
    cols = starts_with("rr_distrib_"), 
    names_to = "simulation_id",                     # column name for the simulation ID
    values_to = "simulated_rr")  %>%         # column name for the simulated incidence values
  mutate(simulation_id = as.numeric(str_remove(simulation_id, "rr_distrib_")))      # simulation ID as a numeric value




################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_walk <- emp_long %>% 
  mutate(week_time = 7*(nbkm_tot_walking + baseline_step*step_length)*60/walk_speed)  %>%       # baseline at 2000 steps  %>% 
  left_join(alt_rr_distrib_table  %>% select(disease, sex, ref)  %>% distinct(disease, sex, .keep_all = TRUE),, by = c("disease" , "sex"))

emp_baseline <- emp_long %>% 
  mutate(week_time = 7*(baseline_step*step_length)*60/walk_speed)  %>% 
  left_join(alt_rr_distrib_table  %>% select(disease, sex, ref)  %>% distinct(disease, sex, .keep_all = TRUE), by = c("disease" , "sex"))



# EMP Dataset per disease
alt_replicate_list <- list() 
for (dis in alt_dis_vec) {
  alt_replicate_list[[dis]] <- emp_walk %>% 
    filter(disease == dis)
}

baseline_replicate_list <- list() 
for (dis in alt_dis_vec) {
  baseline_replicate_list[[dis]] <- emp_baseline %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                 5. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################
##############################################################
#                      PER DISEASES                          #
##############################################################

alt_replicate_list <- calc_alt_HIA(alt_replicate_list, incidence_distrib_table, dep_distrib_table, alt_rr_distrib_table, dw_distrib_table, alt_dis_vec, prop_relapse, duration_recovery, vsl)
# ça fonctionne

# Total of prevented burden of each disease for each simulation 
set.seed(123)
alt_burden_total <- HIA_burden_total(
  data_list = alt_replicate_list,
  function_calc_HIA = calc_alt_HIA,
  incidence_distrib_table = incidence_distrib_table,
  dep_duration_table = dep_distrib_table,
  reduction_risk_distrib_table = alt_rr_distrib_table,
  dw_distrib_table = dw_distrib_table,
  dis_vec = alt_dis_vec,
  prop_relapse = prop_relapse,
  duration_recovery = duration_recovery,
  vsl = vsl,
  group = NULL,
  N = 1
)
 

# Export : Table of HIA outcomes per simulation
export(alt_burden_total, here("output", "RDS", "2019", "Monte Carlo", "HIA_1000replicate.rds"))


################################################################################################################################
#                                            6. HEALTH IMPACT ASSESSMENT - baseline                                            #
################################################################################################################################

