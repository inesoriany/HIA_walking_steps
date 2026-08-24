#################################################
########     SENSITIVITY ANALYSES         #######
#################################################

# Modified parameters - files outputted :
  # Scenario 1 - Former DRF
  # Scenario 2 - Walking speed 5.3 km/h
  # Scenario 3 - Age limit < 75 years
  


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
  scales,
  progress      # Progression bar
)

################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Walkers dataset
emp_long <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Risk reduction distribution table
reduction_risk_distrib_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))

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

# Outcome considered
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")



################################################################################################################################
################################################################################################################################
#                                                       ALTERNATIVE DRF                                                        #
################################################################################################################################
################################################################################################################################

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
      disease = dis,
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
  mutate(week_time = 7*(nbkm_tot_walking + baseline_step*step_length)*60/walk_speed)  %>%       # baseline at 2000 steps 
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
  N = 1000
)
 

# Export : Table of HIA outcomes per simulation
export(alt_burden_total, here("output", "RDS", "Sensitivity analyses", "HIA_DRF_sensitivity_1000replicate.rds"))



################################################################################################################################
#                                            6. HEALTH IMPACT ASSESSMENT - baseline                                            #
################################################################################################################################

##############################################################
#                      PER DISEASES                          #
##############################################################

# Total of prevented burden of each disease for each simulation 
set.seed(123)
baseline_burden_total <- HIA_burden_total(
  data_list = baseline_replicate_list,
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
  N = 1000
)
 

# Export : Table of HIA outcomes per simulation
export(baseline_burden_total, here("output", "RDS", "Sensitivity analyses", "HIA_DRF_baseline_sensitivity_1000replicate.rds"))



################################################################################################################################
#                         7. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

# Import data
alt_burden_total <- import(here("output", "RDS", "Sensitivity analyses", "HIA_DRF_sensitivity_1000replicate.rds"))
baseline_burden_total <- import(here("output", "RDS", "Sensitivity analyses", "HIA_DRF_baseline_sensitivity_1000replicate.rds"))


##############################################################
#                       PER DISEASES                         #
##############################################################

# 2019 burden (compared to a counterfactual scenario of "zero walking")
set.seed(123)
ALT_burden_per_disease <- HIA_burden_IC(alt_burden_total, alt_dis_vec, outcome_vec, calc_replicate_IC) 

# Baseline burden
set.seed(123)
BASE_burden_per_disease <- HIA_burden_IC(baseline_burden_total, alt_dis_vec, outcome_vec, calc_replicate_IC)


# Additional gains of the modal shift scenario
DRF_burden_per_disease <- ALT_burden_per_disease %>% 
    mutate(tot_cases = tot_cases - BASE_burden_per_disease[["tot_cases"]],
           tot_cases_low = tot_cases_low - BASE_burden_per_disease[["tot_cases_low"]], 
           tot_cases_up = tot_cases_up - BASE_burden_per_disease[["tot_cases_up"]],
           tot_daly = tot_daly - BASE_burden_per_disease[["tot_daly"]],
           tot_daly_low = tot_daly_low - BASE_burden_per_disease[["tot_daly_low"]],
           tot_daly_up = tot_daly_up - BASE_burden_per_disease[["tot_daly_up"]])  %>% 
    select(disease, tot_cases, tot_cases_low, tot_cases_up, tot_daly, tot_daly_low, tot_daly_up)

# Total for morbidity
  DRF_burden_morbidity <- DRF_burden %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  
# Total for all diseases
  DRF_burden_global <- DRF_burden %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 

# Gather results
  DRF_burden <- bind_rows(DRF_burden_per_disease, DRF_burden_morbidity, DRF_burden_global)

################################################################################################################################
#                                                      8. EXPORT DATA                                                          #
################################################################################################################################
export(DRF_burden, here("output", "Tables", "Sensitivity analyses", "HIA_DRF_sensitivity_1000replicate.xlsx"))





################################################################################################################################
################################################################################################################################
#                                                        WALKING SPEED                                                         #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                                      1. PARAMETERS                                                           #
################################################################################################################################

# Change of walking speed (Compendium, 2011)
walk_speed <- 5.3  # km/h

# Diseases considered
dis_vec <- c("mort", "cvd", "diab2", "dem", "dep")


################################################################################################################################
#                                                    2. DATA PREPARATION                                                       #
################################################################################################################################
emp_walk_speed <- emp_long  %>% 
  mutate(nbkm_intermodal_walk = intermodal_walk_time * walk_speed / 60,
         nbkm_tot_walking = nbkm_main_walk + nbkm_intermodal_walk,
         step_commute = nbkm_tot_walking / step_length,
         step = pmin(12000, round(step_commute / 100) * 100 + baseline_step))


# EMP Dataset per disease
speed_replicate_list <- list() 

for (dis in dis_vec) {
  speed_replicate_list[[dis]] <- emp_walk_speed %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                 3. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                      PER DISEASES                          #
##############################################################
# Total of prevented burden of each disease for each simulation 
set.seed(123)
speed_burden_total <- HIA_burden_total(speed_replicate_list, calc_HIA_replicate, incidence_distrib_table, dep_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 prop_relapse, duration_recovery, vsl, 
                                 group = NULL,
                                 N = 1000)
 

# Export : Table of HIA outcomes per simulation
export(speed_burden_total, here("output", "RDS", "Sensitivity analyses", "HIA_speed_sensitivity_1000replicate.rds"))



################################################################################################################################
#                         4. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
##############################################################

# Import data
speed_burden_total <- import(here("output", "RDS", "Sensitivity analyses", "HIA_speed_sensitivity_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  speed_burden_per_disease <- HIA_burden_IC(speed_burden_total, dis_vec, outcome_vec, calc_replicate_IC) 


  # Total for morbidity
  speed_burden_morbidity <- speed_burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  
  # Total for all diseases
  speed_burden_global <- speed_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  speed_burden <- bind_rows(speed_burden_per_disease, speed_burden_morbidity, speed_burden_global)



################################################################################################################################
#                                                      5. EXPORT DATA                                                          #
################################################################################################################################
export(speed_burden, here("output", "Tables", "Sensitivity analyses", "HIA_speed_sensitivity_1000replicate.xlsx")) 




################################################################################################################################
################################################################################################################################
#                                                          AGE LIMIT                                                           #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                                    1. DATA PREPARATION                                                       #
################################################################################################################################

emp_age <- emp_long  %>% 
  mutate(step = pmin(12000, round(step_commute / 100) * 100 + baseline_step))  %>% 
  filter(age < 75)       # age limit of 75 years above which physical activity doesn’t decrease the mortality risk


# EMP Dataset per disease
age_replicate_list <- list() 

for (dis in dis_vec) {
  age_replicate_list[[dis]] <- emp_age %>% 
    filter(disease == dis)
}



################################################################################################################################
#                                                 2. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

##############################################################
#                      PER DISEASES                          #
##############################################################
# Total of prevented burden of each disease for each simulation 
set.seed(123)
age_burden_total <- HIA_burden_total(age_replicate_list, calc_HIA_replicate, incidence_distrib_table, dep_distrib_table, reduction_risk_distrib_table, dw_distrib_table, 
                                 dis_vec, 
                                 prop_relapse, duration_recovery, vsl, 
                                 group = NULL,
                                 N = 1000)
 

# Export : Table of HIA outcomes per simulation
export(age_burden_total, here("output", "RDS", "Sensitivity analyses", "HIA_age_sensitivity_1000replicate.rds"))



################################################################################################################################
#                         3. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
##############################################################

# Import data
age_burden_total <- import(here("output", "RDS", "Sensitivity analyses", "HIA_age_sensitivity_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  age_burden_per_disease <- HIA_burden_IC(age_burden_total, dis_vec, outcome_vec, calc_replicate_IC) 


  # Total for morbidity
  age_burden_morbidity <- age_burden_per_disease %>%
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  
  # Total for all diseases
  age_burden_global <- age_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  age_burden <- bind_rows(age_burden_per_disease, age_burden_morbidity, age_burden_global)



################################################################################################################################
#                                                      4. EXPORT DATA                                                          #
################################################################################################################################
export(age_burden, here("output", "Tables", "Sensitivity analyses", "HIA_age_sensitivity_1000replicate.xlsx"))




################################################################################################################################
################################################################################################################################
#                                              SENSITIVITY ANALYSES - Figure                                                   #
################################################################################################################################
################################################################################################################################

################################################################################################################################
#                                                      1. IMPORT DATA                                                          #
################################################################################################################################

# Main analysis
HIA_main <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))  %>% 
  filter(disease != "Morbidity")  %>% 
  mutate(analysis = "main") %>% 
  select(analysis, disease, tot_daly, tot_daly_low, tot_daly_up)


# Alternative DRF
HIA_DRF <- import(here("output", "Tables", "Sensitivity analyses", "HIA_DRF_sensitivity_1000replicate.xlsx")) %>% 
  filter(disease != "Morbidity") %>% 
  mutate(analysis = "sc1") %>% 
  select(analysis, disease, tot_daly, tot_daly_low, tot_daly_up) %>% 
  left_join(HIA_main %>% select(disease, tot_daly_ref = tot_daly), by = "disease") %>% 
  mutate(relative_perc_daly = (tot_daly - tot_daly_ref) * 100 / tot_daly_ref) %>%               # Relative percentage change in total DALYs
  select(-tot_daly_ref)


# Walking speed 5.3 km/h
HIA_speed <- import(here("output", "Tables", "Sensitivity analyses", "HIA_speed_sensitivity_1000replicate.xlsx")) %>% 
  filter(disease == "All")  %>% 
  mutate(analysis = "sc2") %>% 
  select(analysis, disease, tot_daly, tot_daly_low, tot_daly_up) %>% 
  mutate(relative_perc_daly = (tot_daly - (HIA_main %>% filter(disease == "All") %>% pull(tot_daly))) *100 /
                              (HIA_main %>% filter(disease == "All") %>% pull(tot_daly)))      # Relative percentage change in total DALYs


# Age limit < 75 years
HIA_age <- import(here("output", "Tables", "Sensitivity analyses", "HIA_age_sensitivity_1000replicate.xlsx")) %>% 
  filter(disease == "All") %>%
  mutate(analysis = "sc3") %>%
  select(analysis, disease, tot_daly, tot_daly_low, tot_daly_up)%>% 
  mutate(relative_perc_daly = (tot_daly - (HIA_main %>% filter(disease == "All") %>% pull(tot_daly))) *100 /
                              (HIA_main %>% filter(disease == "All") %>% pull(tot_daly)))     # Relative percentage change in total DALYs



################################################################################################################################
#                                                    2. DATA PREPARATION                                                       #
################################################################################################################################
SENSI_data <- bind_rows(HIA_main, HIA_DRF, HIA_speed, HIA_age)  %>% 
  mutate(outcome = "tot_daly")



################################################################################################################################
#                                                     3. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Sensitivity analysis of the number health events prevented under different assumptions

labels_analysis <- c(
  "main" = "Main Analysis",
  "sc1"  = "Alternative DRF",
  "sc2" = "Walking speed 5.3 km/h",
  "sc3" = "Age limit < 75 years"
)


plot_sensi <- ggplot(SENSI_data) +
  geom_point(aes(x = tot_daly,
                 y = disease,
                 color = factor(disease)),
             size = 2.5,
             shape = 18) +
  geom_segment(aes(x = tot_daly_low,
                   xend = tot_daly_up,
                   y = disease,
                   yend = disease,
                   color = factor(disease),
                   group = disease),
               linewidth = 0.6) +
  facet_grid(
    rows = vars(analysis),
    cols = vars(outcome),
    labeller = labeller(
      analysis = labels_analysis,
      outcome = c("tot_daly" = "Prevented DALYs")
    ),
    scales = "free_y",
    space = "free"
  ) +
  scale_color_manual(
    values = colors_disease,
    labels = names_disease
  ) +
  labs(
    title = ,
    y = NULL,
    x = NULL,
    color = "Disease"
  ) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 8),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.6, "lines")
  ) +
  scale_y_discrete(labels = names_disease) +
  scale_x_continuous(labels = scales::label_comma()) +
  geom_vline(data = data.frame(outcome = c("tot_daly"),
                               xintercept = c(50000)),
             aes(xintercept = xintercept),
             linetype = "dashed", color = "black", linewidth = 0.3) +
  guides(color = guide_legend(title = NULL))


plot_sensi 



################################################################################################################################
#                                                        2. EXPORT                                                             #
################################################################################################################################
ggsave(here("output", "Plots", "Sensitivity analyses", "plot_sensitivity.png"), plot = plot_sensi)
export(SENSI_data, here("output", "Tables", "Sensitivity analyses", "sensitivity_analysis.xlsx"))
