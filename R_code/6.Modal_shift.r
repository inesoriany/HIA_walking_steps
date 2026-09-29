#################################################
##############     MODAL SHIFT     ##############
#################################################





###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                         HIA - Modal shift of 50% of short car trips (<2km)                                                              #
###########################################################################################################################################################################
###########################################################################################################################################################################

################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  purrr,        # Loop
  srvyr,        # Survey
  tidyr,        # Table - Data organization, extraction
  tidyverse,    # Data manipulation and visualization
  ggplot2       # Plotting
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Drivers dataset
emp_car_trip <- import(here("data_clean", "EMP_dis_car_trips.xlsx"))

# Walkers dataset
emp_walkers <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Diseases duration distribution table
duration_distrib_table <- import(here("data_clean", "Diseases", "duration_distrib_table.xlsx"))

# Risk reduction distribution table
reduction_risk_distrib_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_distrib_table.xlsx"))

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
dis_vec <- c("mort", "cvd", "cancer", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")


# Modal shift scenario
  # Driven distance shifted to walked (km)
  dist <- 2 

  # Percentage of drivers shifted to walking
  perc <- 0.5

################################################################################################################################
#                                                      3. CO2 EMISSIONS                                                        #
################################################################################################################################
# Adjust CO2 emissions excluding intermodal walk
emp_car_trip <- emp_car_trip  %>% 
  mutate(co2_adjusted = ((nbkm_car - nbkm_intermodal_walk) * (co2_depl / nbkm_car)))


# Total CO2 emissions due to car for each individual 
tot_car_co2 <- emp_car_trip %>% 
  distinct(ident_ind, ident_dep, .keep_all = TRUE) %>% 
  group_by(ident_ind) %>% 
  summarise(co2_all_car = sum(co2_adjusted, na.rm = TRUE), .groups = "drop")


emp_car_trip <- emp_car_trip %>% 
  left_join(tot_car_co2, by = "ident_ind")


################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################
# Initialization
emp_short_trip <- emp_car_trip %>% 
  filter(!is.na(nbkm_car_jour) & nbkm_car_jour <= dist) %>% 
  mutate(co2_diminution = co2_all_car - co2_adjusted,
         co2_prop_reduction = 1 - co2_adjusted / co2_all_car)


# Conversion of short car trips to steps
step_shift_by_ind <- emp_short_trip %>% 
  distinct(ident_ind, ident_dep, .keep_all = TRUE) %>% 
  group_by(ident_ind) %>% 
  summarise(step_shift = sum(step_commute_jour, na.rm = TRUE), .groups = "drop")
  

emp_short_driver <- emp_short_trip %>% 
    left_join(step_shift_by_ind, by = "ident_ind")  %>% 
    mutate(short_car = TRUE)  %>% 
    group_by(ident_ind)  %>% 
    distinct(disease, .keep_all = TRUE)





################################################################################################################################
#                                              4. HEALTH IMPACT ASSESSMENT                                                     #
################################################################################################################################
set.seed(123)

N <- 1000

MODAL_burden_total <- data.frame()

for (i in 1:N) {
  print(paste0("Run ", i))
  
  sampled_ids <- emp_short_driver %>%
    distinct(ident_ind) %>%
    ungroup() %>%
    slice_sample(prop = perc) %>%
    pull(ident_ind)

  emp_driver_sample <- emp_short_driver %>%
      filter(ident_ind %in% sampled_ids)

  emp_walk_drive_sample <- emp_walkers  %>% 
      left_join(emp_driver_sample  %>% select(ident_ind, disease, step_shift, short_car), by = c("ident_ind", "disease"))  %>% 
      replace_na(list(step_shift = 0))  %>%
      mutate(step_total = step_commute_jour + step_shift,
             step = pmin(12000, round(step_total / 100) * 100 + baseline_step))     # Round the number of steps to the nearest hundred and baseline at 2000

  short_trip_list <- list()
  
  for (dis in dis_vec) {
    emp_dis_sample <- emp_walk_drive_sample %>%
      filter(disease == dis)

    short_trip_list[[dis]] <- emp_dis_sample
  }
  
  burden_run <- HIA_burden_total(short_trip_list, 
                calc_HIA_replicate,
                incidence_distrib_table, duration_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
                dis_vec, prop_relapse, dep_recovery, vsl, NULL, 1, FALSE) %>%
    mutate(run = i)
  
  MODAL_burden_total <- bind_rows(MODAL_burden_total, burden_run)
}



# Export HIA outcomes of 1000 replications
export(MODAL_burden_total, here("output", "RDS", "Modal shift", "HIA_modal_shift_1000replicate.rds"))




################################################################################################################################
#                         5. IC and MEDIAN - TOTAL BURDEN: PREVENTED CASES, DALY, MEDICAL, SOCIAL COSTS                        #
################################################################################################################################

##############################################################
#                       PER DISEASES                         #
##############################################################

# Import data
MODAL_burden_total <- import(here("output", "RDS", "Modal shift", "HIA_modal_shift_1000replicate.rds"))


# --------------------------------------
# MONTE-CARLO
# --------------------------------------
# IC95 and median 
  # Per disease
  set.seed(123)
  MODAL_burden_per_disease <- HIA_burden_IC(MODAL_burden_total, dis_vec, outcome_vec, calc_replicate_IC) 


  # Total for morbidity
  MODAL_burden_morbidity <- MODAL_burden_per_disease %>%
    filter(disease != c("mort", "dep")) %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Chronic diseases") %>%
    select(disease, everything()) 
  
  
  # Total for all diseases
  MODAL_burden_global <- MODAL_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  MODAL_burden <- bind_rows(MODAL_burden_per_disease, MODAL_burden_morbidity, MODAL_burden_global)
  


  
# --------------------------------------
# RUBIN'S RULE
# --------------------------------------
  # Per disease
  MODAL_Rubin_burden_per_disease <- HIA_burden_IC(MODAL_burden_total, dis_vec, outcome_vec, calc_IC_Rubin) 

  # Total for morbidity
  MODAL_Rubin_burden_morbidity <- MODAL_Rubin_burden_per_disease %>%
    filter(disease != c("mort", "dep")) %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Chronic diseases") %>%
    select(disease, everything()) 
  
  # Total for all diseases
  MODAL_Rubin_burden_global <- MODAL_Rubin_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  MODAL_Rubin_burden <- bind_rows(MODAL_Rubin_burden_per_disease, MODAL_Rubin_burden_morbidity, MODAL_Rubin_burden_global)




# Import 2019 data
MODAL_burden <- import(here("output", "Tables", "Modal shift", "HIA_modal_shift_1000replicate.xlsx"))
burden_2019 <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))


# Additional gains of the modal shift scenario
MODAL_burden_add <- MODAL_burden %>% 
    mutate(tot_cases = tot_cases - burden_2019[["tot_cases"]],
           tot_cases_low = tot_cases_low - burden_2019[["tot_cases_low"]], 
           tot_cases_up = tot_cases_up - burden_2019[["tot_cases_up"]],
           tot_daly = tot_daly - burden_2019[["tot_daly"]],
           tot_daly_low = tot_daly_low - burden_2019[["tot_daly_low"]],
           tot_daly_up = tot_daly_up - burden_2019[["tot_daly_up"]],
           tot_medic_costs = tot_medic_costs - burden_2019[["tot_medic_costs"]],
           tot_medic_costs_low = tot_medic_costs_low - burden_2019[["tot_medic_costs_low"]],
           tot_medic_costs_up = tot_medic_costs_up - burden_2019[["tot_medic_costs_up"]],
           tot_soc_costs = tot_soc_costs - burden_2019[["tot_soc_costs"]],
           tot_soc_costs_low = tot_soc_costs_low - burden_2019[["tot_soc_costs_low"]],
           tot_soc_costs_up = tot_soc_costs_up - burden_2019[["tot_soc_costs_up"]])  %>% 
    select(disease, tot_cases, tot_cases_low, tot_cases_up, tot_daly, tot_daly_low, tot_daly_up,
           tot_medic_costs, tot_medic_costs_low, tot_medic_costs_up, tot_soc_costs, tot_soc_costs_low, tot_soc_costs_up)



################################################################################################################################
#                                                       6. VISUALIZATION                                                       #
################################################################################################################################

# Import 2019 data
burden_2019 <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))



# Plot : Cases prevented had 50% of short car trips (<2km) been walked, compared to 2019 levels
plot_MODAL_cases_prev <- 
  ggplot() +

  # 2019 baseline
  geom_bar(data = burden_2019 %>%  
           filter(!disease %in% c("All", "Chronic diseases"))  %>% 
           mutate(disease = factor(disease, levels = c("mort", "cvd", "cancer", "diab2", "dem", "dep"))),
           mapping = aes(x = disease, y = tot_cases, fill = disease, alpha = "2019 baseline"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = burden_2019 %>%  
                filter(!disease %in% c("All", "Chronic diseases")),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "2019 baseline"),
                position = position_dodge(0.7),
                width = 0.25) +

  scale_fill_manual(values = colors_disease, labels = names_disease) +
  
  # Modal shift
  geom_bar(data = MODAL_burden %>%  
           filter(!disease %in% c("All", "Chronic diseases")) %>% 
           mutate(disease = factor(disease, levels = c("mort", "cvd", "cancer", "diab2", "dem", "dep"))), 
           mapping = aes(x = disease, y = tot_cases, fill = disease, alpha = "Modal shift scenario"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(values = c("2019 baseline" = 1, "Modal shift scenario" = 0.4), guide = "none") +
  
  geom_errorbar(data = MODAL_burden %>%  
                filter(!disease %in% c("All", "Chronic diseases")),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "Modal shift scenario"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
  scale_y_continuous(labels = label_comma(big.mark = ",", decimal.mark = ".")) +
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal() 

plot_MODAL_cases_prev 



################################################################################################################################
#                                            7. DISTANCE SHIFTED & CO2 EMISSIONS                                               #
################################################################################################################################
# Total km walked with IC and CO2 emissions prevented with IC per year
set.seed(123)

N <- 1000
tot_km_list <- vector("list", N)

for (i in 1:N) {
  print(i)

  tot_sample <- emp_short_driver %>%
    filter(!is.na(pond_jour) & !is.na(co2_depl)) %>%
    distinct(ident_ind, .keep_all = TRUE) %>%
    ungroup() %>%
    slice_sample(prop = perc, replace = TRUE) %>%
    as_survey_design(ids = ident_ind, weights = pond_jour) %>%
    summarise(
      tot_km             = survey_total(nbkm_car, na.rm = TRUE) * 365.25 / 7,
      tot_co2_shift      = survey_total(co2_adjusted, na.rm = TRUE) * 365.25 / 7,
      mean_co2_shift     = survey_mean(co2_adjusted, na.rm = TRUE),
      tot_co2_diminution = survey_total(co2_diminution, na.rm = TRUE) * 365.25 / 7,
      mean_co2_reduction = survey_mean(co2_prop_reduction, na.rm = TRUE)
    )

  tot_km_list[[i]] <- tot_sample
}

tot_km_drivers <- bind_rows(tot_km_list)



set.seed(123)
# Total km shifted
IC_Mkm <- calc_replicate_IC(tot_km_drivers, "tot_km") / 1e6                                       # in million km
tot_Mkm_IC <- data.frame(
  measure = "Total distance shifted (Mkm)",
  value = paste0(round(IC_Mkm["50%"], 3), " (", round(IC_Mkm["2.5%"], 3), " - ", round(IC_Mkm["97.5%"], 3), ")"))
    
IC_Mkm_Rubin <- calc_IC_Rubin (tot_km_drivers, "tot_km") / 1e6                                    # Rubin's rule
tot_Mkm_IC_Rubin <- data.frame(
  measure = "Total distance shifted (Mkm, Rubin)",
  value = paste0(round(IC_Mkm_Rubin[2], 3), " (", round(IC_Mkm_Rubin[1], 3), " - ", round(IC_Mkm_Rubin[3], 3), ")"))


# Total CO2 emissions prevented
IC_kt_co2_prev <- calc_replicate_IC(tot_km_drivers, "tot_co2_shift") *1e-9                                                             # CO2 emissions (in kt CO2)
tot_kt_co2_prev_IC <- data.frame(
  measure = "CO2 emissions prevented (kt CO2)",
  value = paste0(round(IC_kt_co2_prev["50%"], 3), " (", round(IC_kt_co2_prev["2.5%"], 3), " - ", round(IC_kt_co2_prev["97.5%"], 3), ")"))
    
IC_kt_co2_prev_Rubin <- calc_replicate_IC(tot_km_drivers, "tot_co2_shift") * 1e-9                                                      # Rubin's rule
tot_kt_co2_prev_IC_Rubin <- data.frame(
  measure = "CO2 emissions prevented (kt CO2, Rubin)",
  value = paste0(round(IC_kt_co2_prev_Rubin[2], 3), " (", round(IC_kt_co2_prev_Rubin[1], 3), " - ", round(IC_kt_co2_prev_Rubin[3], 3), ")"))


# Mean CO2 emissions prevented
mean_IC_kt_co2_prev <- calc_replicate_IC(tot_km_drivers, "mean_co2_shift") *1e-9                                                             # CO2 emissions (in kt CO2)
mean_kt_co2_prev_IC <- data.frame(
  measure = "Mean CO2 emissions prevented (kt CO2)",
  value = paste0(round(IC_kt_co2_prev["50%"], 3), " (", round(IC_kt_co2_prev["2.5%"], 3), " - ", round(IC_kt_co2_prev["97.5%"], 3), ")"))
    
mean_IC_kt_co2_prev_Rubin <- calc_replicate_IC(tot_km_drivers, "mean_co2_shift") * 1e-9                                                      # Rubin's rule
mean_kt_co2_prev_IC_Rubin <- data.frame(
  measure = "Mean CO2 emissions prevented (kt CO2, Rubin)",
  value = paste0(round(IC_kt_co2_prev_Rubin[2], 3), " (", round(IC_kt_co2_prev_Rubin[1], 3), " - ", round(IC_kt_co2_prev_Rubin[3], 3), ")"))


# Total diminution of CO2 emissions
IC_kt_co2_dim <- calc_replicate_IC(tot_km_drivers, "tot_co2_diminution") *1e-9                                                             # CO2 emissions (in kt CO2)
tot_kt_co2_dim_IC <- data.frame(
  measure = "Total CO2 emissions diminution (kt CO2)",
  value = paste0(round(IC_kt_co2_dim["50%"], 3), " (", round(IC_kt_co2_dim["2.5%"], 3), " - ", round(IC_kt_co2_dim["97.5%"], 3), ")"))
    
IC_kt_co2_dim_Rubin <- calc_replicate_IC(tot_km_drivers, "tot_co2_diminution") * 1e-9                                                      # Rubin's rule
tot_kt_co2_dim_IC_Rubin <- data.frame(
  measure = "Total CO2 emissions diminution (kt CO2, Rubin)",
  value = paste0(round(IC_kt_co2_dim_Rubin[2], 3), " (", round(IC_kt_co2_dim_Rubin[1], 3), " - ", round(IC_kt_co2_dim_Rubin[3], 3), ")"))


# Mean reduction of CO2 emissions
mean_IC_kt_co2_reduc <- calc_replicate_IC(tot_km_drivers, "mean_co2_reduction")                                                           # CO2 emissions (in kt CO2)
mean_kt_co2_reduc_IC <- data.frame(
  measure = "Mean CO2 emissions reduction (kt CO2)",
  value = paste0(round(mean_IC_kt_co2_reduc["50%"], 3), " (", round(mean_IC_kt_co2_reduc["2.5%"], 3), " - ", round(mean_IC_kt_co2_reduc["97.5%"], 3), ")"))
    
mean_IC_kt_co2_reduc_Rubin <- calc_replicate_IC(tot_km_drivers, "mean_co2_reduction")                                                # Rubin's rule
mean_kt_co2_reduc_IC_Rubin <- data.frame(
  measure = "Mean CO2 emissions reduction (kt CO2, Rubin)",
  value = paste0(round(mean_IC_kt_co2_reduc_Rubin[2], 3), " (", round(mean_IC_kt_co2_reduc_Rubin[1], 3), " - ", round(mean_IC_kt_co2_reduc_Rubin[3], 3), ")"))




tot_km_CO2 <- bind_rows(
  tot_Mkm_IC,
  tot_Mkm_IC_Rubin,
  tot_kt_co2_prev_IC,
  tot_kt_co2_prev_IC_Rubin,
  mean_kt_co2_prev_IC,
  mean_kt_co2_prev_IC_Rubin,
  tot_kt_co2_dim_IC,
  tot_kt_co2_dim_IC_Rubin,
  mean_kt_co2_reduc_IC,
  mean_kt_co2_reduc_IC_Rubin
)




################################################################################################################################
#                                                           8. DESCRIPTION                                                     #
################################################################################################################################


##############################################################
#                       DISTANCE DRIVEN                      #
##############################################################
# Total and mean distance driven of short car trips (< 2 km) per year (in km)
short_km_driven <- emp_car_trip  %>% 
  filter(!is.na(nbkm_car_jour) & nbkm_car_jour <= dist,
          !is.na(pond_jour))  %>%
  as_survey_design(ids = ident_ind, weights = pond_jour) %>% 
  summarise(tot_km = survey_total(nbkm_car, na.rm = T) * 365.25 / 7, 
            tot_mean = survey_mean(nbkm_car, na.rm = T))  



##############################################################
#                            DRIVERS                         #
##############################################################
# Calculate number of unique drivers (weighted) by summing the weight per unique `ident_ind`.
# Some individuals may have multiple trip rows; we take one row per `ident_ind` and use their `pond_indc`.
short_drivers <- emp_short_driver %>%
  distinct(ident_ind, .keep_all = TRUE)

nb_short_drivers <- tibble(
  nb_respondents = nrow(short_drivers),
  total = sum(short_drivers[["pond_indc"]], na.rm = TRUE))

nb_short_drivers




################################################################################################################################
#                                                      9. EXPORT DATA                                                          #
################################################################################################################################

## Plots
  ggsave(here("output", "Plots", "Modal shift", "modalshift_cases_prev.png"), plot = plot_MODAL_cases_prev)


## Tables
  # HIA of modal shift
    export(MODAL_burden, here("output", "Tables", "Modal shift", "HIA_modal_shift_1000replicate.xlsx"))
    export(MODAL_Rubin_burden, here("output", "Tables", "Modal shift", "HIA_modal_shift_Rubin_1000replicate.xlsx"))
    export(MODAL_burden_add, here("output", "Tables", "Modal shift", "HIA_modal_shift_added_1000replicate.xlsx"))

  # Total km walked with IC and CO2 emissions prevented with IC
    export(tot_km_CO2, here("output", "Tables", "Modal shift", "modalshift_km_CO2_emit.xlsx"))  
