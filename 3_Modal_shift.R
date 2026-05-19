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
emp_car_trip <- import(here("data_clean", "EMP_dis_car_trips.xlsx"))

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
dis_vec <- c("mort", "cvd", "cancer", "diab2", "dem", "dep")

# HIA outcomes
outcome_vec <- c("tot_cases", "tot_daly", "tot_medic_costs", "tot_soc_costs")

# Modal shift scenario
  # Driven distance shifted to walked (km)
  dist <- 2 

  # Percentage of car trips shifted to walking
  perc <- 0.5




################################################################################################################################
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################
# Initialization
emp_car_trip <- emp_car_trip %>% 
  # Day step and baseline at 2000
  mutate(step = pmin(12000,round(step_commute / 100) * 100 + baseline_step))



################################################################################################################################
#                                              5. HEALTH IMPACT ASSESSMENT                                                     #
################################################################################################################################
set.seed(123)
burden_tot <- data.frame()                                         # Burden for each run and each disease

drivers_dist <- emp_car_trip  %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist)                      # Select the drivers under this distance

for(dis in dis_vec){
  print(paste0("Disease : ", dis))
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
          mutate(run = i,
                disease = dis)
          
    ## ----------------------------------------
    ## 2. Draw RR
    ## ----------------------------------------
    sim_id <- rr_distrib_table %>%
      filter(disease == dis) %>%
      pull(simulation_id) %>%
      unique() %>%
      sample(1)
          
    rr_baseline <- rr_distrib_table %>%
      filter(disease == dis,
             step == baseline_step,
             simulation_id == sim_id) %>%
      pull(rr_interpolated)
          
    rr_step <- rr_distrib_table %>%
      filter(disease == dis,
             simulation_id == sim_id) %>%
      select(step, rr_interpolated) %>%
      mutate(reduction_risk = (rr_baseline - rr_interpolated) / rr_baseline)
          
    burden_sample <- burden_sample %>%
      left_join(rr_step %>% select(step, reduction_risk),
                by = "step") %>%
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
    burden_sample <- daly(burden_sample, dis)
    burden_sample <- medic_costs(burden_sample, dis)
          
    burden_sample <- burden_sample %>%
      mutate(soc_costs = daly * vsl)
          
    ## ----------------------------------------
    ## 5. Aggregate burden by run
    ## ----------------------------------------
    surv_run <- burden_sample %>%
      as_survey_design(ids = ident_ind,
                       weights = pond_indc)
          
    burden_run[[i]] <- surv_run %>%
    summarise(
      tot_cases = survey_total(cases, na.rm = TRUE),
      tot_daly = survey_total(daly, na.rm = TRUE),
      tot_medic_costs = survey_total(medic_costs, na.rm = TRUE),
      tot_soc_costs = survey_total(soc_costs, na.rm = TRUE)) %>%
    mutate(run = i,
           disease = dis,
           distance = dist,
           percentage = perc)
    }
        
  ## ----------------------------------------
  ## 6. Combine runs
  ## ----------------------------------------
  burden_tot <- bind_rows(burden_tot_run, bind_rows(burden_run))
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
  summarise(tot_cases = sum(tot_cases),
            tot_cases_se = sum(tot_cases_se),
            tot_daly = sum(tot_daly),
            tot_daly_se = sum(tot_daly_se),
            tot_medic_costs = sum(tot_medic_costs) * 1e-6,                                # in million €
            tot_medic_costs_se = sum(tot_medic_costs_se) * 1e-6,
            tot_soc_costs = sum(tot_soc_costs * 1e-6),
            tot_soc_costs_se = sum(tot_soc_costs_se * 1e-6))



# IC95 and median per scenario
set.seed(123)
global_shift_IC <- data.frame()
global_cases_IC <- as.data.frame(t(calc_replicate_IC(global_shift, "tot_cases"))) %>% 
  rename(cases_low = "2.5%", cases_mid = "50%", cases_sup = "97.5%") %>% 
  mutate(distance = dist, percentage = perc)
    
global_daly_IC <- as.data.frame(t(calc_replicate_IC(global_shift, "tot_daly"))) %>% 
  rename(daly_low = "2.5%", daly_mid = "50%", daly_sup = "97.5%") %>% 
  mutate(distance = dist, percentage = perc)
    
global_medic_costs_IC <- as.data.frame(t(calc_replicate_IC(global_shift, "tot_medic_costs"))) %>% 
  rename(medic_costs_low = "2.5%", medic_costs_mid = "50%", medic_costs_sup = "97.5%") %>% 
  mutate(distance = dist, percentage = perc)
    
global_soc_costs_IC <- as.data.frame(t(calc_replicate_IC(global_shift, "tot_soc_costs"))) %>% 
  rename(soc_costs_low = "2.5%", soc_costs_mid = "50%", soc_costs_sup = "97.5%") %>% 
  mutate(distance = dist, percentage = perc)
    
global_shift_IC <- left_join(global_cases_IC, global_daly_IC, by = c("distance", "percentage"))
global_shift_IC <- left_join(global_shift_IC, global_medic_costs_IC, by = c("distance", "percentage"))
global_shift_IC <- left_join(global_shift_IC, global_soc_costs_IC, by = c("distance", "percentage")) %>% 
  relocate(distance, percentage)



##############################################################
#                       PER DISEASE                          #
##############################################################
disease_shift <- burden_tot %>% 
  group_by(disease) %>%
  summarise(tot_cases = sum(tot_cases),
            tot_cases_se = sum(tot_cases_se),
            tot_daly = sum(tot_daly),
            tot_daly_se = sum(tot_daly_se),
            tot_medic_costs = sum(tot_medic_costs) * 1e-6,                                # in million €
            tot_medic_costs_se = sum(tot_medic_costs_se) * 1e-6,
            tot_soc_costs = sum(tot_soc_costs * 1e-6),
            tot_soc_costs_se = sum(tot_soc_costs_se * 1e-6))

# IC95 and median per scenario
set.seed(123)
disease_shift_IC <- HIA_burden_IC(disease_shift, dis_vec, outcome_vec, calc_replicate_IC) %>%
  mutate(distance = dist, percentage = perc)  %>% 
  relocate(distance, percentage)


# Re-organize by decreasing order
disease_shift_IC_order <- disease_shift_IC %>% 
    left_join(disease_shift %>% select(disease, TOTAL_mixed = tot_cases), by = "disease") %>% 
    arrange(desc(tot_cases)) %>%                      
    mutate(disease = factor(disease, levels = unique(disease)))  %>% 
    select(-c(TOTAL_mixed))


##############################################################
#                        MORBIDITY                           #
##############################################################
morbidity_shift <- burden_tot %>% 
  filter(disease != "mort") %>% 
  summarise(tot_cases = sum(tot_cases),
            tot_cases_se = sum(tot_cases_se),
            tot_daly = sum(tot_daly),
            tot_daly_se = sum(tot_daly_se),
            tot_medic_costs = sum(tot_medic_costs) / 1e6,
            tot_medic_costs_se = sum(tot_medic_costs_se) / 1e6)

calc_replicate_IC(morbidity_shift, "tot_cases")  
calc_replicate_IC(morbidity_shift, "tot_daly") 




################################################################################################################################
#                                                       6. VISUALIZATION                                                       #
################################################################################################################################

# Import 2019 data
disease_2019 <- import(here("output", "Tables", "2019", "cases_prev_2019.xlsx"))

# Plot : Cases prevented had 50% of short car trips (<2km) been walked, compared to 2019 levels
plot_shift_cases_prev <- 
  ggplot() +

  # 2019 baseline
  geom_bar(data = disease_2019,
           mapping = aes(x = disease, y = tot_cases_mid, fill = disease, alpha = "2019 baseline"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = disease_2019,
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "2019 baseline"),
                position = position_dodge(0.7),
                width = 0.25) +

  scale_fill_manual(values = colors_disease, labels = names_disease) +
  
  # Modal shift
  geom_bar(data = disease_shift_IC_order, 
           mapping = aes(x = disease, y = tot_cases, fill = disease, alpha = "Modal shift scenario"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "Modal shift scenario" = 0.4)) +
  
  geom_errorbar(data = disease_shift_IC_order,
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "Modal shift scenario"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal() 

plot_shift_cases_prev 


################################################################################################################################
#                                            7. DISTANCE SHIFTED & CO2 EMISSIONS                                               #
################################################################################################################################
# Total km walked with IC and CO2 emissions prevented with IC per year
set.seed(123)
N=100
tot_km_CO2 <- data.frame()

drivers_dist <- emp_car_trip %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist)                                       # Select the drivers under this distance

tot_km_drivers <- data.frame()                                                    # 1 dataframe per scenario
for(i in 1:N) {
  print(i)
  tot_sample <- drivers_dist %>%                                                  # Distances shifted for a year for 1 scenario
  filter(pond_jour != "NA") %>% 
  slice_sample(prop = perc) %>% 
  as_survey_design(ids= ident_ind, weights = pond_jour) %>% 
  summarise(tot_km = survey_total(nbkm_car, na.rm = T)*365.25/7)
      
  tot_km_drivers <- bind_rows(tot_km_drivers, tot_sample)
}
    
IC_Mkm <- calc_replicate_IC(tot_km_drivers, "tot_km") / 1e6                                       # in million km
tot_Mkm_IC <- paste0(round(IC_Mkm["50%"], 3), " (", round(IC_Mkm["2.5%"], 3), " - ", round(IC_Mkm["97.5%"], 3), ")")
    
IC_Mkm_Rubin <- calc_IC_Rubin (tot_km_drivers, "tot_km") / 1e6                                    # Rubin's rule
tot_Mkm_IC_Rubin <- paste0(round(IC_Mkm_Rubin[2], 3), " (", round(IC_Mkm_Rubin[1], 3), " - ", round(IC_Mkm_Rubin[3], 3), ")")
    
    
IC_kt <- IC_Mkm * 1e6 * CO2_emit *1e-9                                                             # CO2 emissions (in kt CO2)
tot_kt_IC <- paste0(round(IC_kt["50%"], 3), " (", round(IC_kt["2.5%"], 3), " - ", round(IC_kt["97.5%"], 3), ")")
    
IC_kt_Rubin <- IC_Mkm_Rubin * 1e6 * CO2_emit * 1e-9                                                # Rubin's rule
tot_kt_IC_Rubin <- paste0(round(IC_kt_Rubin[2], 3), " (", round(IC_kt_Rubin[1], 3), " - ", round(IC_kt_Rubin[3], 3), ")")
    
tot_km_CO2 <- bind_rows(tot_km_CO2, data.frame(
  distance = dist,
  percentage = paste0(perc*100, "%"),
  total_millions_km = tot_Mkm_IC,
  Rubin_total_millions_km = tot_Mkm_IC_Rubin,
  CO2_emissions_kt = tot_kt_IC,
  Rubin_CO2_emissions_kt = tot_kt_IC_Rubin))




################################################################################################################################
#                                                           8. DESCRIPTION                                                     #
################################################################################################################################
# Tout re-simplifier

##############################################################
#                       DISTANCE DRIVEN                      #
##############################################################
# Total and mean distance driven of short car trips (< 2 km) per year (in km)
short_km_driven <- emp_car_trip  %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist,
          !is.na(pond_jour))  %>%
  as_survey_design(ids = ident_ind, weights = pond_jour) %>% 
  summarise(tot_km = survey_total(nbkm_car, na.rm = T) * 365.25 / 7, 
            tot_mean = survey_mean(nbkm_car, na.rm = T))  



##############################################################
#                      NUMBER OF DRIVERS                     #
##############################################################
short_drivers <- emp_car_trip  %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist)

sample_drivers <- short_drivers %>% 
  slice_sample(prop = perc) 

nb_short_drivers <-sum(short_drivers$pond_indc)
nb_short_drivers 
nb_shifted_drivers <- sum(sample_drivers$pond_indc)
nb_shifted_drivers 


################################################################################################################################
#                                                      9. EXPORT DATA                                                          #
################################################################################################################################

## Tables
  # HIA of modal shifts
    # Total all diseases
    export(global_shift_IC, here("output", "Tables", "Modal shift", "HIA_modal_shift_100replicate.xlsx"))
    # Disease
    export(disease_shift_IC, here("output", "Tables", "Modal shift", "disease_modal_shift_100replicate.xlsx"))
  
  # Total km walked per scenario with IC and CO2 emissions prevented per scenario with IC
    export(tot_km_CO2, here("output", "Tables", "Modal shift", "modalshift_tot_km_CO2_emit.xlsx"))  

  

## Plots
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_cases_prev.png"), plot = plot_shift_cases_prev)