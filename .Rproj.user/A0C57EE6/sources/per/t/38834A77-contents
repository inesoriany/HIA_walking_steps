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



# IC95 and median per scenario
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



##############################################################
#                       PER DISEASE                          #
##############################################################

for (dis in dis_vec) {
  dis_shift <- burden_tot %>% 
    filter(disease == dis) %>% 
    group_by(distance, percentage) %>% 
    summarise(tot_cases = sum(tot_cases),
              tot_cases_se = sum(tot_cases_se),
              tot_daly = sum(tot_daly),
              tot_daly_se = sum(tot_daly_se),
              tot_medic_costs = sum(tot_medic_costs) / 1e6,
              tot_medic_costs_se = sum(tot_medic_costs_se) / 1e6,
              tot_soc_costs = sum(tot_soc_costs * 1e-6),
              tot_soc_costs_se = sum(tot_soc_costs * 1e-6))
  assign(paste0(dis, "_shift"), dis_shift)
}


# IC95 and median
calc_replicate_IC(mort_shift[10,], "tot_cases")     # Scenario 50% <1km
calc_replicate_IC(mort_shift[10,], "tot_daly")      # Scenario 50% <1km



##############################################################
#                        MORBIDITY                           #
##############################################################
morbidity_shift <- burden_tot %>% 
  filter(disease != "mort") %>% 
  group_by(distance, percentage) %>% 
  summarise(tot_cases = sum(tot_cases),
            tot_cases_se = sum(tot_cases_se),
            tot_daly = sum(tot_daly),
            tot_daly_se = sum(tot_daly_se),
            tot_medic_costs = sum(tot_medic_costs) / 1e6,
            tot_medic_costs_se = sum(tot_medic_costs_se) / 1e6)

calc_replicate_IC(morbidity_shift[10,], "tot_cases")     # Scenario 50% <1km
calc_replicate_IC(morbidity_shift[10,], "tot_daly")      # Scenario 50% <1km




################################################################################################################################
#                                                     6. HIA - HEAT MAPS                                                       #
################################################################################################################################

##############################################################
#                        ALL DISEASES                        #
##############################################################
# All cases prevented per scenario
global_shift_cases <- ggplot(data = global_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_cases),
            color = "white") +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "Prevented morbi-mortality cases depending on different scenarios of car trips shifted to walk trips",
       fill = "Number of cases") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
global_shift_cases



# DALY prevented per scenario
global_shift_daly <- ggplot(data = global_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_daly),
            color = "white") +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "DALY prevented",
       fill = "Number of years") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
global_shift_daly



# Medical costs saved per scenario
global_shift_costs <- ggplot(data = global_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_medic_costs/1e3),                    # in billion €
            color = "white") +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "Direct medical (tangible) costs savings depending on different scenarios of car trips shifted to walk trips",
       fill = "Costs (in billion €)")+
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
global_shift_costs


# Intangible costs saved per scenario
global_shift_soc_costs <- ggplot(data = global_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_soc_costs/1e3),                    # in billion €
            color = "white") +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "Intangible costs savings depending on different scenarios of car trips shifted to walk trips",
       fill = "Costs (in billion €)")+
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
global_shift_soc_costs



##############################################################
#                        MORBIDITY                           #
##############################################################
# Chronic diseases prevented per scenario
morbidity_shift_cases <- ggplot(data = morbidity_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_cases)) +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "Diseases prevented",
       fill = "Number of cases") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
morbidity_shift_cases



##############################################################
#                        MORTALITY                           #
##############################################################
# Premature deaths prevented per scenario
mortality_shift_cases <- ggplot(data = mort_shift) +
  geom_tile(aes(x = distance, y = percentage*100, fill = tot_cases)) +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       title = "Premature deaths prevented",
       fill = "Number of cases") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
mortality_shift_cases



##############################################################
#                     MORBI-MORTALITY                        #
##############################################################
# List of the morbidity and mortality graphs
morbi_mortality_shift_graph <- list(morbidity_shift_cases, mortality_shift_cases)

# Theme for the common graph
common_theme <- theme(
  axis.title = element_text(size = 10, face = "bold"),
  strip.text = element_text(size = 8),
  axis.text.y = element_text(size = 9),
  legend.key.height = grid::unit(1, "cm"),
  legend.key.width = grid::unit(0.6, "cm")
)


# Apply the common theme to each graph
morbi_mortality_shift_graph <- lapply(morbi_mortality_shift_graph, 
                                      function(p) p + theme_minimal() + common_theme +
                                        scale_fill_gradient2(limits = c(0, max(morbidity_shift[["tot_cases"]])),          # Caliber a common scale
                                                             low = "#440154FF",
                                                             mid = "#1F968BFF",
                                                             high = "#FDE725FF",
                                                             midpoint = mean(c(0, max(morbidity_shift[["tot_cases"]]))) ) +
                                        theme(legend.position = "none")
)

# Make the legend appear once
morbi_mortality_shift_graph[[2]] <- morbi_mortality_shift_graph[[2]] + theme(legend.position = "right")

# Combine the graphs into one single figure
morbi_mortality_shift_cases <- wrap_plots(morbi_mortality_shift_graph, ncol = 2)

morbi_mortality_shift_cases




################################################################################################################################
#                                            7. DISTANCE SHIFTED & CO2 EMISSIONS                                               #
################################################################################################################################
# Total km walked per scenario with IC and CO2 emissions prevented per scenario with IC per year
set.seed(123)
N=100
tot_km_CO2 <- data.frame()
tot_km_CO2_scenario <- data.frame()

for (dist in dist_vec) {
  print(paste0("Distance = ", dist))
  drivers_dist <- emp_drive %>% 
    filter(!is.na(mdisttot_fin1) & mdisttot_fin1 <= dist)                             # Select the drivers under this distance
  
  for(perc in perc_vec) {
    print(paste0("Share = ", perc))
    
    tot_km_drivers <- data.frame()                                                    # 1 dataframe per scenario
    
    for(i in 1:N) {
      print(i)
      tot_sample <- drivers_dist %>%                                                  # Distances shifted for a year for 1 scenario
        filter(pond_jour != "NA") %>% 
        slice_sample(prop = perc) %>% 
        as_survey_design(ids= ident_ind, weights = pond_jour) %>% 
        summarise(tot_km = survey_total(mdisttot_fin1, na.rm = T)*365.25/7)
      
      tot_km_drivers <- bind_rows(tot_km_drivers, tot_sample)
    }
    
    IC_km <- calc_replicate_IC(tot_km_drivers, "tot_km") / 1e6                                       # in million km
    tot_km_IC <- paste0(round(IC_km["50%"], 3), " (", round(IC_km["2.5%"], 3), " - ", round(IC_km["97.5%"], 3), ")")
    
    IC_km_Rubin <- calc_IC_Rubin (tot_km_drivers, "tot_km") / 1e6                                    # Rubin's rule
    tot_km_IC_Rubin <- paste0(round(IC_km_Rubin[2], 3), " (", round(IC_km_Rubin[1], 3), " - ", round(IC_km_Rubin[3], 3), ")")
    
    
    IC_kt <- IC_km * CO2_emit *1e-3                                                                  # CO2 emissions (in kt CO2)
    tot_kt_IC <- paste0(round(IC_kt["50%"], 3), " (", round(IC_kt["2.5%"], 3), " - ", round(IC_kt["97.5%"], 3), ")")
    
    IC_kt_Rubin <- IC_km_Rubin* CO2_emit * 1e-3                                                      # Rubin's rule
    tot_kt_IC_Rubin <- paste0(round(IC_kt_Rubin[2], 3), " (", round(IC_kt_Rubin[1], 3), " - ", round(IC_kt_Rubin[3], 3), ")")
    
    
    tot_km_CO2_scenario <- bind_rows(tot_km_CO2_scenario, data.frame(
      distance = dist,
      percentage = perc,
      km = IC_km[["50%"]],
      km_low = IC_km[["2.5%"]],
      km_sup = IC_km[["97.5%"]], 
      Rubin_km = IC_km_Rubin[2],
      Rubin_km_low = IC_km_Rubin[1],
      Rubin_km_sup = IC_km_Rubin[3],
      CO2_kt = IC_kt[["50%"]],
      CO2_kt_low = IC_kt[["2.5%"]],
      CO2_kt_sup = IC_kt[["97.5%"]],
      Rubin_kt = IC_kt_Rubin[2],
      Rubin_kt_low = IC_kt_Rubin[1],
      Rubin_kt_sup = IC_kt_Rubin[3] ))
    
    tot_km_CO2 <- bind_rows(tot_km_CO2, data.frame(
      distance = dist,
      percentage = paste0(perc*100, "%"),
      total_millions_km = tot_km_IC,
      Rubin_total_millions_km = tot_km_IC_Rubin,
      CO2_emissions_kt = tot_kt_IC,
      Rubin_CO2_emissions_kt = tot_kt_IC_Rubin))
    
  }
}


# Export replications - Total km walked per scenario
export(tot_km_CO2_scenario, here("output", "RDS", "Modal shift", "modalshift_tot_km_CO2.rds"))



##############################################################
#                   HEATMAP - CO2 EMISSIONS                  #
##############################################################

# Import replications - Total km walked per scenario with IC and CO2 emissions prevented per scenario with IC
tot_km_CO2_scenario <- import(here("output", "RDS", "Modal shift", "modalshift_tot_km_CO2.rds"))


# Heatmap (in kilotons)
CO2_shift <- ggplot(data = tot_km_CO2_scenario) +
  geom_tile(aes(x = distance, y = percentage*100, fill = CO2_kt),
            color = "white") +
  scale_fill_viridis() +
  labs(x = "Distances of car trips shifted (km)",
       y = "Share shifted (%)", 
       fill = "CO2 emissions (kto)") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks = element_line(linewidth = 0.4),
        axis.title = element_text(size = 12, face = "bold"),
        
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  theme_minimal()
CO2_shift




################################################################################################################################
#                                                           8. DESCRIPTION                                                     #
################################################################################################################################
##############################################################
#                       DISTANCE DRIVEN                      #
##############################################################
# Total distance driven of all drivers for each scenario distance per year (in km)
km_driven <- data.frame()
for(dist in dist_vec) {
  drivers_dist <- emp_drive %>% 
    filter(!is.na(mdisttot_fin1) & mdisttot_fin1 <= dist)
  km_driven_dist <- drivers_dist %>% 
    filter(pond_jour != "NA") %>% 
    as_survey_design(ids = ident_ind, weights = pond_jour) %>% 
    summarise(tot_km = survey_total(mdisttot_fin1, na.rm = T) * 365.25 / 7,
              tot_mean = survey_mean(mdisttot_fin1, na.rm = T)) %>% 
    mutate(distance = dist)
  
  km_driven <- bind_rows(km_driven, km_driven_dist)
}



##############################################################
#                      NUMBER OF DRIVERS                     #
##############################################################
nb_short_drivers <-  data.frame()

for (dist in dist_vec) {
  drivers_dist <- emp_drive %>% 
    filter(!is.na(mdisttot_fin1) & mdisttot_fin1 <= dist)                             # Select the drivers under this distance
  
  tot_short_drivers <- sum(drivers_dist$pond_indc)
  
  for (perc in perc_vec) {
    sample_drivers <- drivers_dist %>% 
      slice_sample(prop = perc)
    
    sample_short_drivers <- sum(sample_drivers$pond_indc)
    
    prop_short_drivers <- sample_short_drivers / tot_short_drivers
    
    nb_short_drivers <- bind_rows(nb_short_drivers, 
                                  data.frame(
                                    distance = dist,
                                    percentage = perc,
                                    nb_shifted_drivers = sample_short_drivers,
                                    prop_short_drivers = prop_short_drivers))
  }
}




################################################################################################################################
#                                                      9. EXPORT DATA                                                          #
################################################################################################################################

## Tables
  # HIA of modal shifts
    # Total all diseases
    export(global_shift_IC, here("output", "Tables", "Modal shift", "HIA_modal_shift_100replicate.xlsx"))
    # Mortality
    export(mort_shift, here("output", "Tables", "Modal shift", "Mortality_modal_shift_100replicate.xlsx"))
    # Morbidity
    export(morbidity_shift, here("output", "Tables", "Modal shift", "Morbidity_modal_shift_100replicate.xlsx"))
  
  # Total km walked per scenario with IC and CO2 emissions prevented per scenario with IC
    export(tot_km_CO2, here("output", "Tables", "Modal shift", "modalshift_tot_km_CO2_emit.xlsx"))  

  

## Plots
  # Total all diseases
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_cases_prevented.png"), plot = global_shift_cases)
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_daly_prevented.png"), plot = global_shift_daly)
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_costs_saved.png"),plot = global_shift_costs)
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_soc_costs_saved.png"),plot = global_shift_soc_costs)
  # Morbidity
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_morbidity_prevented.png"), plot = morbidity_shift_cases)
  # Mortality
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_mortality_prevented.png"), plot = mortality_shift_cases)
  # Morbi-mortality (combined)
  ggsave(here("output", "Plots", "Modal shift", "plot_modalshift_morbi_mortality_prevented.png"), plot = morbi_mortality_shift_cases)

  # CO2 emissions (in kilotons)
  ggsave(here("output", "Plots", "Modal shift", "modalshift_CO2_emit.png"), plot = CO2_shift)




