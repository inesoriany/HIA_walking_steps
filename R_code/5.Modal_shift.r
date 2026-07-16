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

# Import drivers dataset
emp_car_trip <- import(here("data_clean", "EMP_dis_car_trips.xlsx"))


# Incidence distribution table
incidence_distrib_table <- import(here("data_clean", "Diseases", "incidence_distrib_table.xlsx"))

# Depression duration distribution table
dep_distrib_table <- import(here("data_clean", "Diseases", "dep_duration_distrib_table.xlsx"))

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

  # Percentage of car trips shifted to walking
  perc <- 0.5


################################################################################################################################
#                                                    3. DATA PREPARATION                                                       #
################################################################################################################################

# Initialization
emp_short_trip <- emp_car_trip %>% 
    filter(!is.na(nbkm_car) & nbkm_car <= dist)  %>% 
    # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000,round(step_commute / 100) * 100 + baseline_step))



################################################################################################################################
#                                              4. HEALTH IMPACT ASSESSMENT                                                     #
################################################################################################################################
set.seed(123)

N <- 1000

MODAL_burden_total <- data.frame()

for (i in 1:N) {
    print(paste0("Run ", i))

    short_trip_list <- list()

    for (dis in dis_vec) {
      emp_short_sample <- emp_short_trip %>%
      filter(disease == dis)  %>% 
      slice_sample(prop = perc)  %>% 
      mutate (disease = dis)
        
      short_trip_list[[dis]] <- emp_short_sample
    }

    burden_run <- HIA_burden_total(short_trip_list, 
             incidence_distrib_table, dep_distrib_table, reduction_risk_distrib_table, dw_distrib_table,
             dis_vec, prop_relapse, duration_recovery, vsl, NULL, 1, FALSE)  %>% 
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
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
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
    filter(disease != "mort") %>% 
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "Morbidity") %>%
    select(disease, everything()) 
  
  # Total for all diseases
  MODAL_Rubin_burden_global <- MODAL_Rubin_burden_per_disease %>%
    summarise(across(where(is.numeric), 
                     ~ sum(.x, na.rm = TRUE) )) %>%
    mutate(disease = "All") %>%
    select(disease, everything()) 
  
  # Gather results
  MODAL_Rubin_burden <- bind_rows(MODAL_Rubin_burden_per_disease, MODAL_Rubin_burden_morbidity, MODAL_Rubin_burden_global)




################################################################################################################################
#                                                       6. VISUALIZATION                                                       #
################################################################################################################################

# Import 2019 data
burden_2019 <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))


# Data preparation
MODAL_burden_add <- MODAL_burden  %>% 
    mutate(tot_cases = tot_cases + burden_2019[["tot_cases"]],
           tot_cases_low = tot_cases_low + burden_2019[["tot_cases_low"]], 
           tot_cases_up = tot_cases_up + burden_2019[["tot_cases_up"]],
           tot_daly = tot_daly + burden_2019[["tot_daly"]],
           tot_daly_low = tot_daly_low + burden_2019[["tot_daly_low"]],
           tot_daly_up = tot_daly_up + burden_2019[["tot_daly_up"]])  %>% 
    select(disease, tot_cases, tot_cases_low, tot_cases_up, tot_daly, tot_daly_low, tot_daly_up)


# Plot : Cases prevented had 50% of short car trips (<2km) been walked, compared to 2019 levels
plot_MODAL_cases_prev <- 
  ggplot() +

  # 2019 baseline
  geom_bar(data = burden_2019 %>%  
           filter(!disease %in% c("All", "Morbidity"))  %>% 
           mutate(disease = factor(disease, levels = c("mort", "cancer", "cvd", "diab2", "dem", "dep"))),
           mapping = aes(x = disease, y = tot_cases, fill = disease, alpha = "2019 baseline"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = burden_2019 %>%  
                filter(!disease %in% c("All", "Morbidity")),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "2019 baseline"),
                position = position_dodge(0.7),
                width = 0.25) +

  scale_fill_manual(values = colors_disease, labels = names_disease) +
  
  # Modal shift
  geom_bar(data = MODAL_burden_add %>%  
           filter(!disease %in% c("All", "Morbidity")), 
           mapping = aes(x = disease, y = tot_cases, fill = disease, alpha = "Modal shift scenario"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "Modal shift scenario" = 0.4)) +
  
  geom_errorbar(data = MODAL_burden_add %>%  
                filter(!disease %in% c("All", "Morbidity")),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, alpha = "Modal shift scenario"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
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

tot_km_CO2 <- data.frame()

emp_short_sample <- emp_car_trip %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist)                                       # Select the drivers under this distance

tot_km_drivers <- data.frame()                                                    # 1 dataframe per scenario
for(i in 1:N) {
  print(i)
  tot_sample <- emp_short_sample %>%                                                  # Distances shifted for a year for 1 scenario
  filter(pond_jour != "NA") %>% 
  slice_sample(prop = perc) %>% 
  as_survey_design(ids= ident_ind, weights = pond_jour) %>% 
  summarise(tot_km = survey_total(nbkm_car, na.rm = T)*365.25/7)
      
  tot_km_drivers <- bind_rows(tot_km_drivers, tot_sample)
}


set.seed(123)
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
#                            DRIVERS                         #
##############################################################
# Number of drivers with short car trips (< 2 km) per year
short_drivers <- emp_car_trip  %>% 
  filter(!is.na(nbkm_car) & nbkm_car <= dist)

sample_drivers <- short_drivers %>% 
  slice_sample(prop = perc) 

# Calculate number of unique drivers (weighted) by summing the weight per unique `ident_ind`.
# Some individuals may have multiple trip rows; we take one row per `ident_ind` and use their `pond_indc`.
nb_short_drivers <- short_drivers %>%
  distinct(ident_ind, .keep_all = TRUE) %>%
  summarise(total = sum(pond_indc, na.rm = TRUE)) %>%
  pull(total)

nb_short_drivers

# For the sampled (shifted) drivers, also count unique individuals and sum their weights
nb_shifted_drivers <- sample_drivers %>%
  distinct(ident_ind, .keep_all = TRUE) %>%
  summarise(total = sum(pond_indc, na.rm = TRUE)) %>%
  pull(total)

nb_shifted_drivers




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
    export(tot_km_CO2, here("output", "Tables", "Modal shift", "modalshift_tot_km_CO2_emit.xlsx"))  
