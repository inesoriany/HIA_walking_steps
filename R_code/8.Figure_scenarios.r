###############################################
##########  SIMULATION SCENARIOS  #############
###############################################




################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  ggplot2      # Data visualization
)



################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# 2019 baseline 
burden_2019 <- import(here("output", "Tables", "2019", "HIA_per_disease.xlsx"))

# Best case scenario - 7000 steps
BEST_burden <- import(here("output", "Tables", "7000 steps", "HIA_7000steps_1000replicate.xlsx"))

# Modal shift scenario
MODAL_burden <- import(here("output", "Tables", "Modal shift", "HIA_modal_shift_1000replicate.xlsx"))

# Best practice scenario
PRACT_burden <- import(here("output", "Tables", "Best practice", "Residence", "HIA_best_practice_1000replicate_RES.xlsx"))


################################################################################################################################
#                                                   4. DATA PREPARATION                                                        #
################################################################################################################################

# Recap burden of all scenarios
burden_scenarios <- 
    bind_rows(burden_2019 %>% filter(disease == "All") %>% mutate(scenario = "2019 baseline"),
              BEST_burden %>% filter(disease == "All") %>% mutate(scenario = "Best case scenario"),
              MODAL_burden %>% filter(disease == "All") %>% mutate(scenario = "Modal shift scenario"),
              PRACT_burden %>% filter(disease == "All", area_type == "All") %>% mutate(scenario = "Best practice scenario")) %>%
  mutate(scenario = factor(scenario, levels = c("2019 baseline", "Best case scenario", "Modal shift scenario", "Best practice scenario")))


# Best practice scenario
PRACT_area_burden <- PRACT_burden  %>% 
    filter(disease == "All", area_type != "All")  %>% 
    mutate(scenario = factor("Best practice scenario")) %>% 
    mutate(area_type = factor(area_type, levels = c("urban", "periurban", "rural"))) %>%
    arrange(desc(area_type)) %>%
    mutate(tot_cases_cum = cumsum(tot_cases),
           tot_cases_low_cum = cumsum(tot_cases_low),
           tot_cases_up_cum = cumsum(tot_cases_up))


################################################################################################################################
#                                                    5. VISUALIZATION                                                          #
################################################################################################################################

# DALY comparison (with area type stratification for the best practice scenario)
plot_scenarios_daly <- 
  ggplot() +

  geom_bar(data = burden_scenarios,
           mapping = aes(x = scenario, y = tot_cases),
           fill = "#999999", width = 0.7, stat = "identity") +

  geom_errorbar(data = burden_scenarios,
                mapping = aes(x = scenario, ymin = tot_cases_low, ymax = tot_cases_up),
                width = 0.25) +

  ylab("DALY prevented") +
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

plot_scenarios_daly




# DALY comparison (with area type stratification for the best practice scenario)
plot_scenarios_daly_area <- 
  ggplot() +

  geom_bar(data = burden_scenarios,
           mapping = aes(x = scenario, y = tot_cases),
           fill = "#999999", width = 0.7, stat = "identity") +

  geom_errorbar(data = burden_scenarios,
                mapping = aes(x = scenario, ymin = tot_cases_low, ymax = tot_cases_up),
                width = 0.25) +
  

  # Best practice scenario - stratification
  geom_bar(data = PRACT_area_burden,
           mapping = aes(x = scenario, y = tot_cases, fill = area_type),
           width = 0.7, stat = "identity", position = "stack") +
    scale_fill_manual(values = c(urban = "#1b9e77",
                                 periurban = "#d95f02",
                                 rural = "#7570b3"),
                      name = "Area type") +

  geom_errorbar(data = PRACT_area_burden,
                mapping = aes(x = scenario, ymin = tot_cases_low_cum, ymax = tot_cases_up_cum),
                width = 0.25) +
  

  ylab("DALY prevented") +
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

plot_scenarios_daly_area



################################################################################################################################
#                                                  6. INCREMENTAL BENEFIT                                                      #
################################################################################################################################

# best case : pas incremental
MODAL_burden_add <- import(here("output", "Tables", "Modal shift", "HIA_modal_shift_added_1000replicate.xlsx"))



################################################################################################################################
#                                                      7. EXPORT DATA                                                          #
################################################################################################################################

# Recap table
  export(burden_scenarios, here("output", "Tables", "Scenarios", "HIA_scenarios.xlsx"))


# Plot
  ggsave(here("output", "Plots", "Scenarios", "scenarios_DALY_prev.png"), plot = plot_scenarios_daly)
  ggsave(here("output", "Plots", "Scenarios", "scenarios_DALY_area_prev.png"), plot = plot_scenarios_daly_area)


