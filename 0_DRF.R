#################################################
################ RR DISTRIBUTIONS ###############
#################################################

# Files needed:
  # rr_table_quanti.xlsx


# Files outputted:
  # rr_central_interpolated.rds
  # rr_sim_interpolated.rds : 1000 simulations every 10 interpolated points
  # rr_interpolated_mean.csv : mean, IC95%
  # DRF plots




################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################

pacman::p_load(
  rio,                 # file import/export
  here,                # file path management
  dplyr,               # data manipulation
  tidyr,               # data manipulation
  tidyverse,           # Data management, ggplot included
  patchwork,           # graphs combination
  splines              # interpolation
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# RR table for each outcome

  # Central values
  rr_table_mid <- import(here("data", "rr_table_quanti.xlsx"), sheet = "Central")
  
  # IC95 lower
  rr_table_low <- import(here("data", "rr_table_quanti.xlsx"), sheet = "Lower")
  
  # IC95 upper
  rr_table_up <- import(here("data", "rr_table_quanti.xlsx"), sheet = "Upper") 


# Import functions
source(here("0_Functions.R"))  


################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Number of simulations
n <- 1000
  
# Diseases considered
dis_vec = c("mort", "cvd", "cancer", "diab2", "dem", "dep")


################################################################################################################################
#                                                     4. PREPARE DATA                                                          #
################################################################################################################################
  
# One unique table 
rr_table <- rr_table_mid %>% 
  left_join(rr_table_low, by = "step") %>%
  left_join(rr_table_up, by = "step") %>% 
  mutate(across(where(is.character), as.numeric))

  
for (dis in dis_vec) {
  rr_table <- rr_table %>% 
    rename(
      !!paste0(dis, "_mid") := !!sym(paste0(dis, ".x")),
      !!paste0(dis, "_low") := !!sym(paste0(dis, ".y")),
      !!paste0(dis, "_up")  := !!sym(paste0(dis))
    )
}
  


# Pivot longer
rr_table_long <- rr_table %>%
  pivot_longer(
    cols = matches(".*_(mid|low|up)$"),
    names_to = c("disease", ".value"),
    names_sep = "_"
  )



################################################################################################################################
#                                             5. GENERATE RR NORMAL DISTRIBUTIONS                                              #
################################################################################################################################

# Generate RR normal distributions
set.seed(123)
rr_table_long <- rr_table_long %>%
  rowwise() %>%
  mutate(rr_distrib = list(
    generate_RR_distrib(mid, low, up, n)
  ))



################################################################################################################################
#                                         6. SORT THE RR DISTRIBUTIONS IN ASCENDING ORDER                                      #
################################################################################################################################

rr_table_long <- rr_table_long %>% 
  rowwise() %>% 
  mutate(rr_distrib = list(sort(unlist(rr_distrib)))) %>%
  ungroup()


################################################################################################################################
#                                                       7. INTERPOLATION                                                       #
################################################################################################################################


##############################################
#       RR Central values interpolation      #
##############################################

  # Interpolation
rr_central_list <- list()

set.seed(456)
for (dis in dis_vec) {
  for (metric in c("mid", "low", "up")) {
    rr_central_list[[paste0(dis, "_", metric)]] <- interpolate_rr(rr_table, dis, metric)
  }
}

rr_central_full <- bind_rows(rr_central_list)

rr_central_table <- rr_central_full %>%
  pivot_wider(
    id_cols = c(step, disease),
    names_from = metric,          
    values_from = rr_interpolated
  ) %>%
  arrange(disease, step) 




##############################################
#    RR normal distributions interpolation   #
##############################################


rr_table_long <- rr_table_long %>% 
  unnest_wider(rr_distrib, names_sep = "_") %>%     # separate the rr_distrib column into multiple columns
  pivot_longer(
    cols = starts_with("rr_distrib_"), 
    names_to = "simulation_id",                     # column name for the simulation ID
    values_to = "simulated_rr"                      # column name for the simulated RR values
  ) %>% 
  mutate(simulation_id = as.numeric(str_remove(simulation_id, "rr_distrib_")))      # simulation ID as a numeric value



# Interpolation
rr_table_interpolated <- rr_table_long %>% 
  group_by(disease, simulation_id) %>% 
  complete(step = seq(0, 12000, by = 100)) %>% 
  arrange(step) %>% 
  mutate(rr_interpolated = case_when(
    # quadratic (degree = 2)
    disease %in% c("mort", "cvd")   ~ if_else(is.na(simulated_rr),spline(step, simulated_rr, xout = step, method = "fmm")$y, simulated_rr),
    # cubic (degree = 3)
    disease %in% c("dem")           ~ if_else(is.na(simulated_rr), spline(step, simulated_rr, xout = step, method = "natural")$y, simulated_rr),
    # linear
    TRUE                            ~ if_else(is.na(simulated_rr), approx(step, simulated_rr, xout = step, method = "linear", rule = 1)$y, simulated_rr)
  )) %>%
  # $y, rr, take the interpolated values in y from the approx function and assign them to rr
  mutate(rr_interpolated = if_else(step > max(step[!is.na(simulated_rr)]), NA_real_, rr_interpolated)) %>%
  ungroup() %>% 
  select("simulation_id", "disease", "step", "rr_interpolated")



# Random mix
set.seed(123)
rr_table_interpolated <- rr_table_interpolated %>%
  group_by(disease) %>%
  mutate(simulation_id = match(simulation_id, unique(simulation_id)),
         simulation_id = sample(unique(simulation_id))[simulation_id]) %>%
  ungroup()


# Mean and IC
ic95_rr <- rr_table_interpolated %>%
  group_by(disease, step) %>%
  summarise(
    rr_mean = mean(rr_interpolated, na.rm = TRUE),
    rr_lci = quantile(rr_interpolated, 0.025, na.rm = TRUE), # Lower limit of the 95% CI
    rr_uci = quantile(rr_interpolated, 0.975, na.rm = TRUE)  # Upper limit of the 95% CI
  ) %>% 
  ungroup()



################################################################################################################################
#                                                 7. VISUALIZATION OF DRF                                                      #
################################################################################################################################

# Disease name
names_disease <- c(
  "cancer" = "Cancer incidence",
  "cvd" = "Cardiovascular disease incidence",
  "dem" = "Dementia",
  "diab2" = "Type 2 diabetes",
  "dep" = "Depressive symptoms",
  "mort" = "All-cause mortality"
)

# Disease colour
colors_disease <- c(
  "cancer" = "darkorange",      
  "cvd" = "gold" ,
  "dem" = "pink" ,
  "diab2" = "palegreen3",
  "dep" = "slateblue",
  "mort" = "steelblue"
)



# All-cause mortality
  # Central value
  graph_drf_central_mort <- graph_DRF("mort", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_mort)
  
  # Simulated RR normal distributions
  graph_drf_sim_mort <- graph_sim_DRF("mort", rr_table_interpolated)
  plot(graph_drf_sim_mort)

  # RR mean and IC95 points
  graph_drf_mort <- graph_DRF("mort", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_mort)




# Cardiovascular disease incidence
  # Central value
  graph_drf_central_cvd <- graph_DRF("cvd", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_cvd)
  
  # Simulated RR normal distributions
  graph_drf_sim_cvd <- graph_sim_DRF ("cvd", rr_table_interpolated)
  plot(graph_drf_sim_cvd)

  # RR mean and IC95 points
  graph_drf_cvd <- graph_DRF ("cvd", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_cvd)




# Cancer incidence
  # Central value
  graph_drf_central_cancer <- graph_DRF("cancer", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_cancer)
  
  # Simulated RR normal distributions
  graph_drf_sim_cancer <- graph_sim_DRF ("cancer", rr_table_interpolated)
  plot(graph_drf_sim_cancer)

  # RR mean and IC95 points
  graph_drf_cancer <- graph_DRF ("cancer", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_cancer)




# Type 2 diabetes
  # Central value
  graph_drf_central_diab2 <- graph_DRF("diab2", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_diab2)
  
  # Simulated RR normal distributions
  graph_drf_sim_diab2 <- graph_sim_DRF ("diab2", rr_table_interpolated)
  plot(graph_drf_sim_diab2)

  # RR mean and IC95 points
  graph_drf_diab2 <- graph_DRF ("diab2", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_diab2)



# Dementia
  # Central value
  graph_drf_central_dem <- graph_DRF("dem", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_dem)
  
  # Simulated RR normal distributions
  graph_drf_sim_dem <- graph_sim_DRF ("dem", rr_table_interpolated)
  plot(graph_drf_sim_dem)

  # RR mean and IC95 points
  graph_drf_dem <- graph_DRF ("dem", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_dem)



# Depressive symptoms
  # Central value
  graph_drf_central_dep <- graph_DRF("dep", rr_central_table, "mid", "low", "up")
  plot(graph_drf_central_dep)
  
  # Simulated RR normal distributions
  graph_drf_sim_dep <- graph_sim_DRF ("dep", rr_table_interpolated)
  plot(graph_drf_sim_dep)

  # RR mean and IC95 points
  graph_drf_dep <- graph_DRF ("dep", ic95_rr, "rr_mean", "rr_lci", "rr_uci")
  plot(graph_drf_dep)




# All DRF in one plot
  # All DRF central values curves
  list_central_drf <- list(graph_drf_central_mort, graph_drf_central_cvd, graph_drf_central_cancer, graph_drf_central_diab2, graph_drf_central_dem,
                   graph_drf_central_dep)
  
  common_theme <- theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
  
  list_central_drf <- lapply(list_central_drf, function(p) p + common_theme)
  
  combined_plot_central_drf <- reduce(list_central_drf, `+`) + plot_layout(ncol = 3)
  
  print(combined_plot_central_drf) 
  
  
  
  # All DRF curves simulated
  list_drf <- list(graph_drf_sim_mort, graph_drf_sim_cvd, graph_drf_sim_cancer, graph_drf_sim_diab2, graph_drf_sim_dem,
                  graph_drf_sim_dep)
  
  common_theme <- theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
  
  list_drf <- lapply(list_drf, function(p) p + common_theme)
  
  combined_plot_drf <- reduce(list_drf, `+`) + plot_layout(ncol = 3)
  
  print(combined_plot_drf) 


  # Mean + IC95
  list_mean_drf <- list(graph_drf_mort, graph_drf_cvd, graph_drf_cancer, graph_drf_diab2, graph_drf_dem,
                   graph_drf_dep)
  
  common_theme <- theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )
  
  list_mean_drf <- lapply(list_mean_drf, function(p) p + common_theme)
  
  combined_plot_mean_drf <- reduce(list_mean_drf, `+`) + plot_layout(ncol = 3)
  
  print(combined_plot_mean_drf) 




################################################################################################################################
#                                                   8. DATA EXPORTATION                                                        #
################################################################################################################################

# Table of DRF simulated
export(rr_central_table, here("data_clean", "DRF", "rr_central_interpolated.rds"))
export(rr_table_interpolated, here("data_clean", "DRF", "rr_sim_interpolated.rds"))
export(ic95_rr, here("data_clean", "DRF", "rr_interpolated_mean.csv"))



# Graphs (all central DRF curves)
ggsave(here("output", "DRF", "drf_central_mort.png"), plot = graph_drf_central_mort)
ggsave(here("output", "DRF", "drf_central_cvd.png"), plot = graph_drf_central_cvd)
ggsave(here("output", "DRF", "drf_central_cancer.png"), plot = graph_drf_central_cancer)
ggsave(here("output", "DRF", "drf_central_diab2.png"), plot = graph_drf_central_diab2)
ggsave(here("output", "DRF", "drf_central_dem.png"), plot = graph_drf_central_dem)
ggsave(here("output", "DRF", "drf_central_dep.png"), plot = graph_drf_central_dep)

ggsave(here("output", "DRF", "central_drf_all.png"), plot = combined_plot_central_drf)



# Graphs (all DRF curves simulated)
ggsave(here("output", "DRF", "drf_mort.png"), plot = graph_drf_sim_mort)
ggsave(here("output", "DRF", "drf_cvd.png"), plot = graph_drf_sim_cvd)
ggsave(here("output", "DRF", "drf_cancer.png"), plot = graph_drf_sim_cancer)
ggsave(here("output", "DRF", "drf_diab2.png"), plot = graph_drf_sim_diab2)
ggsave(here("output", "DRF", "drf_dem.png"), plot = graph_drf_sim_dem)
ggsave(here("output", "DRF", "drf_dep.png"), plot = graph_drf_sim_dep)

ggsave(here("output", "DRF", "drf_all.png"), plot = combined_plot_drf)



# Graphs (mean + IC95)
ggsave(here("output", "DRF", "drf_mort_mean.png"), plot = graph_drf_mort)
ggsave(here("output", "DRF", "drf_cvd_mean.png"), plot = graph_drf_cvd)
ggsave(here("output", "DRF", "drf_cancer_mean.png"), plot = graph_drf_cancer)
ggsave(here("output", "DRF", "drf_diab2_mean.png"), plot = graph_drf_diab2)
ggsave(here("output", "DRF", "drf_dem_mean.png"), plot = graph_drf_dem)
ggsave(here("output", "DRF", "drf_dep_mean.png"), plot = graph_drf_dep)

ggsave(here("output", "DRF", "drf_all_mean.png"), plot = combined_plot_mean_drf)

