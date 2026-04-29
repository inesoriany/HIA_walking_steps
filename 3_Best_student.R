#################################################
#############     BEST STUDENTS     #############
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
  ggplot2       # Plotting
)


################################################################################################################################
#                                                    2. IMPORT DATA                                                            #
################################################################################################################################

# Walking trips dataset
emp_walk_trip <- import(here("data_clean", "EMP_dis_walking_trips.xlsx"))


# Risk reductions
reduc_rr_table <- import(here("data_clean", "Diseases", "DRF", "reduction_risk_central.xlsx"))


# Disability weights
dw_table <- import(here("data", "dw_table.xlsx"))



# Import functions
source(here("0_Functions.R"))


################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "bc", "cvd", "cancer", "diab2", "dem", "dep")

# Bound
bound_vec <- c("mid", "low", "up")



################################################################################################################################
#                                                4. MAIN WALKING TRIPS DATASET                                                 #
################################################################################################################################
# Only exclusively walking trips (no intermodal walk)
  emp_main_walk_trip <- emp_walk_trip %>% 
    filter(intermodal_walk_time == 0)

# Add population by age group and area type
emp_main_walk_trip <- emp_main_walk_trip %>%
  left_join(emp_main_walk_trip %>%
      select(ident_ind, age_grp10, area_type, pond_indc) %>%
      distinct() %>%
      group_by(age_grp10, area_type) %>%
      summarise(pop_age_area = sum(pond_indc, na.rm = TRUE), .groups = "drop"), 
      by = c("age_grp10", "area_type"))


# Sum steps by individual, area type, and disease
steps_by_individual_area_disease <- emp_main_walk_trip %>% 
  group_by(ident_ind, area_type, disease) %>% 
  summarise(
    total_steps = sum(step_commute, na.rm = TRUE),
    n_trips = n_distinct(ident_dep),
    pond_indc = first(pond_indc),
    .groups = "drop"
  )

# Merge aggregated data back to original trips dataset
emp_main_walk_trip <- emp_main_walk_trip %>% 
  left_join(
    steps_by_individual_area_disease %>% 
      select(ident_ind, area_type, disease, total_steps, n_trips),
    by = c("ident_ind", "area_type", "disease")
  )


################################################################################################################################
#                                                   5. AGE DISTRIBUTION                                                        #
################################################################################################################################
# Survey design ponderated by day
main_jour <- emp_main_walk_trip %>% 
  filter(pond_jour != "NA") %>% 
  as_survey_design(ids = ident_ind,
                   weights = pond_jour,
                   strata = c(sex, age_grp10),
                   nest = TRUE)



# Walking pyramid : Age distribution of walking volume for each territory (in steps)
distrib_main_walk_EMP2019 <- main_jour %>% 
  group_by(age_grp10, area_type) %>% 
  summarise(mean_ind = survey_mean(step_commute, na.rm = TRUE))


# Distribution coefficient by area type
distrib_main_walk_EMP2019 <- distrib_main_walk_EMP2019 %>% 
  group_by(area_type) %>% 
  mutate(rho = mean_ind / mean_ind[age_grp10 == "20-29"]) %>% 
  ungroup()

emp_main_walk_trip <- emp_main_walk_trip %>% 
  left_join(distrib_main_walk_EMP2019 %>% 
  select(age_grp10, rho), by = "age_grp10")


# Calculate target to reach for each individual depending on their area type and age group
#(steps per person by age group and area type)
emp_main_walk_trip <- emp_main_walk_trip %>% 
  group_by(area_type) %>%
  mutate(
    target_area = case_when(
      area_type == "urban" ~ urban_target,
      area_type == "periurban" ~ periurban_target,
      area_type == "rural" ~ rural_target),
    
    # structural weights
    weight = rho * pop_age_area,
    
    # average distribution coefficient for the area type
    rho_bar = sum(weight, na.rm = TRUE) / sum(pop_age_area, na.rm = TRUE),
    
    # target steps per person for each individual, adjusted by the distribution coefficient and population structure
    target_pp = target_area * (rho / rho_bar)) %>%

  ungroup()


################################################################################################################################
#                                                   5. DATA PREPARATION                                                        #
################################################################################################################################

# Filter individuals below targets
walk_below_targets <- emp_main_walk_trip  %>% 
  filter(total_steps < target_pp)


# Initialization
walk_below_targets <- walk_below_targets %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000, round((target_pp - total_steps)/ 10) * 10 + 2000))



# EMP Dataset per disease and bound
walk_trips_target_list <- list()

for(bound in bound_vec) {
  bound_list <- list()
  
  for(dis in dis_vec) {
    dis_walkers <- walk_below_targets %>% 
    filter(disease == dis)
    
    bound_list[[dis]] <- dis_walkers
  }
  
  walk_trips_target_list[[bound]] <- bound_list
}



################################################################################################################################
#                                                6. HEALTH IMPACT ASSESSMENT                                                   #
################################################################################################################################
# Calculate for each individual the number of prevented cases, DALY and costs
HIA_target_list <- calc_HIA(data_list = walk_trips_target_list,
                       rr_table = reduc_rr_table,
                       dw_table = dw_table,
                       dis_vec = dis_vec,
                       bound_vec = bound_vec)



################################################################################################################################
#                                   7. HIA OUTCOMES : Total of prevented cases, DALY, costs                                    #
################################################################################################################################

##############################################################
#                        PER DISEASE                         #
##############################################################
# Total of prevented burden of each disease
burden_target <- burden_prevented(data_list = HIA_target_list, 
                                  dis_vec = dis_vec,
                                  bound_vec = bound_vec,
                                  group = NULL)


##############################################################
#                       PER AREA TYPE                        #
##############################################################
# Prevented diseases per area type
burden_target_area <- burden_prevented(data_list = HIA_target_list, 
                                       dis_vec = dis_vec,
                                       bound_vec = bound_vec,
                                       group = "area_type")


# Re-organize by decreasing order
burden_target_area_order <- burden_target_area %>% 
    left_join(burden_target %>% select(disease, TOTAL_mixed = tot_cases_mid), by = "disease") %>% 
    arrange(desc(tot_cases_mid)) %>%                      
    mutate(disease = factor(disease, levels = unique(disease)))




##############################################################
#                 PER SEX, AGE AND AREA TYPE                 #
##############################################################
# Prevented diseases per sex and age categories
burden_target_area_sex_age <- burden_prevented(data_list = HIA_target_list, 
                                               dis_vec = dis_vec,
                                               bound_vec = bound_vec,
                                               group = c("age_grp10", "sex", "area_type"))



##############################################################
#                    PER AGE AND AREA TYPE                   #
##############################################################
# Prevented diseases per age groups and area type
burden_target_area_age <- burden_prevented(data_list = HIA_target_list, 
                                               dis_vec = dis_vec,
                                               bound_vec = bound_vec,
                                               group = c("age_grp10", "area_type"))



##############################################################
#                    PER SEX AND AREA TYPE                   #
##############################################################
# Prevented diseases per sex and area type
burden_target_area_sex <- burden_prevented(data_list = HIA_target_list, 
                                           dis_vec = dis_vec,
                                           bound_vec = bound_vec,
                                           group = c("sex", "area_type"))



################################################################################################################################
#                                                     7. VISUALIZATION                                                         #
################################################################################################################################

##############################################################
#                       PER AREA TYPE                        #
##############################################################
# Plot : Total cases prevented per area type
plot_cases_prev_by_area <- burden_target_area_order %>% filter (disease != "bc")  %>% 
  ggplot(aes(x = area_type, y = tot_cases_mid, fill = area_type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("urban" = "#1f78b4", "periurban" = "#33a02c", "rural" = "#e31a1c")) +
  labs(
    title = "Total cases prevented by walking - best student scenario",
    x = "Territory type",
    y = "Cases prevented (mid estimate)",
    fill = "Territory type"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

plot_cases_prev_by_area



# Plot : Cases prevented by disease and area type
plot_cases_prev_by_area_disease <- burden_target_area_order  %>% filter (disease != "bc")  %>% 
  ggplot(aes(x = disease, y = tot_cases_mid, fill = area_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("urban" = "#1f78b4", "periurban" = "#33a02c", "rural" = "#e31a1c")) +
  scale_x_discrete(labels = names_disease) +
  labs(
    title = "Cases prevented by disease and territory type",
    x = "Disease",
    y = "Cases prevented (mid estimate)",
    fill = "Territory type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

plot_cases_prev_by_area_disease


##############################################################
#                    PER AGE AND AREA TYPE                   #
##############################################################
# Plot : Cases prevented by age groups, separated by area type
plot_cases_prev_by_age_area <- burden_target_area_age %>% filter (disease != "bc")  %>% 
  filter(disease != "bc") %>% 
  ggplot(aes(x = age_grp10, y = tot_cases_mid, fill = age_grp10)) +
  geom_col() +
  facet_wrap(~ area_type, scales = "free_x") +
  scale_fill_brewer(palette = "Spectral") +
  labs(
    title = "Cases prevented by walking - best student scenario",
    subtitle = "By age group and territory type",
    x = "Age group",
    y = "Cases prevented (mid estimate)",
    fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"))

plot_cases_prev_by_age_area



# Plot : Cases prevented by age group and disease, separated by area type
plot_cases_prev_by_age_area_disease <- burden_target_area_age %>% 
  filter(disease != "bc") %>% 
  ggplot(aes(x = age_grp10, y = tot_cases_mid, fill = disease)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  facet_wrap(~ area_type, scales = "free_x") +
  scale_fill_manual(values = colors_disease) +
  labs(
    title = "Cases prevented by age group, disease and territory type",
    x = "Age group",
    y = "Cases prevented (mid estimate)",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

plot_cases_prev_by_age_area_disease



##############################################################
#                    PER SEX AND AREA TYPE                   #
##############################################################
# Plot : Cases prevented by sex, separated by area type
plot_cases_prev_by_sex_area <- burden_target_area_sex %>% filter (disease != "bc")  %>% 
  filter(disease != "bc") %>% 
  ggplot(aes(x = sex, y = tot_cases_mid, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ area_type) +
  scale_fill_manual(values = colors_sex) +
  labs(
    title = "Cases prevented by walking - best student scenario",
    subtitle = "By sex and territory type",
    x = "Sex",
    y = "Cases prevented (mid estimate)",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

plot_cases_prev_by_sex_area



# Plot : Cases prevented by sex and disease, separated by area type
plot_cases_prev_by_sex_area_disease <- burden_target_area_sex %>% 
  filter(disease != "bc") %>% 
  ggplot(aes(x = sex, y = tot_cases_mid, fill = disease)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  facet_wrap(~ area_type) +
  scale_fill_manual(values = colors_disease) +
  labs(
    title = "Cases prevented by sex, disease and territory type",
    x = "Sex",
    y = "Cases prevented (mid estimate)",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

plot_cases_prev_by_sex_area_disease




################################################################################################################################
#                                                      8. DESCRIPTION                                                          #
################################################################################################################################
# Proportion

# Description : how many people below the tragets and age/sex distribution below targets gy area type
# Filter individuals below targets ? and do HIA (to know the added value)
  ################## See how to integrate that later ----
  # Check if individuals meet the targets
  below_targets <- emp_main_walk_trip %>% 
    mutate(
      target = case_when(
        area_type == "urban" ~ urban_target,
        area_type == "periurban" ~ periurban_target,
        area_type == "rural" ~ rural_target,
        TRUE ~ NA_real_
      ),
      meets_target = total_steps >= target
    ) %>% 
    filter(!meets_target)

  cat("\n=== PERSONNES EN DESSOUS DES CIBLES ===\n")
  cat("Total (pondéré) :", sum(below_targets$pond_indc, na.rm = TRUE), "personnes\n")
  cat("Total (non pondéré) :", nrow(below_targets), "personnes\n\n")

  by_area <- below_targets %>% 
    group_by(area_type) %>% 
    summarise(
      n_below_weighted = sum(pond_indc, na.rm = TRUE),
      n_below_unweighted = n(),
      mean_steps = mean(total_steps, na.rm = TRUE),
      .groups = "drop"
    )
  print(by_area)


# same for above targets
# Filter individuals above targets
walk_above_targets <- emp_main_walk_trip  %>% 
  filter(total_steps > target_pp)




################################################################################################################################
#                                                      9. EXPORT DATA                                                          #
################################################################################################################################




