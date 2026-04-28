#################################################
############ HEALTH IMPACT ASSESSMENT ###########
#################################################



###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                          HIA - 2019                                                                                     #
###########################################################################################################################################################################
###########################################################################################################################################################################


################################################################################################################################
#                                                    1. LOAD PACKAGES                                                          #
################################################################################################################################
pacman :: p_load(
  rio,          # Data importation
  here,         # Localization of files 
  dplyr,        # Data management
  tidyr,        # Pivot tables
  srvyr,        # Survey
  survey,
  ggplot2,      # Data visualization
  patchwork,    # Plot combination
  ggbreak       # Break 
)




################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################
# Walkers dataset
emp_walk <- import(here("data_clean", "EMP_dis_walkers.xlsx"))


# RR by step, simulated dose-response relationships
rr_central_table <- import(here("data_clean", "Diseases", "DRF", "rr_central_interpolated.rds"))


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
#                                                    4. DATA PREPARATION                                                       #
################################################################################################################################

##############################################################
#                     WALKERS DATASET                        #
##############################################################

# Initialization
emp_walk <- emp_walk %>% 
  # Round the number of steps to the nearest hundred and baseline at 2000
  mutate(step = pmin(12000, round(step_commute / 10) * 10 + 2000))



# EMP Dataset per disease and bound
walkers_list <- list()

for(bound in bound_vec) {
  bound_list <- list()
  
  for(dis in dis_vec) {
    dis_walkers <- emp_walk %>% filter(disease == dis)
    
    bound_list[[dis]] <- dis_walkers
  }
  
  walkers_list[[bound]] <- bound_list
}





##############################################################
#                DISEASE REDUCTION RISK                      #
##############################################################

## ---------------------------------
## REDUCTION RISK FOR ALL DISEASES 
## ---------------------------------

# Baseline step 2000
rr_baseline <- rr_central_table %>%
  filter(step == 2000) %>%
  select(disease, rr2000_mid = mid, rr2000_low = low, rr2000_up = up)

rr_central_table <- rr_central_table %>%
  left_join(rr_baseline, by = "disease")


# Reduction risk for all diseases
for (bound in bound_vec) {
  rr_central_table <- rr_central_table %>%
    mutate(!!paste0("reduction_risk_", bound) := 
             (.data[[paste0("rr2000_", bound)]] - .data[[bound]]) /.data[[paste0("rr2000_", bound)]]) 
}
# Rename column: To calculate the upper bound of reduction of the relative risk, use RR lower bound because the decrease will be higher i.e. the person exposed (walking) is less likely to have the disease 
rr_central_table <- rr_central_table %>% 
  rename(reduction_risk_low = reduction_risk_up,
         reduction_risk_up = reduction_risk_low)



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################
# Calculate for each individual the number of prevented cases, DALY and costs
HIA_list <- calc_HIA(data_list = walkers_list,
                       rr_table = rr_central_table,
                       dw_table = dw_table,
                       week_base = week_base,
                       dis_vec = dis_vec,
                       bound_vec = bound_vec)





################################################################################################################################
#                                   5. HIA OUTCOMES : Total of prevented cases, DALY, costs                                    #
################################################################################################################################

##############################################################
#                        PER DISEASE                         #
##############################################################
# Total of prevented burden of each disease
burden <- burden_prevented(data_list = HIA_list, 
                           dis_vec = dis_vec,
                           bound_vec,
                           group = NULL)



##############################################################
#                      PER SEX AND AGE                       #
##############################################################
# Total of prevented burden per sex and age categories
burden_sex_age <- burden_prevented(data_list = HIA_list, 
                                   dis_vec = dis_vec,
                                   bound_vec,
                                   group = c("age_grp10", "sex"))



##############################################################
#                          PER SEX                           #
##############################################################
# Total of prevented burden per sex
burden_sex <- burden_prevented(data_list = HIA_list, 
                                   dis_vec = dis_vec,
                                   bound_vec,
                                   group = "sex")




  # Re-organize by decreasing order
  burden_sex_order <- burden_sex %>% 
    left_join(burden %>% select(disease, TOTAL_mixed = tot_cases_mid), by = "disease") %>% 
    arrange(desc(tot_cases_mid)) %>%                      
    mutate(disease = factor(disease, levels = unique(disease)))



################################################################################################################################
#                                                     6. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Cases prevented by walking in 2019 according to sex 
plot_cases_prev <- burden_sex_order %>% filter (disease != "bc") %>% 
  ggplot(aes(x = disease, y = tot_cases_mid, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x.top = element_blank(),      # delete labels X at the top
    axis.ticks.x.top = element_blank()      # delete ticks X at the top
  )

plot_cases_prev




# Plot : Cases prevented (EXCEPT DEPRESSION)
plot_no_dep_prev <- burden_sex_order %>%  filter(!disease %in% c("bc", "dep")) %>% 
  ggplot(aes(x = disease, y = tot_cases_mid, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal()

plot_no_dep_prev


# Plot : DEPRESSION
plot_dep_prev <- burden_sex_order %>% filter(disease == "dep") %>% 
  ggplot(aes(x = disease, y = tot_cases_mid, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal()

plot_dep_prev


# Plot : Combine plot 

combined_plot_dep <- plot_no_dep_prev + plot_dep_prev + plot_layout(ncol = 2)

combined_plot_dep





################################################################################################################################
#                                                      7. DESCRIPTION                                                         #
################################################################################################################################




################################################################################################################################
#                                                      8. EXPORT DATA                                                          #
################################################################################################################################
# Tables
export(rr_central_table, here("data_clean", "Diseases", "DRF", "reduction_risk_central.xlsx"))
export(burden, here("output", "Tables", "2019", "cases_prev_2019.xlsx"))
export(burden_sex_age, here("output", "Tables", "2019", "cases_prev_2019_sex_age.xlsx"))
export(burden_sex, here("output", "Tables", "2019", "cases_prev_2019_sex.xlsx"))

# Plot 
ggsave(here("output", "Plots", "2019", "cases_prevented.png"), plot = plot_cases_prev)
ggsave(here("output", "Plots", "2019", "dep_cases_prevented.png"), plot = combined_plot_dep)
