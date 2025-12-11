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
rr_central_table <- import(here("data_clean", "DRF", "rr_central_interpolated.rds"))


# Import functions
source(here("0_Functions.R"))



################################################################################################################################
#                                                      3. PARAMETERS                                                           #
################################################################################################################################

# Import parameters
source(here("0_Parameters.R"))

# Diseases considered
dis_vec = c("mort", "cvd", "cancer", "diab2", "dem", "dep")

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

## -------------------------------------------------------
## REDUCTION RISK FOR ALL DISEASES (EXCEPT BREAST CANCER)
## -------------------------------------------------------

# Baseline step 2000
rr_baseline <- rr_central_table %>%
  filter(step == 2000) %>%
  select(disease, rr2000_mid = mid, rr2000_low = low, rr2000_up = up)

rr_central_table <- rr_central_table %>%
  left_join(rr_baseline, by = "disease")


# Reduction risk for all diseases except breast cancer
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
cases_list <- calc_cases(data_list = walkers_list,
                         rr_table = rr_central_table,
                         week_base = week_base,
                         dis_vec = dis_vec,
                         bound_vec = bound_vec)









################################################################################################################################
#                                              5. TOTAL OF PREVENTED CASES                                                     #
################################################################################################################################


##############################################################
#                        PER DISEASE                         #
##############################################################
# Total of prevented cases
cases <- cases_prevented (data = cases_list,
                               bound_vec = c("low", "mid", "up"),
                               dis_vec = dis_vec,
                               group = NULL)



# Gather results with IC
IC_cases <- cases$mid %>% 
  mutate(tot_cases_low = cases$low[,"tot_cases"], 
         tot_cases_up = cases$up[,"tot_cases"]) %>% 
  select(disease, tot_cases, tot_cases_se, tot_cases_low, tot_cases_up) 



##############################################################
#                      PER SEX AND AGE                       #
##############################################################
# Total of prevented cases per sex and age categories
cases_sex_age <- cases_prevented (data = cases_list,
                               bound_vec = c("low", "mid", "up"),
                               dis_vec = dis_vec,
                               group = c("age_grp10", "sex")) 



# Gather results with IC
IC_cases_sex_age <- cases_sex_age$mid %>% 
  mutate(tot_cases_low = cases_sex_age$low[,"tot_cases"], 
         tot_cases_up = cases_sex_age$up[,"tot_cases"]) %>% 
  select(disease, sex, age_grp10, tot_cases, tot_cases_se, tot_cases_low, tot_cases_up)



##############################################################
#                          PER SEX                           #
##############################################################
# Total of prevented cases per sex
cases_sex <- cases_prevented (data = cases_list,
                                  bound_vec = c("low", "mid", "up"),
                                  dis_vec = dis_vec,
                                  group = "sex")



# Gather results with IC
IC_cases_sex <- cases_sex$mid %>% 
  mutate(tot_cases_low = cases_sex$low[,"tot_cases"], 
         tot_cases_up = cases_sex$up[,"tot_cases"]) %>% 
  select(disease, sex, tot_cases, tot_cases_se, tot_cases_low, tot_cases_up)

  # Re-organize by decreasing order
  IC_tot_cases_sex <- IC_cases_sex %>% 
    left_join(IC_cases %>% select(disease, TOTAL_mixed = tot_cases), by = "disease") %>% 
    arrange(desc(tot_cases)) %>%                      
    mutate(disease = factor(disease, levels = unique(disease))) 




################################################################################################################################
#                                                     6. VISUALIZATION                                                         #
################################################################################################################################

# Plot : Cases prevented by walking in 2019 according to sex 
plot_cases_prev <-
  ggplot(IC_tot_cases_sex, aes(x = disease, y = tot_cases, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  scale_y_break(c(60000, 100000)) +
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x.top = element_blank(),      # supprime labels X en haut
    axis.ticks.x.top = element_blank()      # supprime les ticks X en haut
  )

plot_cases_prev




# Plot : Cases prevented (EXCEPT DEPRESSION)
plot_no_dep_prev <- 
  ggplot(IC_tot_cases_sex, aes(x = disease, y = tot_cases, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
  geom_bar(width = 0.7, position = position_dodge2(.7), stat = "identity")  +
  geom_errorbar(position = position_dodge(.7), width = .25) +
  scale_fill_manual(values = colors_sex) +
  scale_x_discrete(labels = names_disease) + 
  ylab ("Cases prevented") +
  xlab("Disease") +
  theme_minimal()

plot_no_dep_prev


# Plot : DEPRESSION
plot_dep_prev <- IC_tot_cases_sex %>% filter(disease == "dep") %>% 
  ggplot(aes(x = disease, y = tot_cases, ymin = tot_cases_low, ymax = tot_cases_up, fill = sex)) +
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
export(rr_central_table, here("data_clean", "DRF", "reduction_risk_central.xlsx"))
export(IC_cases, here("output", "Tables", "2019", "cases_prev_2019.xlsx"))
export(IC_cases_sex_age, here("output", "Tables", "2019", "cases_prev_2019_sex_age.xlsx"))
export(IC_cases_sex, here("output", "Tables", "2019", "cases_prev_2019_sex.xlsx"))

# Plot 
ggsave(here("output", "Plots", "2019", "cases_prevented.png"), plot = plot_cases_prev)
ggsave(here("output", "Plots", "2019", "dep_cases_prevented.png"), plot = combined_plot_dep)



