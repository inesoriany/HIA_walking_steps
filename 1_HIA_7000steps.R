#################################################
############ HEALTH IMPACT ASSESSMENT ###########
#################################################


###########################################################################################################################################################################
###########################################################################################################################################################################
#                                                                       HIA - 7 000 STEPS                                                                                 #
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
  ggplot2       # Data visualization
)


################################################################################################################################
#                                                     2. IMPORT DATA                                                           #
################################################################################################################################

# Walkers dataset
emp_step <- import(here("data_clean", "EMP_dis_walkers.xlsx"))

# Risk reductions
reduc_rr_table <- import(here("data_clean", "DRF", "reduction_risk_central.xlsx"))

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

# Initialization
emp_step <- emp_step %>% 
  # Recommendation of 7000 steps and baseline at 2000
  mutate(step = 7000 + baseline_step)


# EMP Dataset per disease
RECO_walkers_list <- list()

for(bound in bound_vec) {
  bound_list <- list()
  
  for(dis in dis_vec) {
    dis_walkers <- emp_step %>%
      filter(disease == dis)
    
    bound_list[[dis]] <- dis_walkers
  }
  
  RECO_walkers_list[[bound]] <- bound_list
}



################################################################################################################################
#                                                 4. HEALTH IMPACT ASSESSMENT                                                  #
################################################################################################################################

RECO_HIA_list <- calc_HIA(data_list = RECO_walkers_list,
                          rr_table = reduc_rr_table,
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
RECO_burden <- burden_prevented(data_list = RECO_HIA_list, 
                                dis_vec = dis_vec,
                                bound_vec,
                                group = NULL)


##############################################################
#                          PER SEX                           #
##############################################################
RECO_burden_sex <- burden_prevented(data_list = RECO_HIA_list, 
                                    dis_vec = dis_vec,
                                    bound_vec,
                                    group = "sex")





## -------------------------------------------------------
## DECREASING ORDER
## -------------------------------------------------------
RECO_burden_sex_order <- RECO_burden_sex %>% 
  left_join(RECO_burden %>% select(disease, TOTAL_mixed = tot_cases_mid), by = "disease") %>% 
  arrange(desc(tot_cases_mid)) %>%                      
  mutate(disease = factor(disease, levels = unique(disease))) 



##############################################################
#                          GLOBAL                            #
##############################################################
# Total of prevented cases
RECO_burden <- burden_prevented(data_list = RECO_HIA_list, 
                                    dis_vec = dis_vec,
                                    bound_vec,
                                    group = NULL)



## -------------------------------------------------------
## ADDITIONNAL GAINS
## -------------------------------------------------------
# Import 2019 burden prevented
burden_2019 <- import(here("output", "Tables", "2019", "cases_prev_2019.xlsx"))

# Data preparation
burden_2019_row <- burden_2019 %>% 
  rename_with(.fn = ~ paste0(.x, "_2019"),
              .cols = -c (disease))

RECO_burden_row <- RECO_burden %>% 
  rename_with(.fn = ~ paste0(.x, "_RECO"),
              .cols = -c (disease))


# Additional prevented cases for each disease according to sex
add_RECO_burden <- burden_2019_row %>%
  left_join(RECO_burden_row, by = "disease", suffix = c("_2019", "_RECO")) %>%
  mutate(across(
    ends_with("_RECO"),
    ~ . - get(sub("_RECO$", "_2019", cur_column())),
    .names = "{.col}_diff"
  ))



## -------------------------------------------------------
## MORBIDITY
## -------------------------------------------------------
RECO_burden %>% 
  filter(disease != "mort") %>% 
  summarise(tot_cases_mid = sum(tot_cases_mid, na.rm = TRUE),
            tot_cases_low = sum(tot_cases_low, na.rm = TRUE),
            tot_cases_up = sum(tot_cases_up, na.rm = TRUE))



## -------------------------------------------------------
## ALL DISEASES
## -------------------------------------------------------
RECO_burden %>% 
  summarise(tot_cases_mid = sum(tot_cases_mid, na.rm = TRUE),
            tot_cases_low = sum(tot_cases_low, na.rm = TRUE),
            tot_cases_up = sum(tot_cases_up, na.rm = TRUE))



################################################################################################################################
#                                                     6. VISUALIZATION                                                         #
################################################################################################################################
# Import 2019 data
burden_sex_2019 <- import(here("output", "Tables", "2019", "cases_prev_2019_sex.xlsx"))

# Plot : Cases prevented had the 2019 population walked 7000 steps according to sex, compared to 2019 levels
plot_RECO_cases_prev <- 
  ggplot() +
  geom_bar(data = burden_sex_2019 %>%  filter(disease != "bc"),
           mapping = aes(x = disease, y = tot_cases_mid, fill = sex, alpha = "2019 baseline"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  
  geom_errorbar(data = burden_sex_2019 %>%  filter(disease != "bc"),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "2019 baseline"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_fill_manual(values = colors_sex) +
  
  
  geom_bar(data = RECO_burden_sex_order %>%  filter(disease != "bc"), 
           mapping = aes(x = disease, y = tot_cases_mid, fill = sex, alpha = "7000 steps recommendation"),
           width = 0.7,
           position = position_dodge2(0.7),
           stat = "identity") +
  scale_alpha_manual(name   = "Scenario",
                     values = c("2019 baseline" = 1, "7000 steps recommendation" = 0.4)) +
  
  geom_errorbar(data = RECO_burden_sex_order %>%  filter(disease != "bc"),
                mapping = aes(x = disease, ymin = tot_cases_low, ymax = tot_cases_up, group = sex, alpha = "7000 steps recommendation"),
                position = position_dodge(0.7),
                width = 0.25) +
  
  scale_x_discrete(labels = names_disease) + 
  ylab("Cases prevented") +
  xlab("Disease") +
  theme_minimal() 

plot_RECO_cases_prev 






################################################################################################################################
#                                                      8. EXPORT DATA                                                          #
################################################################################################################################
# Tables
export(RECO_burden_sex, here("output", "Tables", "7000 steps", "cases_prev_7000steps.xlsx"))
export(add_RECO_burden_sex, here("output", "Tables", "7000 steps", "add_cases_prev_7000steps.xlsx"))

# Plot 
ggsave(here("output", "Plots", "7000 steps", "cases_prev_7000steps.png"), plot = plot_RECO_cases_prev)










